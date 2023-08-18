#include "volume.h"

#include "argparse.h"

#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <PltFileManager.H>

#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/MultiResGrid.h>


#include <array>
#include <span>


using FloatMultiGrid    = openvdb::tools::MultiResGrid<openvdb::FloatTree>;
using FloatMultiGridPtr = openvdb::tools::MultiResGrid<openvdb::FloatTree>::Ptr;
using FloatGrid         = openvdb::FloatGrid;
using FloatGridPtr      = openvdb::FloatGrid::Ptr;

struct VolumeConfig {
    int                   max_level = -1;
    std::set<std::string> variables;
};

struct SampledGrid {
    std::string       name;
    FloatMultiGridPtr grid;
};

struct Result {
    openvdb::BBoxd           bbox;
    std::vector<SampledGrid> grids;
};

struct LPlt : public pele::physics::pltfilemanager::PltFileManager {
public:
    LPlt(std::filesystem::path path)
        : pele::physics::pltfilemanager::PltFileManager(path) { }

    ~LPlt() { std::cout << "Destroying pltfile state" << std::endl; }

    LPlt(LPlt const&)                  = delete;
    LPlt const& operator=(LPlt const&) = delete;
    LPlt(LPlt&&)                       = delete;
    LPlt const& operator=(LPlt&&)      = delete;

    auto const& all_data() const { return m_data; }

    auto const& get_data(int nlvl) { return m_data.at(nlvl); }

    auto const& get_dmap(int nlvl) { return m_dmaps.at(nlvl); }
};

void loop_box(amrex::Box const&                       bx,
              amrex::Array4<amrex::Real const> const& a,
              std::span<openvdb::FloatGrid::Accessor> accessors,
              std::span<int const>                    var_indices) {

    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    // std::cout << "BOX " << lo << hi << std::endl;

    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {

                for (int offset = 0; offset < var_indices.size(); offset++) {
                    auto index = var_indices[offset];

                    float value = *a.ptr(i, j, k, index);

                    accessors[offset].setValue({ i, j, k }, value);

                    assert(accessors[offset].getValue({ i, j, k }) == value);
                }
            }
        }
    }
}

struct AMRState {
    VolumeConfig  config;
    amrex::AMReX* amr;

    std::unique_ptr<LPlt> plt_data;

    std::vector<std::string> var_names;
    std::vector<int>         var_ids;

    AMRState(VolumeConfig c) : config(c) {
        // dummy params
        std::string          exe_name = "amr2vdb";
        std::array<char*, 1> argv_src = { exe_name.data() };
        int                  argc     = argv_src.size();
        char**               argv     = argv_src.data();

        amr = amrex::Initialize(argc, argv);
    }

    AMRState(AMRState const&)                  = delete;
    AMRState const& operator=(AMRState const&) = delete;
    AMRState(AMRState&&)                       = delete;
    AMRState const& operator=(AMRState&&)      = delete;

    bool init(std::filesystem::path path) {
        if (!std::filesystem::exists(path)) {
            std::cerr << "Path does not exist: " << path << std::endl;
            return false;
        }

        plt_data      = std::make_unique<LPlt>(path);
        auto plt_vars = plt_data->getVariableList();
        // int  nvars    = plt_vars.size();

        {
            std::set<std::string> extant_names(plt_vars.begin(),
                                               plt_vars.end());

            for (auto const& vname : config.variables) {
                if (!extant_names.contains(vname)) {
                    std::cerr << "Variable " << vname
                              << " does not exist in pltfile.\n";
                    return false;
                }
            }
        }

        std::cout << "Loading from " << path << std::endl;

        for (long i = 0; i < plt_vars.size(); i++) {
            auto const& vname = plt_vars[i];
            if (config.variables.contains(vname)) {
                std::cout << "Variable " << vname << " at " << i << std::endl;
                var_names.push_back(vname);
                var_ids.push_back(i);
            }
        }

        if (config.max_level < 0) {
            config.max_level = plt_data->getNlev() - 1;
        }

        if (config.max_level >= plt_data->getNlev()) {
            std::cerr << "Asking for refinement level " << config.max_level
                      << " but pltfile only has " << plt_data->getNlev()
                      << std::endl;
            return false;
        }

        plt_data->readPlotFileData();

        std::cout << "Done." << std::endl;

        return true;
    }

    Result write_to_vdbs() const {
        auto geom = plt_data->getGeom(config.max_level);

        auto domain = geom.Domain();
        auto larray = domain.smallEnd().toArray();
        auto harray = domain.bigEnd().toArray();

        auto grid = amrex::BoxArray(domain);

        std::cout << "Level:    " << config.max_level << std::endl;
        std::cout << "Domain L: " << larray[0] << " " << larray[1] << " "
                  << larray[2] << std::endl;
        std::cout << "Domain H: " << harray[0] << " " << harray[1] << " "
                  << harray[2] << std::endl;

        std::cout << "Vars:     " << plt_data->getVariableList().size()
                  << std::endl;

        std::vector<FloatMultiGrid::Ptr> vdb_grids;

        for (size_t i = 0; i < var_ids.size(); i++) {
            auto g =
                std::make_shared<FloatMultiGrid>(config.max_level + 1, 0.0f);
            vdb_grids.push_back(g);
        }


        // note the <= here
        for (int current_level = 0; current_level <= config.max_level;
             current_level++) {
            std::cout << "Working on level " << current_level << std::endl;
            // also note that VDB does inverse sampling. 0 is the finest.
            // so we want to map, say 3 -> 0 and 0 -> 3

            size_t vdb_mapped_level = (config.max_level - current_level);

            std::cout << "Mapping to VDB level " << vdb_mapped_level
                      << std::endl;

            // now check if we have a coarser grid to copy into this one

            //            if (vdb_mapped_level != config.max_level) {
            //                std::cout << "Interpolation..." << std::endl;
            //                for (auto& mg : vdb_grids) {
            //                    openvdb::tools::resampleToMatch<openvdb::tools::BoxSampler>(
            //                        *(mg->grid(vdb_mapped_level + 1)),
            //                        *(mg->grid(vdb_mapped_level)));
            //                }
            //            }


            // for each level, get each vars grid
            std::vector<FloatGridPtr>        per_var_grid;
            std::vector<FloatGrid::Accessor> per_var_accessor;

            per_var_grid.reserve(vdb_grids.size());
            per_var_accessor.reserve(vdb_grids.size());

            for (auto& mg : vdb_grids) {
                auto grid = mg->grid(vdb_mapped_level);

                assert(!!grid);

                per_var_grid.push_back(grid);
                per_var_accessor.push_back(grid->getAccessor());
            }

            auto const& ba = plt_data->getGrid(current_level);

            auto const& mf = plt_data->get_data(current_level);

            auto iter = amrex::MFIter(mf);

            std::cout << "Iterating blocks..." << std::endl;

            for (; iter.isValid(); ++iter) {
                auto const& box = iter.validbox();

                amrex::FArrayBox const& fab = mf[iter];

                auto const& a = fab.array();

                loop_box(box, a, per_var_accessor, var_ids);
            }
        }

        std::cout << "Sampling complete" << std::endl;

        std::vector<SampledGrid> ret_grid;

        for (int i = 0; i < vdb_grids.size(); i++) {
            ret_grid.emplace_back(SampledGrid { var_names[i], vdb_grids[i] });
        }

        Result ret;

        ret.grids = std::move(ret_grid);
        ret.bbox  = {
            { (double)larray[0], (double)larray[1], (double)larray[2] },
            { (double)harray[0], (double)harray[1], (double)harray[2] }
        };


        return ret;
    }

    ~AMRState() { amrex::Finalize(amr); }
};


Result load_file(std::filesystem::path path, VolumeConfig const& c) {
    auto state = AMRState(c);

    if (!state.init(path)) { return {}; }

    return state.write_to_vdbs();
}

FloatGridPtr resample(SampledGrid source, openvdb::BBoxd box) {
    // create new grid
    auto new_grid =
        std::make_shared<FloatMultiGrid>(source.grid->numLevels(), 0.0f);
    new_grid->setName(source.name);

    // interpolate from coarse to fine to this grid

    int max_levels = source.grid->numLevels();
    int level      = max_levels;

    while (level-- > 0) {

        if (level != max_levels - 1) {
            std::cout << "Interpolation " << level << std::endl;
            openvdb::tools::resampleToMatch<openvdb::tools::QuadraticSampler>(
                *(new_grid->grid(level + 1)), *(new_grid->grid(level)));
        }

        // copy over new data
        openvdb::tools::compReplace(*new_grid->grid(level),
                                    *source.grid->grid(level));
    }

    auto ret_grid = new_grid->grid(0);

    std::cout << "Final " << ret_grid->activeVoxelCount() << std::endl;

    ret_grid->clipGrid(box);
    ret_grid->setName(source.name);
    return ret_grid;
}

int amr_to_volume(Arguments const& c) {
    if (c.positional.size() < 3) return EXIT_FAILURE;

    std::string  source_plt = c.positional.at(0);
    std::string  dest_vdb   = c.positional.at(1);
    VolumeConfig config;


    for (auto iter = c.positional.begin() + 2; iter != c.positional.end();
         ++iter) {
        config.variables.insert(*iter);
    }

    auto amr = load_file(source_plt, config);


    openvdb::GridPtrVec upsampled;

    for (auto& multires : amr.grids) {
        auto new_grid = resample(multires, amr.bbox);
        multires.grid.reset(); // try to minimize mem usage
        upsampled.push_back(new_grid);
    }

    openvdb::io::File file(dest_vdb);
    file.write(upsampled);
    file.close();

    return EXIT_SUCCESS;
}
