#include "amr.h"

#include "config.h"

#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <PltFileManager.H>

#include <openvdb/openvdb.h>
#include <openvdb/tools/MultiResGrid.h>

#include <array>
#include <span>

struct LPlt : public pele::physics::pltfilemanager::PltFileManager {
public:
    using pele::physics::pltfilemanager::PltFileManager::PltFileManager;

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
    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {

                for (auto index : var_indices) {
                    accessors[index].setValue({ i, j, k },
                                              a.ptr(i, j, k)[index]);
                }
            }
        }
    }
}

struct AMRState {
    Config        config;
    amrex::AMReX* amr;

    std::unique_ptr<LPlt> plt_data;

    std::vector<std::string> var_names;
    std::vector<int>         var_ids;

    AMRState(Config c) : config(c) {
        // dummy params
        std::string          exe_name = "amr2vdb";
        std::array<char*, 1> argv_src = { exe_name.data() };
        int                  argc     = argv_src.size();
        char**               argv     = argv_src.data();

        amr = amrex::Initialize(argc, argv);
    }

    bool init(std::filesystem::path path) {
        std::cout << path << std::endl;
        if (!std::filesystem::exists(path)) {
            std::cerr << "Path does not exist: " << path << std::endl;
            return false;
        }

        plt_data      = std::make_unique<LPlt>(path);
        auto plt_vars = plt_data->getVariableList();
        int  nvars    = plt_vars.size();

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

        for (long i = 0; i < plt_vars.size(); i++) {
            auto const& vname = plt_vars[i];
            if (config.variables.contains(vname)) {
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

        return true;
    }

    std::vector<SampledGrid> write_to_vdbs() const {
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

        // Index by variable id
        std::vector<FloatMultiGrid::Ptr> vdb_grids;

        for (size_t i = 0; i < var_ids.size(); i++) {
            auto g = std::make_shared<FloatMultiGrid>(4, 0.0f);
            vdb_grids.push_back(g);
        }


        // note the <= here
        for (int current_level = 0; current_level <= config.max_level;
             current_level++) {
            std::cout << "Working on level " << current_level << std::endl;
            // for each level, get each vars grid
            std::vector<FloatGridPtr>        per_var_grid;
            std::vector<FloatGrid::Accessor> per_var_accessor;

            // also note that VDB does inverse sampling. 0 is the finest.
            // so we want to map, say 3 -> 0 and 0 -> 3

            size_t vdb_res = (config.max_level - current_level);

            std::cout << "Mapping to VDB level " << vdb_res << std::endl;

            for (auto& mg : vdb_grids) {
                // MAGIC 1 here is for trilinear sampling

                per_var_grid.push_back(mg->createGrid<1>(vdb_res));
                per_var_accessor.push_back(per_var_grid.back()->getAccessor());
            }

            auto const& ba = plt_data->getGrid(current_level);

            auto const& mf = plt_data->get_data(current_level);

            std::cout << "Count of mf: " << mf.size() << std::endl;

            auto iter = amrex::MFIter(mf);

            for (; iter.isValid(); ++iter) {
                auto const& box = iter.validbox();

                amrex::FArrayBox const& fab = mf[iter];

                auto const& a = fab.array();

                loop_box(box, a, per_var_accessor, var_ids);
            }
        }

        std::cout << "Sampling complete" << std::endl;

        std::vector<SampledGrid> ret;

        for (int i = 0; i < vdb_grids.size(); i++) {
            ret.emplace_back(SampledGrid { var_names[i], vdb_grids[i] });
        }
        return ret;
    }

    ~AMRState() { amrex::Finalize(amr); }
};

std::shared_ptr<AMRState> load_file(std::filesystem::path path,
                                    Config const&         c) {
    auto ret = std::make_shared<AMRState>(c);

    if (ret->init(path)) { return ret; }

    return {};
}


std::vector<SampledGrid> write_to_vdbs(AMRState const& state) {
    return state.write_to_vdbs();
}
