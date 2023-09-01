#include "volume.h"

#include "amr_common.h"
#include "argparse.h"
#include "lzs3d.h"
#include "postprocess.h"
#include "tricubic.h"

#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <PltFileManager.H>

#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/Filter.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/MultiResGrid.h>


#include <array>
#include <span>


enum class SampleType {
    QUAD,
    LANC,
    TRIC,
};

const char* to_string(SampleType type) {
    switch (type) {
    case SampleType::QUAD: return "Triquadratic";
    case SampleType::LANC: return "Lanczos";
    case SampleType::TRIC: return "Tricubic";
    }
}


struct SampledGrid {
    std::string       name;
    FloatMultiGridPtr grid;
};

struct Result {
    openvdb::BBoxd           bbox;
    std::vector<SampledGrid> grids;
};

static void loop_box(amrex::Box const&                       bx,
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

struct ConversionState {
    VolumeConfig  config;

    std::unique_ptr<LocalPlotFile> plt_data;

    std::vector<std::string> var_names;
    std::vector<int>         var_ids;

    ConversionState(VolumeConfig c) : config(c) { }


    bool init(std::filesystem::path path) {
        if (!std::filesystem::exists(path)) {
            std::cerr << "Path does not exist: " << path << std::endl;
            return false;
        }

        plt_data      = std::make_unique<LocalPlotFile>(path);
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
};


static Result load_file(std::filesystem::path path, VolumeConfig const& c) {
    auto state = AMRState();

    auto c_state = ConversionState(c);

    if (!c_state.init(path)) { return {}; }

    return c_state.write_to_vdbs();
}

struct BlurArgs {
    enum Strat {
        NO_BLUR,
        BLUR_AFTER_EVERY_LEVEL,
        BLUR_AFTER_LAST_LEVEL,
        BLUR_AT_END,
    } strategy = NO_BLUR;

    int blur_radius     = 1;
    int blur_iterations = 1;

    void blur_fld(FloatGridPtr ptr) const {
        if (blur_radius <= 0) return;

        std::cout << "Blur radius " << blur_radius << " x" << blur_iterations
                  << std::endl;

        auto filter = openvdb::tools::Filter(*ptr);

        // filter.gaussian(blur_radius, blur_iterations);
        filter.mean(blur_radius, blur_iterations);
    }

    static Strat string_to_strat(std::string text) {
        static std::unordered_map<std::string, Strat> mapper = {
            { "after_every", BLUR_AFTER_EVERY_LEVEL },
            { "after_last", BLUR_AFTER_LAST_LEVEL },
            { "at_end", BLUR_AT_END },
        };

        for (auto& c : text) {
            c = std::tolower(c);
        }

        if (auto iter = mapper.find(text); iter != mapper.end()) {
            return iter->second;
        }

        return NO_BLUR;
    }
};

std::vector<FloatGridPtr> make_grids(int level) {
    std::vector<FloatGridPtr> ret;
    assert(level > 0);

    for (int i = 0; i < level; i++) {
        auto g = FloatGrid::create(0.0);

        auto xform = openvdb::math::Transform::createLinearTransform();

        if (level > 0) xform->preScale(1 << i);

        g->setTransform(xform);

        ret.push_back(g);
    }

    return ret;
}

FloatGridPtr resample(SampledGrid    source,
                      openvdb::BBoxd box,
                      SampleType     type,
                      BlurArgs       opts) {
    std::cout << "Flattening to fine grid...\n";

    // interpolate from coarse to fine to this grid
    int max_levels = source.grid->numLevels();
    int level      = max_levels;

    // create new grid
    auto new_multi_grid = make_grids(max_levels);

    while (level-- > 0) {
        std::cout << "Packing " << level << std::endl;

        auto in_grid  = source.grid->grid(level);
        auto out_grid = new_multi_grid.at(level);

        std::cout << "In " << in_grid->activeVoxelCount() << std::endl;

        if (level != max_levels - 1) {
            auto previous_grid = new_multi_grid.at(level + 1);

            std::cout << "Interpolate " << (level + 1) << " to " << level
                      << std::endl;

            std::cout << "Previous has " << previous_grid->activeVoxelCount()
                      << std::endl;

            switch (type) {
            case SampleType::QUAD:
                openvdb::tools::resampleToMatch<
                    openvdb::tools::QuadraticSampler>(*previous_grid,
                                                      *out_grid);
                break;
            case SampleType::LANC:
                openvdb::tools::resampleToMatch<LZS3D>(*previous_grid,
                                                       *out_grid);
                break;
            case SampleType::TRIC:
                openvdb::tools::resampleToMatch<Tricubic>(*previous_grid,
                                                          *out_grid);
                break;
            }

            std::cout << "Out prepared with " << out_grid->activeVoxelCount()
                      << std::endl;
        }

        /*
        // at end blur
        if (level == 0 and opts.blur_radius > 0 and opts.blur_iterations >= 1) {
            std::cout << "Blur radius " << opts.blur_radius << " x"
                      << opts.blur_iterations << std::endl;
            auto in_grid = new_grid->grid(level + 1);

            auto filter = openvdb::tools::Filter(*in_grid);

            filter.gaussian(opts.blur_radius, opts.blur_iterations);
        }
        */


        // copy over new data
        openvdb::tools::compReplace(*out_grid, *source.grid->grid(level));

        // smoother but overall more shitty.
        //        openvdb::tools::resampleToMatch<openvdb::tools::BoxSampler>(
        //            *source.grid->grid(level), *out_grid);

        std::cout << "Storing to " << level << std::endl;
        std::cout << "Out resampled with " << out_grid->activeVoxelCount()
                  << std::endl;


        if (level != 0 and opts.strategy == BlurArgs::BLUR_AFTER_EVERY_LEVEL) {

            opts.blur_fld(out_grid);
        } else if (level == 1 and
                   opts.strategy == BlurArgs::BLUR_AFTER_LAST_LEVEL) {
            opts.blur_fld(out_grid);
        }
    }

    auto ret_grid = new_multi_grid.at(0);

    if (opts.strategy == BlurArgs::BLUR_AT_END) { opts.blur_fld(ret_grid); }

    std::cout << "Final " << ret_grid->activeVoxelCount() << std::endl;

    ret_grid->clipGrid(box);
    ret_grid->setName(source.name);
    return ret_grid;
}

int amr_to_volume(Arguments const& c) {
    std::string source_plt = toml::find<std::string>(c.root, "input");
    std::string dest_vdb   = toml::find<std::string>(c.root, "output");

    auto amr_config_node = toml::find(c.root, "amr");

    VolumeConfig config;

    auto source_vars = toml::find(amr_config_node, "variables").as_array();

    for (auto const& var : source_vars) {
        config.variables.insert(var.as_string());
    }

    auto amr = load_file(source_plt, config);

    auto sample_type =
        toml::find_or<std::string>(amr_config_node, "sample", "triquad");

    SampleType type = [&] {
        if (sample_type == "lanczos") return SampleType::LANC;
        if (sample_type == "tricubic") return SampleType::TRIC;
        return SampleType::QUAD;
    }();

    std::cout << "Using sample method " << to_string(type) << std::endl;

    auto blur_config_node =
        toml::find_or(amr_config_node, "blur", toml::value());

    BlurArgs opts {
        .strategy = BlurArgs::string_to_strat(
            toml::find_or(blur_config_node, "strategy", "none")),
        .blur_radius = toml::find_or<int>(blur_config_node, "radius", -1),
        .blur_iterations =
            toml::find_or<int>(blur_config_node, "iterations", 0),
    };

    GridMap completed;

    for (auto& multires : amr.grids) {
        auto new_grid = resample(multires, amr.bbox, type, opts);
        multires.grid.reset(); // try to minimize mem usage
        completed[new_grid->getName()] = new_grid;
    }

    if (c.root.contains("post")) {
        postprocess(PostProcessOptions::from_toml(c.root), completed);
    }

    openvdb::GridPtrVec to_save;

    {
        for (auto const& [k, v] : completed) {
            to_save.push_back(v);
        }
    }

    openvdb::io::File file(dest_vdb);
    file.write(to_save);
    file.close();

    return EXIT_SUCCESS;
}
