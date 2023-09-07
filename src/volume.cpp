#include "volume.h"

#include "amr_common.h"
#include "argparse.h"
#include "lzs3d.h"
#include "postprocess.h"
#include "tricubic.h"

#include "spdlog/spdlog.h"

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
    __builtin_unreachable();
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
            spdlog::error("Path does not exist: {}", path.string());
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
                    spdlog::error("Variable {} does not exist in pltfile.",
                                  vname);
                    return false;
                }
            }
        }

        spdlog::info("Loading from {}", path.string());

        for (long i = 0; i < plt_vars.size(); i++) {
            auto const& vname = plt_vars[i];
            if (config.variables.contains(vname)) {
                spdlog::info("Variable {} is at index {}", vname, i);
                var_names.push_back(vname);
                var_ids.push_back(i);
            }
        }

        if (config.max_level < 0) {
            config.max_level = plt_data->getNlev() - 1;
        }

        if (config.max_level >= plt_data->getNlev()) {
            spdlog::error(
                "Asking for refinement level {} but pltfile only has {}",
                config.max_level,
                plt_data->getNlev());
            return false;
        }

        plt_data->readPlotFileData();

        spdlog::info("Done.");

        return true;
    }

    Result write_to_vdbs() const {
        spdlog::info("Extracting AMR to Multires VDB");
        auto geom = plt_data->getGeom(config.max_level);

        auto domain = geom.Domain();
        auto larray = domain.smallEnd().toArray();
        auto harray = domain.bigEnd().toArray();

        auto grid = amrex::BoxArray(domain);

        spdlog::info("Level:    {}", config.max_level);
        spdlog::info("Domain L: {} {} {}", larray[0], larray[1], larray[2]);
        spdlog::info("Domain H: {} {} {}", harray[0], harray[1], harray[2]);
        spdlog::info("Var Num:  {}", plt_data->getVariableList().size());

        std::vector<FloatMultiGrid::Ptr> vdb_grids;

        for (size_t i = 0; i < var_ids.size(); i++) {
            auto g =
                std::make_shared<FloatMultiGrid>(config.max_level + 1, 0.0f);
            vdb_grids.push_back(g);
        }


        // note the <= here
        for (int current_level = 0; current_level <= config.max_level;
             current_level++) {
            spdlog::info("Working on AMR level {}", current_level);
            // also note that VDB does inverse sampling. 0 is the finest.
            // so we want to map, say 3 -> 0 and 0 -> 3

            size_t vdb_mapped_level = (config.max_level - current_level);

            spdlog::debug("Mapping to VDB level {}", vdb_mapped_level);

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

            spdlog::debug("Iterating blocks...");

            for (; iter.isValid(); ++iter) {
                auto const& box = iter.validbox();

                amrex::FArrayBox const& fab = mf[iter];

                auto const& a = fab.array();

                loop_box(box, a, per_var_accessor, var_ids);
            }
        }

        spdlog::info("Sampling complete");

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

    static BlurArgs from_toml(toml::value const& val) {
        return BlurArgs {
            .strategy = BlurArgs::string_to_strat(
                toml::find_or(val, "strategy", "none")),
            .blur_radius     = toml::find_or<int>(val, "radius", -1),
            .blur_iterations = toml::find_or<int>(val, "iterations", 0),
        };
    }

    void blur_fld(FloatGridPtr ptr) const {
        if (blur_radius <= 0) return;

        spdlog::info("Blur radius {}x{}", blur_radius, blur_iterations);

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
    spdlog::info("Resampling {} to fine grid", source.name);

    // interpolate from coarse to fine to this grid
    int max_levels = source.grid->numLevels();
    int level      = max_levels;

    // create new grid
    auto new_multi_grid = make_grids(max_levels);

    while (level-- > 0) {
        spdlog::info("Packing {}", level);

        auto in_grid  = source.grid->grid(level);
        auto out_grid = new_multi_grid.at(level);

        spdlog::debug("In {}", in_grid->activeVoxelCount());

        if (level != max_levels - 1) {
            auto previous_grid = new_multi_grid.at(level + 1);

            spdlog::debug("Interpolate {} to {}", (level + 1), level);

            spdlog::debug("Previous has {}", previous_grid->activeVoxelCount());

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

            spdlog::debug("Out prepared with {}", out_grid->activeVoxelCount());
        }


        // copy over new data
        openvdb::tools::compReplace(*out_grid, *source.grid->grid(level));

        // smoother but overall more shitty.
        //        openvdb::tools::resampleToMatch<openvdb::tools::BoxSampler>(
        //            *source.grid->grid(level), *out_grid);

        spdlog::debug("Storing to level {}", level);
        spdlog::debug("Out resampled with {}", out_grid->activeVoxelCount());


        if (level != 0 and opts.strategy == BlurArgs::BLUR_AFTER_EVERY_LEVEL) {

            opts.blur_fld(out_grid);
        } else if (level == 1 and
                   opts.strategy == BlurArgs::BLUR_AFTER_LAST_LEVEL) {
            opts.blur_fld(out_grid);
        }
    }

    auto ret_grid = new_multi_grid.at(0);

    if (opts.strategy == BlurArgs::BLUR_AT_END) { opts.blur_fld(ret_grid); }

    spdlog::debug("Final ", ret_grid->activeVoxelCount());

    ret_grid->clipGrid(box);
    ret_grid->setName(source.name);
    ret_grid->setGridClass(openvdb::GridClass::GRID_FOG_VOLUME);

    openvdb::tools::prune(ret_grid->tree());

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

    spdlog::info("Using sample method {}", to_string(type));

    auto blur_config_node =
        toml::find_or(amr_config_node, "blur", toml::value());

    BlurArgs opts = BlurArgs::from_toml(blur_config_node);

    GridMap completed;

    for (auto& multires : amr.grids) {
        // check if there is a per-variable override

        BlurArgs local_blur_args = opts;

        if (amr_config_node.contains(multires.name)) {
            spdlog::info("Using custom options for variable {}", multires.name);

            auto local_config = toml::find(amr_config_node, multires.name);

            auto local_blur_config =
                toml::find_or(local_config, "blur", toml::value());

            if (local_blur_config.is_table()) {
                auto base = blur_config_node;
                merge_values(base, local_blur_config);

                local_blur_args = BlurArgs::from_toml(base);
            }
        }

        auto new_grid = resample(multires, amr.bbox, type, local_blur_args);
        multires.grid.reset(); // try to minimize mem usage
        completed[new_grid->getName()] = new_grid;
    }

    if (c.root.contains("post")) { postprocess(c.root, completed); }

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
