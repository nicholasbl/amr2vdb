#include "volume.h"

#include "amr_basic.h"
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

namespace volume {

/// Type of sampling to use during interpolation
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

/// Intermediate multiresolution grid, ready for interpolation
struct SampledGrid {
    std::string       name;
    FloatMultiGridPtr multi_grid;
    FloatGridPtr      plain_grid;
};

/// Extraction result
struct Result {
    openvdb::BBoxd bbox;
    SampledGrid    grid;
};

/// Take an AMR box, and copy selected variables of that box to an AMR accessor
///
/// \param bx AMR box
/// \param a AMR box data
/// \param accessors List of destination accessors
/// \param var_indices List of variable ids to extract
///
static void loop_box(amrex::Box const&                       bx,
                     amrex::Array4<amrex::Real const> const& a,
                     openvdb::FloatGrid::Accessor            accessor,
                     int const                               var_index) {

    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                float value = *a.ptr(i, j, k, var_index);

                accessor.setValue({ i, j, k }, value);

                assert(accessors[offset].getValue({ i, j, k }) == value);
            }
        }
    }
}

/// Take an AMR box, and, using the box bounds, write a constant to a VDB grid
///
/// \param bx AMR box
/// \param value Constant to write
/// \param accessors VDB grid destination
///
static void loop_box_constant(amrex::Box const&             bx,
                              float                         value,
                              openvdb::FloatGrid::Accessor& accessors) {

    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {

                accessors.setValue({ i, j, k }, value);
            }
        }
    }
}

///
/// The ConversionState class models the settings and state for a conversion
///
struct ConversionState {
    VolumeConfig config;

    amr_basic::LocalPlotFilePtr plt_data;

    /// Requested variable name
    std::string var_name;
    /// Corresponding variable id in the AMR
    int var_id = -1;

    ConversionState(VolumeConfig c) : config(c) { }

    ///
    /// Attempt to load an AMR file
    /// \param path Pltfile path
    /// \return True if Pltfile could be opened
    ///
    bool init(amr_basic::LocalPlotFilePtr ptr) {
        spdlog::info("Initializing extraction");
        plt_data      = ptr;
        auto plt_vars = amr_basic::get_variable_list(plt_data);

        // Figure out which var we were asked for, and where they are in the plt
        for (long i = 0; i < plt_vars.size(); i++) {
            auto const& vname = plt_vars[i];
            if (config.variable == vname) {
                spdlog::info("Variable {} is at index {}", vname, i);
                var_name = vname;
                var_id   = i;
            }
        }

        if (var_id < 0) {
            spdlog::error("Unable to find variable {}", config.variable);
            return false;
        }

        if (config.max_level < 0) {
            config.max_level = amr_basic::get_manager(plt_data).getNlev() - 1;
        }

        if (config.max_level >= amr_basic::get_manager(plt_data).getNlev()) {
            spdlog::error(
                "Asking for refinement level {} but pltfile only has {}",
                config.max_level,
                amr_basic::get_manager(plt_data).getNlev());
            return false;
        }

        // Pull plt data in core. The API suggests that this pulls in ALL the
        // data, so watch out
        amr_basic::get_manager(plt_data).readPlotFileData();

        spdlog::info("Done.");

        return true;
    }

    void extract_level(int current_level, FloatGridPtr var_grid) {

        // Get all the var-vdbs at the appropriate level, along with
        // accessors
        auto per_var_accessor = var_grid->getAccessor();

        auto const& ba =
            amr_basic::get_manager(plt_data).getGrid(current_level);

        auto const& mf = amr_basic::get_multifab_at(plt_data, current_level);

        // We can now iterate the AMR at this level
        auto iter = amrex::MFIter(mf);

        spdlog::debug("Iterating blocks...");

        for (; iter.isValid(); ++iter) {
            auto const& box = iter.validbox();

            amrex::FArrayBox const& fab = mf[iter];

            auto const& a = fab.array();

            loop_box(box, a, per_var_accessor, var_id);
        }
    }

    /// Execute a copy from AMR to multires VDBs
    Result write_to_vdbs() {
        spdlog::info("Extracting AMR to Multires VDB");
        auto geom = get_manager(plt_data).getGeom(config.max_level);

        auto domain = geom.Domain();
        auto larray = domain.smallEnd().toArray();
        auto harray = domain.bigEnd().toArray();

        auto grid = amrex::BoxArray(domain);

        spdlog::info("Level:    {}", config.max_level);
        spdlog::info("Domain L: {} {} {}", larray[0], larray[1], larray[2]);
        spdlog::info("Domain H: {} {} {}", harray[0], harray[1], harray[2]);
        spdlog::info("Var Num:  {}", get_variable_list(plt_data).size());

        SampledGrid ret_grid;

        if (config.max_level > 0) {
            // we have to do some resampling

            auto multivdb_grid =
                std::make_shared<FloatMultiGrid>(config.max_level + 1, 0.0f);

            // note the <= here
            for (int current_level = 0; current_level <= config.max_level;
                 current_level++) {
                spdlog::info("Working on AMR level {}", current_level);
                // also note that VDB does inverse resolution order. 0 is the
                // finest. so we want to map, say 3 -> 0 and 0 -> 3
                size_t vdb_mapped_level = (config.max_level - current_level);
                spdlog::debug("Mapping to VDB level {}", vdb_mapped_level);

                FloatGridPtr local_grid = multivdb_grid->grid(vdb_mapped_level);

                extract_level(current_level, local_grid);
            }

            spdlog::info("Sampling complete");

            // Transform output a bit to simplify things

            ret_grid = SampledGrid {
                var_name,
                multivdb_grid,
            };
        } else {
            // we shall be doing no resampling

            auto local_grid = FloatGrid::create(0.0f);

            extract_level(0, local_grid);

            ret_grid = SampledGrid {
                .name       = var_name,
                .plain_grid = local_grid,
            };
        }

        // If the user asked for an AMR structure grid, execute that and add it
        // in

        // if (config.save_amr) { ret_grid.push_back(save_amr()); }

        Result ret;

        ret.grid = ret_grid;
        ret.bbox = {
            { (double)larray[0], (double)larray[1], (double)larray[2] },
            { (double)harray[0], (double)harray[1], (double)harray[2] }
        };


        return ret;
    }

    /// Save the AMR struture. This creates a float multi_grid VDB, where every
    /// voxel
    /// is the finest level of the AMR refinement at that spot
    SampledGrid save_amr() {
        spdlog::info("Building AMR structure grid...");
        auto g = std::make_shared<FloatMultiGrid>(config.max_level + 1, 0.0f);

        for (int current_level = 0; current_level <= config.max_level;
             current_level++) {
            spdlog::info("Working on AMR level {}", current_level);

            size_t vdb_mapped_level = (config.max_level - current_level);

            spdlog::debug("Mapping to VDB level {}", vdb_mapped_level);

            FloatGridPtr        per_var_grid     = g->grid(vdb_mapped_level);
            FloatGrid::Accessor per_var_accessor = per_var_grid->getAccessor();

            auto const& ba = get_manager(plt_data).getGrid(current_level);

            auto const& mf = get_multifab_at(plt_data, current_level);

            auto iter = amrex::MFIter(mf);

            spdlog::debug("Iterating blocks...");

            for (; iter.isValid(); ++iter) {
                auto const& box = iter.validbox();

                amrex::FArrayBox const& fab = mf[iter];

                auto const& a = fab.array();

                loop_box_constant(box, vdb_mapped_level + 1, per_var_accessor);
            }

            openvdb::tools::prune(per_var_grid->tree());
        }

        return {
            .name       = config.save_amr.value(),
            .multi_grid = g,
        };
    }
};

///
/// Wrapper function to pull in a PLT file
///
static Result load_file(std::filesystem::path path, VolumeConfig const& c) {
    auto state = AMRState();

    auto c_state = ConversionState(c);

    if (!c_state.init(path)) { return {}; }

    return c_state.write_to_vdbs();
}

///
/// The BlurArgs class holds smoothing/blurring options
///
struct BlurArgs {
    enum Strat {
        NO_BLUR,                // No blur requested
        BLUR_AFTER_EVERY_LEVEL, // As we process levels, blur each time
        BLUR_AFTER_LAST_LEVEL,  // Only blur the last coarsest level
        BLUR_AT_END,            // Blur the whole thing at the end
    } strategy = NO_BLUR;

    int blur_radius     = 1;
    int blur_iterations = 1;

    /// Read blur options from a TOML node
    static BlurArgs from_toml(toml::value const& val) {
        return BlurArgs {
            .strategy = BlurArgs::string_to_strat(
                toml::find_or(val, "strategy", "none")),
            .blur_radius     = toml::find_or<int>(val, "radius", -1),
            .blur_iterations = toml::find_or<int>(val, "iterations", 0),
        };
    }

    /// Execute a blur on a given grid
    void blur_fld(FloatGridPtr ptr) const {
        if (blur_radius <= 0) return;

        spdlog::info("Blur radius {}x{}", blur_radius, blur_iterations);

        auto filter = openvdb::tools::Filter(*ptr);

        // filter.gaussian(blur_radius, blur_iterations);
        filter.mean(blur_radius, blur_iterations);
    }

    /// Turing a string into a strategy
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

/// Make a multiresolution VDB (which is just a list of grids with different
/// transforms)
std::vector<FloatGridPtr> make_grids(int level) {
    std::vector<FloatGridPtr> ret;
    assert(level > 0);

    for (int i = 0; i < level; i++) {
        auto g = FloatGrid::create(0.0);

        auto xform = openvdb::math::Transform::createLinearTransform();

        // each level is a power of two refinement
        if (level > 0) xform->preScale(1 << i);

        g->setTransform(xform);

        ret.push_back(g);
    }

    return ret;
}


/// Take an AMR styled VDB multi_grid and resample into a single uniform
/// multi_grid.
///
/// This requires care. AMR only has the blocks that matter at a refinement
/// level. VDB expects their multires data to be like a mipmap, where all the
/// coarse levels are mirrors (at lower res) of the finest level. So to make a
/// uniform multi_grid with the tools VDB has, we go level by level,
/// interpolating coarser information first, then copying over the AMR data at
/// that level.
///
/// This function can also do some smoothing; coarse data can get pretty blocky,
/// so we can optionally smooth levels or at the end. Smoothing does change the
/// domain, so a clip box is provided to restrict the output VDB
///
/// \param source Multires multi_grid to pack into a uniform multi_grid
/// \param box Box to clip the data to
/// \param type Interpolation sampling type
/// \param opts Options for smoothing/blurring
/// \return A uniform multi_grid at the finest resolution of the input
///
FloatGridPtr resample(SampledGrid    source,
                      openvdb::BBoxd box,
                      SampleType     type,
                      BlurArgs       opts) {
    spdlog::info("Resampling {} to fine grid", source.name);

    // interpolate from coarse to fine to this multi_grid
    int max_levels = source.multi_grid->numLevels();
    int level      = max_levels;

    // create new grid
    auto new_multi_grid = make_grids(max_levels);

    while (level-- > 0) {
        spdlog::info("Packing {}", level);

        auto in_grid  = source.multi_grid->grid(level);
        auto out_grid = new_multi_grid.at(level);

        spdlog::debug("In {}", in_grid->activeVoxelCount());

        // We first check if there is a coarser grid, and if so, interpolate
        // that into this grid. This a bit like a painter's algorithm approach;
        // interpolate first, then copy in data for this level

        if (level != max_levels - 1) {
            // Take from a coarser grid, and interpolate into this grid
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
        openvdb::tools::compReplace(*out_grid, *source.multi_grid->grid(level));

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

    if (amr_config_node.contains("save_amr")) {
        config.save_amr = toml::find<std::string>(amr_config_node, "save_amr");
    }

    for (auto const& var : source_vars) {
        config.variables.insert(var.as_string());
    }

    auto loaded_amr_grids = load_file(source_plt, config);

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

    for (auto& multires : loaded_amr_grids.grids) {
        // check if there is a per-variable override

        if (multires.multi_grid) {

            BlurArgs local_blur_args = opts;

            if (amr_config_node.contains(multires.name)) {
                spdlog::info("Using custom options for variable {}",
                             multires.name);

                auto local_config = toml::find(amr_config_node, multires.name);

                auto local_blur_config =
                    toml::find_or(local_config, "blur", toml::value());

                if (local_blur_config.is_table()) {
                    auto base = blur_config_node;
                    merge_values(base, local_blur_config);

                    local_blur_args = BlurArgs::from_toml(base);
                }
            }

            auto new_grid = resample(
                multires, loaded_amr_grids.bbox, type, local_blur_args);
            multires.multi_grid.reset(); // try to minimize mem usage
            completed[new_grid->getName()] = new_grid;
        }

        if (multires.plain_grid) {
            multires.plain_grid->setName(multires.name);
            completed[multires.name] = multires.plain_grid;
        }
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

void register_functions(sol::state&) { }

} // namespace volume
