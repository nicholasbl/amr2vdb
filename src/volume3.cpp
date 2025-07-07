#include "volume.h"

#include "amr_common.h"
#include "argparse.h"
#include "lzs3d.h"
#include "pltfilereader.h"
#include "postprocess.h"
#include "tricubic.h"

#include "spdlog/spdlog.h"

#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <openvdb/openvdb.h>
#include <openvdb/tools/Clip.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/Filter.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/Mask.h>
#include <openvdb/tools/MultiResGrid.h>


#include <array>
#include <span>

// IMPROVEMENTS:
// Extract levels in parallel. We are almost idle on IO.

/// This is a single quantity
struct SampledGrid {
    std::string       name;
    FloatMultiGridPtr multi_grid = nullptr;
};

/// Extraction result
struct Result {
    openvdb::BBoxd           bbox;
    std::vector<SampledGrid> grids;
};

/// Take an AMR box, and copy selected variables of that box to a VDB accessor
///
/// \param bx AMR box
/// \param a AMR box data
/// \param accessors List of destination accessors
/// \param var_indices List of variable ids to extract
///
static void loop_box(amrex::Box const&                       bx,
                     amrex::Array4<amrex::Real const> const& a,
                     std::span<openvdb::FloatGrid::Accessor> accessors,
                     std::span<int const>                    var_indices,
                     openvdb::CoordBBox const&               keep) {

    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    float last_value = 0.0;

    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {

                if (!keep.isInside(openvdb::Coord(i, j, k))) continue;

                for (int offset = 0; offset < var_indices.size(); offset++) {
                    auto index = var_indices[offset];

                    float value = *a.ptr(i, j, k, index);

                    accessors[offset].setValue({ i, j, k }, value);

                    bool check_inf = !std::isfinite(value);
                    bool check_pre =
                        std::abs(accessors[offset].getValue({ i, j, k }) -
                                 value) > 1e-4f;

                    // spdlog::debug("Values: {}", value);

                    if (check_inf or check_pre) {
                        spdlog::error(
                            "Value not written correctly: {} ({}) stored {}",
                            value,
                            check_inf,
                            accessors[offset].getValue({ i, j, k }));
                        accessors[offset].setValue({ i, j, k }, last_value);
                    } else {
                        // only save good values
                        last_value = value;
                    }
                }
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

    std::unique_ptr<PltFileReader> plt_data;

    /// Requested variable names
    std::vector<std::string> var_names;
    /// Corresponding variable ids in the AMR
    std::vector<int> var_ids;

    ConversionState(VolumeConfig c) : config(c) { }

    ///
    /// Attempt to load an AMR file
    /// \param path Pltfile path
    /// \return True if Pltfile could be opened
    ///
    bool init(std::filesystem::path path) {
        if (!std::filesystem::exists(path)) {
            spdlog::error("Path does not exist: {}", path.string());
            return false;
        }

        plt_data      = std::make_unique<PltFileReader>(path);
        auto plt_vars = plt_data->getVariableList();

        // Figure out which vars were asked for, and where they are in the plt

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

        // Pull plt data in core. The API suggests that this pulls in ALL the
        // data, so watch out
        plt_data->readPlotFileData();

        spdlog::info("Done.");

        return true;
    }

    void extract_level(int                     current_level,
                       std::span<FloatGridPtr> per_var_grid) const {

        // Get all the var-vdbs at the appropriate level, along with
        // accessors
        std::vector<FloatGrid::Accessor> per_var_accessor;

        per_var_accessor.reserve(per_var_grid.size());

        for (auto& mg : per_var_grid) {
            per_var_accessor.push_back(mg->getAccessor());
        }

        auto const& ba = plt_data->getGrid(current_level);

        auto const& mf = plt_data->get_data(current_level);

        // We can now iterate the AMR at this level
        auto iter = amrex::MFIter(mf);

        spdlog::debug("Iterating blocks for {}", current_level);

        auto level_bb = config.bounding_box.value_or(openvdb::CoordBBox::inf());

        // only do this if there is a bounding box. 'inf' uses real values that
        // will overflow.
        if (current_level > 0 and config.bounding_box.has_value()) {
            openvdb::Vec3I l0 = level_bb.min().asVec3I() * (1 << current_level);
            openvdb::Vec3I l1 = level_bb.max().asVec3I() * (1 << current_level);


            level_bb =
                openvdb::CoordBBox(openvdb::Coord(l0), openvdb::Coord(l1));
        }

        spdlog::debug("Level BB: {} {} {} - {} {} {}",
                      level_bb.min().x(),
                      level_bb.min().y(),
                      level_bb.min().z(),
                      level_bb.max().x(),
                      level_bb.max().y(),
                      level_bb.max().z());


        for (; iter.isValid(); ++iter) {
            auto const& box = iter.validbox();

            amrex::FArrayBox const& fab = mf[iter];

            auto const& a = fab.array();

            AMREX_ALWAYS_ASSERT(fab.box().contains(box));

            loop_box(box, a, per_var_accessor, var_ids, level_bb);
        }
    }

    /// Execute a copy from AMR to multires VDBs
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


        // This is a bit silly, but we just make a big list of vdbs that
        // correspond to our list of desired variables
        std::vector<FloatMultiGrid::Ptr> per_var_multivdb_grids;

        for (size_t i = 0; i < var_ids.size(); i++) {
            auto g =
                std::make_shared<FloatMultiGrid>(config.max_level + 1, 0.0f);
            per_var_multivdb_grids.push_back(g);
        }


        // note the <= here
        for (int current_level = 0; current_level <= config.max_level;
             current_level++) {
            spdlog::info("Working on AMR level {}", current_level);
            // also note that VDB does inverse resolution order. 0 is the
            // finest. so we want to map, say 3 -> 0 and 0 -> 3
            size_t vdb_mapped_level = (config.max_level - current_level);
            spdlog::debug("Mapping to VDB level {}", vdb_mapped_level);

            std::vector<FloatGridPtr> local_grids;

            for (auto& g : per_var_multivdb_grids) {
                local_grids.push_back(g->grid(vdb_mapped_level));
            }

            extract_level(current_level, local_grids);
        }

        spdlog::info("Sampling complete");

        // Transform output a bit to simplify things
        std::vector<SampledGrid> ret_grid;

        for (int i = 0; i < per_var_multivdb_grids.size(); i++) {
            ret_grid.emplace_back(SampledGrid {
                .name       = var_names[i],
                .multi_grid = per_var_multivdb_grids[i],
            });
        }

        // If the user asked for an AMR structure grid, execute that and add it
        // in

        if (config.save_amr) { ret_grid.push_back(save_amr()); }

        Result ret;

        ret.grids = std::move(ret_grid);
        ret.bbox  = {
            { (double)larray[0], (double)larray[1], (double)larray[2] },
            { (double)harray[0], (double)harray[1], (double)harray[2] }
        };


        return ret;
    }

    /// Save the AMR struture. This creates a float multi_grid VDB, where every
    /// voxel
    /// is the finest level of the AMR refinement at that spot
    SampledGrid save_amr() const {
        spdlog::info("Building AMR structure grid...");
        auto g = std::make_shared<FloatMultiGrid>(config.max_level + 1, 0.0f);

        for (int current_level = 0; current_level <= config.max_level;
             current_level++) {
            spdlog::info("Working on AMR level {}", current_level);

            size_t vdb_mapped_level = (config.max_level - current_level);

            spdlog::debug("Mapping to VDB level {}", vdb_mapped_level);

            FloatGridPtr        per_var_grid     = g->grid(vdb_mapped_level);
            FloatGrid::Accessor per_var_accessor = per_var_grid->getAccessor();

            auto const& ba = plt_data->getGrid(current_level);

            auto const& mf = plt_data->get_data(current_level);

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

// static void mask_grid(FloatGridPtr mask_grid, FloatGridPtr to_mask) {
//     // make a mask

//     auto resampled = openvdb::BoolGrid::create(false);
//     {
//         auto bool_grid = openvdb::BoolGrid::create(false);
//         bool_grid->setName("mask");
//         bool_grid->setTransform(mask_grid->transform().copy());
//         // bool_grid->topologyUnion(*mask_grid);
//         auto accessor = bool_grid->getAccessor();
//         for (auto iter = mask_grid->cbeginValueOn(); iter; ++iter) {
//             accessor.setValueOn(iter.getCoord());
//         }

//         openvdb::io::File("mask_before.vdb").write({ bool_grid });


//         resampled->setTransform(to_mask->transform().copy());
//         openvdb::tools::resampleToMatch<openvdb::tools::PointSampler>(
//             *bool_grid, *resampled);

//         resampled->setName("mask");
//         openvdb::io::File("mask_after.vdb").write({ resampled });
//     }


//     to_mask->topologyDifference(*resampled);
// }

static void mask_grid(FloatGridPtr mask_grid, FloatGridPtr to_mask) {
    // make a mask
    auto mask = openvdb::BoolGrid::create(false);
    mask->setTransform(mask_grid->transform().copy());
    mask->setGridClass(openvdb::GRID_FOG_VOLUME);

    {
        auto accessor = mask->getAccessor();
        for (auto iter = mask_grid->cbeginValueOn(); iter; ++iter) {
            accessor.setValue(iter.getCoord(), true);
        }
    }

    mask->pruneGrid();

    auto resampled_mask = openvdb::BoolGrid::create(false);
    resampled_mask->setTransform(to_mask->transform().copy());
    openvdb::tools::resampleToMatch<openvdb::tools::BoxSampler>(
        *mask, *resampled_mask);

    // save some memory
    mask->clear();

    int erosion_voxels = 2;
    openvdb::tools::erodeActiveValues(resampled_mask->tree(), erosion_voxels);

    // have the same transforms, so coords line up
    auto mask_accessor = resampled_mask->getAccessor();


    // Deactivate voxels
    auto accessor = to_mask->getAccessor();
    for (auto iter = resampled_mask->cbeginValueOn(); iter.test(); ++iter) {
        auto mask_value = iter.getValue();
        auto mask_coord = iter.getCoord();


        if (mask_value) {
            accessor.setValue(mask_coord, 0.0);
            accessor.setValueOff(mask_coord);
        }
    }

    // openvdb::tools::compMul(*to_mask, *resampled_mask);

    to_mask->pruneGrid(0.0);
}

// Create a sphere SDF as test data
FloatGridPtr
createSphere(float radius, openvdb::Vec3f center, float voxelSize) {
    auto grid      = openvdb::FloatGrid::create(0);
    auto transform = openvdb::math::Transform::createLinearTransform(voxelSize);
    grid->setTransform(transform);

    int  dim = static_cast<int>(radius * 2.2f / voxelSize);
    auto lo  = openvdb::Coord(-dim, -dim, -dim);
    auto hi  = openvdb::Coord(dim, dim, dim);

    auto accessor = openvdb::FloatGrid::Accessor(grid->tree());

    spdlog::info(
        "{} {} {} < {} {} {}", lo.x(), lo.y(), lo.z(), hi.x(), hi.y(), hi.z());

    for (int k = lo.z(); k <= hi.z(); ++k) {
        for (int j = lo.y(); j <= hi.y(); ++j) {
            for (int i = lo.x(); i <= hi.x(); ++i) {
                openvdb::Coord ijk(i, j, k);
                openvdb::Vec3d pos = transform->indexToWorld(ijk);

                double dist = (pos - center).length();
                if (dist < radius) { accessor.setValue(ijk, 1.0); }
            }
        }
    }

    // grid->pruneGrid();

    grid->setGridClass(openvdb::GRID_FOG_VOLUME);
    return grid;
}

int amr_to_volume_sets(Arguments const& c) {
    spdlog::info("Using volume sets");

#if 0
    {
        spdlog::info("Testing...");

        // Create a coarse sphere
        auto coarseGrid = createSphere(
            /*radius=*/50.0f,
            openvdb::Vec3f(0, 0, 0),
            /*voxelSize=*/1.0f);

        // Create a high-res sphere overlapping partially
        auto highResGrid = createSphere(
            /*radius=*/20.0f,
            openvdb::Vec3f(10, 0, 0),
            /*voxelSize=*/0.5f);

        // Apply your masking
        mask_grid(highResGrid, coarseGrid);

        // write grids
        coarseGrid->setName("Coarse");
        highResGrid->setName("High");

        openvdb::GridPtrVec to_save;
        to_save.push_back(coarseGrid);
        to_save.push_back(highResGrid);

        openvdb::io::File file("test.vdb");
        file.write(to_save);
        file.close();

        return EXIT_SUCCESS;
    }
#endif

    std::string source_plt = toml::find<std::string>(c.root, "input");
    std::string dest_vdb   = toml::find<std::string>(c.root, "output");

    auto amr_config_node = toml::find(c.root, "amr3");

    VolumeConfig config;

    auto source_vars = toml::find(amr_config_node, "variables").as_array();

    if (amr_config_node.contains("save_amr")) {
        config.save_amr = toml::find<std::string>(amr_config_node, "save_amr");
    }

    for (auto const& var : source_vars) {
        config.variables.insert(var.as_string());
    }

    if (amr_config_node.contains("clip_coarse")) {
        auto clip_node = toml::find(amr_config_node, "clip_coarse");

        auto min = toml::find(clip_node, "min").as_array();
        auto max = toml::find(clip_node, "max").as_array();

        auto min_v = openvdb::Coord(min.at(0).as_integer(),
                                    min.at(1).as_integer(),
                                    min.at(2).as_integer());

        auto max_v = openvdb::Coord(max.at(0).as_integer(),
                                    max.at(1).as_integer(),
                                    max.at(2).as_integer());

        config.bounding_box = openvdb::CoordBBox(min_v, max_v);
    }

    auto loaded_amr_grids = load_file(source_plt, config);

    spdlog::info("Read: {}", loaded_amr_grids.grids.size());

    GridMap completed;

    for (auto& sampled_grid : loaded_amr_grids.grids) {

        if (!sampled_grid.multi_grid) {
            spdlog::error("No grids for {}", sampled_grid.name);
            spdlog::error("Details: {}", loaded_amr_grids.grids.size());
            return EXIT_FAILURE;
        }

        auto multi = sampled_grid.multi_grid;

        auto max_levels = multi->numLevels();

        if (max_levels < 1) {
            spdlog::error("Unable to handle zero sized levels!");
            exit(1);
        }

        // remember that levels are inverted. the higher you go, the coarser
        // we want to mask coarser with finer
        // now this means we want to process things in a certain order.

        for (int level_i = max_levels - 1; level_i >= 0; --level_i) {
            auto this_grid = multi->grid(level_i);
            auto name      = sampled_grid.name + "_" + std::to_string(level_i);

            auto bb = this_grid->evalActiveVoxelBoundingBox();

            spdlog::info("Creating grid {} with {}: {} {} {} - {} {} {}",
                         name,
                         this_grid->activeVoxelCount(),
                         bb.min().x(),
                         bb.min().y(),
                         bb.min().z(),
                         bb.max().x(),
                         bb.max().y(),
                         bb.max().z());

            if (level_i != 0) {
                mask_grid(multi->grid(level_i - 1), this_grid);

                auto bb = this_grid->evalActiveVoxelBoundingBox();

                spdlog::info("Masked {} to: {} {} {} - {} {} {}",
                             name,
                             this_grid->evalActiveVoxelBoundingBox(),
                             this_grid->activeVoxelCount(),
                             bb.min().x(),
                             bb.min().y(),
                             bb.min().z(),
                             bb.max().x(),
                             bb.max().y(),
                             bb.max().z());
            }

            this_grid->setName(name);

            completed[name] = this_grid;
        }
    }

    if (c.root.contains("post")) { postprocess(c.root, completed); }

    openvdb::GridPtrVec to_save;

    {
        for (auto const& [k, v] : completed) {
            // v->setGridClass(openvdb::v12_0::GRID_FOG_VOLUME);
            to_save.push_back(v);
        }
    }

    openvdb::io::File file(dest_vdb);
    file.write(to_save);
    file.close();

    return EXIT_SUCCESS;
}
