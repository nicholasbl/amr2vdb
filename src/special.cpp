#include "special.h"

#include "argparse.h"
#include "spdlog/spdlog.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Activate.h>
#include <openvdb/tools/Clip.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/VolumeToMesh.h>

#include <filesystem>
#include <span>
#include <unordered_set>

namespace fs = std::filesystem;

auto extract_iso_cells(openvdb::FloatGrid::Ptr grid,
                       float                   isovalue,
                       float                   tolerance) {
    spdlog::debug("Extracting isovalue {}", isovalue);

    auto xfrmr_keep_range = [&](openvdb::FloatGrid::ValueOnCIter const& iter,
                                openvdb::BoolGrid::Accessor& accessor) {
        auto value = *iter;

        auto out_value = openvdb::math::Abs(value - isovalue) <= tolerance;

        if (iter.isVoxelValue()) {
            accessor.setValue(iter.getCoord(), out_value);
        } else {
            openvdb::CoordBBox bbox;
            iter.getBoundingBox(bbox);
            accessor.getTree()->fill(bbox, out_value);
        }
    };

    auto mask = openvdb::BoolGrid::create();

    mask->setTransform(grid->transformPtr());

    spdlog::debug("Transforming values");

    openvdb::tools::transformValues(
        grid->cbeginValueOn(), *mask, xfrmr_keep_range);

    spdlog::debug("Cleaning up...");

    // Deactivate false
    openvdb::tools::deactivate(*mask, false);

    // Prune inactive
    openvdb::tools::pruneInactive(mask->tree());

    grid->transform().print();
    mask->transform().print();

    spdlog::debug("Clip");
    // Clip grid
    return openvdb::tools::clip(*grid, *mask);
}

int all_iso_merge(Arguments& c) {
    auto fname   = fs::path(toml::find<std::string>(c.root, "input"));
    auto outname = fs::path(toml::find<std::string>(c.root, "output"));

    auto config_node = toml::find(c.root, "all_iso_merge");

    auto min_value = toml::find<int>(config_node, "from");
    auto max_value = toml::find<int>(config_node, "to");

    auto quant_name = toml::find<std::string>(config_node, "quantity");
    auto excludes   = toml::find(config_node, "exclude");

    auto adapt = toml::find_or(config_node, "adaptivity", 0.0);


    auto excludes_array = excludes.as_array();

    std::unordered_set<int> to_exclude;

    for (auto const& value : excludes_array) {
        int v = value.as_integer();
        to_exclude.insert(v);
    }

    openvdb::io::File file(fname);

    file.open();

    auto grid = file.readGrid(quant_name);

    auto float_grid = openvdb::gridPtrCast<openvdb::FloatGrid>(grid);

    if (!float_grid) {
        spdlog::error("Only supports floating-point grids at this time");
        return EXIT_FAILURE;
    }

    spdlog::debug("Min and max: {} {}", min_value, max_value);

    if (min_value >= max_value) return EXIT_SUCCESS;

    std::vector<openvdb::Vec3f> global_points;
    std::vector<float>          global_uv;
    std::vector<openvdb::Vec3I> global_tris;
    std::vector<openvdb::Vec4I> global_quads;

    std::vector<openvdb::Vec3f> local_points;
    std::vector<openvdb::Vec3I> local_tris;
    std::vector<openvdb::Vec4I> local_quads;

    // isosurfaces alone arent going to do the job here.

    for (int iso = min_value; iso <= max_value; iso++) {
        spdlog::debug("Working on level {}", iso);
        if (to_exclude.contains(iso)) continue;

        int offset = global_points.size() + 1;

        // Get clipped grid
        auto clipped = extract_iso_cells(float_grid, iso, .5);

        spdlog::debug("Clipped {} to {}",
                      float_grid->activeVoxelCount(),
                      clipped->activeVoxelCount());

        openvdb::tools::volumeToMesh(
            *clipped, local_points, local_tris, local_quads, iso, adapt);

        spdlog::debug("Extracted {} points, {} tris, {} quads",
                      local_points.size(),
                      local_tris.size(),
                      local_quads.size());

        for (auto& p : local_tris) {
            p.x() += offset;
            p.y() += offset;
            p.z() += offset;

            global_tris.push_back(p);
        }

        for (auto& p : local_quads) {
            p.x() += offset;
            p.y() += offset;
            p.z() += offset;
            p.w() += offset;

            global_quads.push_back(p);
        }

        global_points.insert(
            global_points.end(), local_points.begin(), local_points.end());

        float uv = float(iso - min_value) / float(max_value - min_value);

        for (auto const& p : local_points) {
            global_uv.push_back(uv);
        }
    }

    std::ofstream stream(outname);

    stream << "o exported\n";

    for (auto const& p : global_points) {
        stream << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    }

    for (auto const& p : global_uv) {
        stream << "vt " << p << " " << p << std::endl;
    }

    for (auto const& f : global_tris) {
        stream << "f "                         //
               << f.x() << "/" << f.x() << " " //
               << f.y() << "/" << f.y() << " " //
               << f.z() << "/" << f.z()        //
               << std::endl;
    }

    for (auto const& f : global_quads) {
        stream << "f "                         //
               << f.x() << "/" << f.x() << " " //
               << f.y() << "/" << f.y() << " " //
               << f.z() << "/" << f.z() << " " //
               << f.w() << "/" << f.w()        //
               << std::endl;
    }

    return EXIT_SUCCESS;
}
