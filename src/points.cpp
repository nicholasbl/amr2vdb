#include "points.h"
#include "argparse.h"
#include "spdlog/spdlog.h"

#include <openvdb/openvdb.h>
#include <openvdb/points/PointConversion.h>
#include <openvdb/points/PointCount.h>
#include <openvdb/tools/LevelSetUtil.h>
#include <openvdb/tools/ParticlesToLevelSet.h>
#include <toml.hpp>

#include <cstdlib>

namespace vdb = openvdb;
namespace pts = openvdb::points;

// wrapper class
struct Points {
    using PosType = vdb::Vec3R;

    Points(const std::vector<vdb::Vec3s>& vtx) : mPoints(vtx) { }
    size_t size() const { return mPoints.size(); }
    void   getPos(size_t n, PosType& p) const { p = mPoints[n]; }

    const std::vector<vdb::Vec3s>& mPoints;
};

template <class Codec, class ValueType, class GridType, class IndexGrid>
void add_attrib_to_grid(GridType               grid,
                        IndexGrid              point_index_grid,
                        std::vector<ValueType> list,
                        std::string            name) {
    pts::appendAttribute(
        grid->tree(),
        name,
        pts::TypedAttributeArray<ValueType, Codec>::attributeType());

    auto id_wrapper = pts::PointAttributeVector(list);

    pts::populateAttribute(
        grid->tree(), point_index_grid->tree(), name, id_wrapper);
}

int txt_to_points(Arguments const& c) {
    std::string dest_vdb = toml::find<std::string>(c.root, "output");

    auto point_config_node = toml::find(c.root, "points");

    auto file_name = toml::find(point_config_node, "input").as_string();
    auto radius    = toml::find(point_config_node, "radius");

    std::string radius_var;
    float       radius_in_min  = 0;
    float       radius_in_max  = 1;
    float       radius_out_min = 0;
    float       radius_out_max = 1;

    if (radius.is_string()) {
        radius_var = radius.as_string();
    } else if (radius.is_table()) {
        radius_var     = toml::find<std::string>(radius, "attribute");
        radius_in_min  = toml::find_or<float>(radius, "radius_in_min", 0.0);
        radius_in_max  = toml::find_or<float>(radius, "radius_in_max", 1.0);
        radius_out_min = toml::find_or<float>(radius, "radius_out_min", 0.0);
        radius_out_max = toml::find_or<float>(radius, "radius_out_max", 1.0);
    }

    std::ifstream file(file_name);

    std::string cache_line;
    std::string cache_token;

    // consume header

    // header is in form of
    // NUM_POINTS EXTRA_ATTRIB_NAME_1 EXTRA_ATTRIB_NAME_2 EXTRA_ATTRIB_NAME_3
    // id x y z a1 a2 a3

    uint64_t num_points = 0;

    std::vector<std::string> attribs;


    {
        std::getline(file, cache_line);

        std::istringstream ss(cache_line);

        std::getline(ss, cache_token, ' ');

        num_points = std::stoull(cache_token);

        while (std::getline(ss, cache_token, ' ')) {
            if (cache_token.empty()) continue;

            attribs.push_back(cache_token);
        }
    }

    const size_t attrib_count = attribs.size();

    spdlog::info("Reading {} points with {} attribs", num_points, attrib_count);

    std::vector<int>            ids;
    std::vector<openvdb::Vec3f> positions;

    std::vector<std::vector<float>> attrib_data(attrib_count);

    ids.reserve(num_points);
    positions.reserve(num_points);

    for (auto& a : attrib_data) {
        a.reserve(num_points);
    }

    while (std::getline(file, cache_line)) {
        std::istringstream ss(cache_line);

        int   id;
        float x, y, z;
        ss >> id >> x >> y >> z;

        ids.push_back(id);
        positions.push_back({ x, y, z });

        for (auto i = 0; i < attrib_count; i++) {
            float value = 0.0;
            ss >> value;
            attrib_data[i].push_back(value);
        }
    }

    if (!radius_var.empty()) {

        auto iter = std::find(attribs.begin(), attribs.end(), radius_var);

        if (iter != attribs.end()) {
            auto src_index = std::distance(attribs.begin(), iter);

            auto const& src_data = attrib_data[src_index];

            std::vector<float> rads;

            rads.reserve(src_data.size());

            auto delta = (radius_out_max - radius_out_min) /
                         (radius_in_max - radius_in_min);

            for (auto& d : src_data) {
                rads.push_back(radius_out_min + (d - radius_in_min) * (delta));
            }

            attribs.push_back("radius");
            attrib_data.emplace_back(rads);
        }
    }

    // points read, create a point vdb

    pts::PointAttributeVector<openvdb::Vec3f> positions_wrapper(positions);

    int   points_per_voxel = 8;
    float voxel_size =
        pts::computeVoxelSize(positions_wrapper, points_per_voxel);

    spdlog::info("Using a voxel size of {}", voxel_size);

    auto transform =
        openvdb::math::Transform::createLinearTransform(voxel_size);

    auto point_index_grid =
        openvdb::tools::createPointIndexGrid<openvdb::tools::PointIndexGrid>(
            positions_wrapper, *transform);

    auto grid = pts::createPointDataGrid<pts::NullCodec, pts::PointDataGrid>(
        *point_index_grid, positions_wrapper, *transform);

    grid->setName("Points");

    add_attrib_to_grid<pts::NullCodec>(grid, point_index_grid, ids, "id");

    for (size_t i = 0; i < attrib_count; i++) {
        auto const& name = attribs[i];
        auto const& data = attrib_data[i];

        add_attrib_to_grid<pts::TruncateCodec>(
            grid, point_index_grid, attrib_data[i], name);
    }

    auto level_set_info =
        toml::find_or(point_config_node, "level", toml::value());

    if (level_set_info.is_table()) {
        std::string attrib =
            toml::find<std::string>(level_set_info, "attribute");
        float voxel_size =
            toml::find_or<float>(level_set_info, "voxel_size", -1);
        int dim = toml::find_or<int>(level_set_info, "voxel_size", 256);

        float width = toml::find_or<int>(level_set_info, "width", 3.0);


        if (voxel_size <= 0) {
            // use dims

            auto box = grid->evalActiveVoxelBoundingBox();

            auto bb_min = box.min().asVec3s();
            auto bb_max = box.max().asVec3s();

            spdlog::info("Points bounding box: {} {} {} <-> {} {} {}",
                         bb_min.x(),
                         bb_min.y(),
                         bb_min.z(),
                         bb_max.x(),
                         bb_max.y(),
                         bb_max.z());

            auto d = box.extents()[box.maxExtent()];

            voxel_size = (float)d / (float)(dim - 2.0f * width);
        }

        auto ls_grid =
            openvdb::createLevelSet<openvdb::FloatGrid>(voxel_size, width);

        ls_grid->setName(grid->getName() + '_' + attrib);


        openvdb::tools::particlesToSdf(Points(positions), *ls_grid, voxel_size);

        // openvdb::tools::sdfToFogVolume(*ls_grid);

        openvdb::io::File(dest_vdb).write({ ls_grid });

        spdlog::info("Writing level set");

        return EXIT_SUCCESS;
    }


    spdlog::info("Writing point set");
    openvdb::io::File(dest_vdb).write({ grid });

    return EXIT_SUCCESS;
}
