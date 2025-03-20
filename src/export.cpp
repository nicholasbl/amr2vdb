#include "export.h"
#include "argparse.h"
#include "spdlog/spdlog.h"

#include <openvdb/openvdb.h>

namespace vdb = openvdb;

static_assert(sizeof(vdb::Vec3f) == sizeof(float) * 3);
static_assert(sizeof(vdb::Vec3i) == sizeof(int) * 3);

void write_meta(vdb::GridBase::Ptr ptr, std::ofstream& out) {
    auto bb = ptr->evalActiveVoxelBoundingBox();

    auto min = bb.min().asVec3i();
    auto max = bb.max().asVec3i();

    out.write(reinterpret_cast<char*>(&min), sizeof(min));
    out.write(reinterpret_cast<char*>(&max), sizeof(max));
}

template <class T>
int write_volume(vdb::GridBase::Ptr ptr,
                 std::string const& type,
                 std::string const& outname) {
    auto cast_grid = vdb::gridPtrCast<T>(ptr);
    auto accessor  = cast_grid->getAccessor();

    if (!cast_grid) {
        spdlog::error("Cannot cast grid of type {} to {}", ptr->type(), type);
        return EXIT_FAILURE;
    }

    auto bbox = cast_grid->evalActiveVoxelBoundingBox();
    auto dims = cast_grid->evalActiveVoxelDim();

    auto min = bbox.min();
    auto max = bbox.max();

    spdlog::info(
        "Writing grid with dimensions {} {} {}", dims.x(), dims.y(), dims.z());


    auto outfile = std::ofstream(outname, std::ios::binary);

    if (!outfile) {
        spdlog::error("Cannot open outfile {}", outname);
        return EXIT_FAILURE;
    }

    write_meta(ptr, outfile);

    constexpr auto BUFFER_SIZE = 1024 * 1024;

    using QType =
        std::remove_cvref_t<decltype(accessor.getValue(vdb::Coord(0, 0, 0)))>;

    std::vector<QType> buffer;
    buffer.resize(BUFFER_SIZE);


    size_t buffer_index = 0;


    auto flush = [&outfile, &buffer, &buffer_index]() {
        outfile.write(reinterpret_cast<char*>(buffer.data()),
                      buffer_index * sizeof(float));
        buffer_index = 0;
    };

    for (int z = min.z(); z <= max.z(); ++z) {
        for (int y = min.y(); y <= max.y(); ++y) {
            for (int x = min.x(); x <= max.x(); ++x) {
                openvdb::Coord coord(x, y, z);

                auto value = accessor.isValueOn(coord)
                                 ? accessor.getValue(coord)
                                 : QType();

                buffer[buffer_index++] = value;

                if (buffer_index >= BUFFER_SIZE) { flush(); }
            }
        }
    }

    if (buffer_index > 0) { flush(); }

    outfile.close();

    spdlog::info("Completed export.");

    return EXIT_SUCCESS;
}

int export_quantity(Arguments const& c) {
    auto export_node = toml::find(c.root, "export");

    auto filename = toml::find(export_node, "input").as_string();
    auto quantity = toml::find(export_node, "quantity").as_string();
    auto type     = toml::find(export_node, "as").as_string();
    auto outname  = toml::find(export_node, "to").as_string();

    vdb::io::File source(filename);

    source.open();

    if (!source.isOpen()) {
        spdlog::error("Unable to open {}", filename.str);
        return EXIT_FAILURE;
    }

    auto ptr = source.readGrid(quantity);

    source.close();

    if (!ptr) {
        spdlog::error("Missing quantity {}", quantity.str);
        return EXIT_FAILURE;
    }

    if (type == "f32") {
        return write_volume<vdb::FloatGrid>(ptr, type.str, outname.str);
    } else if (type == "vec3") {
        return write_volume<vdb::Vec3SGrid>(ptr, type.str, outname.str);
    }

    spdlog::error("Unknown type {}", type.str);

    return EXIT_FAILURE;
}
