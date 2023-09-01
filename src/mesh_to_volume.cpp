#include "mesh_to_volume.h"

#include "argparse.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/MeshToVolume.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <filesystem>
#include <span>

namespace fs = std::filesystem;
namespace pt = boost::property_tree;

using openvdb::Vec3I;
using openvdb::Vec3s;
using openvdb::Vec4I;

struct Mesh {
    std::vector<Vec3s> points;
    std::vector<Vec3I> triangles;
    std::vector<Vec4I> quads;

    // MeshDataAdapter

    size_t polygonCount() const { return triangles.size(); }
    size_t pointCount() const { return points.size(); }
    size_t vertexCount(size_t n) const { return 3; }

    void getIndexSpacePoint(size_t n, size_t v, openvdb::Vec3d& pos) const {
        pos = points[triangles[n][v]] * 100;
    }
};

inline void require(bool value, std::string_view expr) {
    if (!value) {
        std::cerr << expr << std::endl;
        exit(EXIT_FAILURE);
    }
}

using ByteVector = std::vector<std::byte>;

struct XMFFile {
    std::string name;
    fs::path    source_file;

    Mesh mesh;

    // absolute path to bytes
    std::unordered_map<std::string, std::shared_ptr<ByteVector>>
        data_source_map;

    std::shared_ptr<ByteVector> fetch(std::string path) {
        auto iter = data_source_map.find(path);

        if (iter != data_source_map.end()) {
            std::cout << "File cached\n";
            return iter->second;
        }

        std::cout << "File uncached, loading...\n";

        std::ifstream stream(path);

        auto bv = std::make_shared<ByteVector>();

        bv->resize(fs::file_size(path));

        stream.read((char*)bv->data(), bv->size());

        data_source_map[path] = bv;

        return bv;
    }

    std::shared_ptr<ByteVector> get_path(std::string path) {
        std::cout << "Finding referenced file: " << path << std::endl;

        auto p = fs::path(path);

        if (fs::exists(p)) {
            // cool lets use it.
            return fetch(p.string());
        }

        // assume it is relative to us

        assert(source_file.is_absolute());

        auto parent = source_file.parent_path();

        for (auto const& dir_entry : fs::recursive_directory_iterator(parent)) {
            std::cout << "Search: " << dir_entry << '\n';
            if (dir_entry.path().filename() == p) {
                // this is it
                std::cout << "Matched.\n";
                return fetch(fs::absolute(dir_entry).string());
            }
        }

        // cant find it. bail

        std::cerr << "Unable to find referenced file: " << path << std::endl;
        exit(EXIT_FAILURE);
    }
};

template <class U, class T>
    requires(std::is_standard_layout_v<T> and std::is_standard_layout_v<U>)
std::span<U> span_cast(std::span<T> s) {
    static_assert(
        std::max(sizeof(T), sizeof(U)) % std::min(sizeof(T), sizeof(U)) == 0);

    return std::span<U>(reinterpret_cast<U*>(s.data()),
                        std::as_bytes(s).size() / sizeof(U));
}

template <class U, class T>
    requires(std::is_standard_layout_v<T> and std::is_standard_layout_v<U>)
std::span<U const> span_cast(std::span<T const> s) {
    static_assert(
        std::max(sizeof(T), sizeof(U)) % std::min(sizeof(T), sizeof(U)) == 0);

    return std::span<U const>(reinterpret_cast<U const*>(s.data()),
                              std::as_bytes(s).size() / sizeof(U));
}

template <class T>
auto safe_subspan(std::span<T> sp,
                  size_t       offset,
                  size_t       count = std::dynamic_extent) {
    if (offset >= sp.size()) return std::span<T>();

    if (count != std::dynamic_extent) {
        if ((offset + count) >= sp.size()) { return sp.subspan(offset); }
    }

    return sp.subspan(offset, count);
}

struct DataItem {
    std::string name;
    std::string format;
    int         precision = -1;
    std::string datatype;
    size_t      seek = 0;
    size_t      dims = 0;

    std::shared_ptr<ByteVector> data;

    std::span<std::byte> get_bytes() {
        if (!data) return {};

        auto ret = std::span<std::byte>(*data);

        if (ret.size() < seek) {
            std::cerr << "Unable to get data slice" << std::endl;
            exit(EXIT_FAILURE);
        }

        return safe_subspan(ret, seek);
    }

    template <class F>
    void read_as_int(F&& handler) {
        require(datatype == "Int" && (precision == -1 or precision == 4),
                "Unable to read as int list");
        auto arr = span_cast<int>(get_bytes());
        arr      = safe_subspan(arr, 0, dims);
        handler(arr);
    }

    template <class F>
    void read_as_double(F&& handler) {
        require(datatype == "Float" && precision == 8,
                "Unable to read as float list");
        auto arr = span_cast<double>(get_bytes());
        arr      = safe_subspan(arr, 0, dims);
        handler(arr);
    }
};

void print_tree(const pt::ptree& pt, int level) {
    const std::string sep(2 * level, ' ');
    for (const pt::ptree::value_type& v : pt) {
        std::cout << sep << v.first << " : " << v.second.data() << "\n";
        print_tree(v.second, level + 1);
    }
}

pt::ptree const* get_attribs(pt::ptree const& tree) {
    for (auto const& node : tree) {
        if (node.first == "<xmlattr>") { return &node.second; }
    }
    return nullptr;
}

DataItem import_data_item(pt::ptree const& data_tree, XMFFile& xmf) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    auto n = get_attribs(data_tree);

    if (!n) return {};

    auto maybe_name      = n->get_optional<std::string>("Name");
    auto format          = n->get<std::string>("Format");
    auto maybe_precision = n->get_optional<int8_t>("Precision");
    auto data_type       = n->get<std::string>("DataType");
    auto maybe_seek      = n->get_optional<uint64_t>("Seek");
    auto dimensions      = n->get<uint64_t>("Dimensions");

    return DataItem {
        .name      = maybe_name.value_or(""),
        .format    = format,
        .precision = maybe_precision.value_or(-1),
        .datatype  = data_type,
        .seek      = maybe_seek.value_or(0),
        .dims      = dimensions,
        .data      = xmf.get_path(data_tree.data()),
    };
}

void import_xmf_topo(pt::ptree const& tree, XMFFile& xmf) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    auto n = get_attribs(tree);

    if (!n) return;

    auto type = n->get_optional<std::string>("TopologyType")
                    .get_value_or(std::string {});
    auto count = n->get_optional<int64_t>("NumberOfElements").get_value_or(0);

    require(type == "Triangle", "Unsupported topology type");

    auto data = import_data_item(tree.get_child("DataItem"), xmf);

    data.read_as_int([&](std::span<int> source) {
        xmf.mesh.triangles.reserve(source.size() / 3);

        auto* p = source.data();
        while (p < (source.data() + source.size())) {
            xmf.mesh.triangles.push_back(
                openvdb::Vec3I(*p, *(p + 1), *(p + 2)));
            p += 3;
        }
    });
}

void import_xmf_geom(pt::ptree const& tree, XMFFile& xmf) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    auto n = get_attribs(tree);

    if (!n) return;

    auto type = n->get<std::string>("GeometryType");

    require(type == "XYZ", "Unsupported geometry type");

    auto data = import_data_item(tree.get_child("DataItem"), xmf);

    data.read_as_double([&](std::span<double> source) {
        static_assert(sizeof(openvdb::Vec3s) == (3 * 4));
        auto& ref = xmf.mesh.points;

        ref.reserve(source.size() / 3);

        auto* p = source.data();
        while (p < (source.data() + source.size())) {
            ref.push_back(openvdb::Vec3s(*p, *(p + 1), *(p + 2)));
            p += 3;
        }
    });
}

void import_xmf_attrib(pt::ptree const& tree, XMFFile& xmf) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
}

void import_xmf_time(pt::ptree const& tree, XMFFile& xmf) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    for (auto const& node : tree) {
        std::cerr << node.first << std::endl;
    }
}

void import_xmf_grid(pt::ptree const& gtree, XMFFile& xmf) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    if (auto n = get_attribs(gtree); n) {
        xmf.name = n->get_optional<std::string>("Name").value_or("Unknown");
    }
    for (auto const& node : gtree) {
        std::cerr << node.first << std::endl;
        if (node.first == "Time") { import_xmf_time(node.second, xmf); }
        if (node.first.starts_with("Topology")) {
            import_xmf_topo(node.second, xmf);
        }
        if (node.first.starts_with("Geometry")) {
            import_xmf_geom(node.second, xmf);
        }
        if (node.first.starts_with("Attribute")) {
            import_xmf_attrib(node.second, xmf);
        }
    }
}

XMFFile import_xmf_domain(pt::ptree const& dtree, fs::path const& path) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    XMFFile ret { .source_file = path };
    for (auto const& node : dtree) {
        std::cerr << node.first << std::endl;
        if (node.first == "Grid") { import_xmf_grid(node.second, ret); }
    }
    return ret;
}

std::vector<XMFFile> import_xmf_root(pt::ptree const& dtree,
                                     fs::path const&  path) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::vector<XMFFile> ret;

    for (auto const& node : dtree) {
        std::cerr << node.first << std::endl;
        if (node.first == "Domain") {
            ret.emplace_back(import_xmf_domain(node.second, path));
        }
    }

    return ret;
}


void burn_to_vdb(std::span<XMFFile> meshes, fs::path outpath) {
    openvdb::GridPtrVec to_write;

    for (auto const& f : meshes) {
        std::cout << "INFO " << f.mesh.points.size() << " "
                  << f.mesh.triangles.size() << std::endl;

        openvdb::math::Transform tf;

        // tf.preScale(10);

        //        auto grid = openvdb::tools::meshToVolume<openvdb::FloatGrid>(
        //            f.mesh, tf, 1, 1,
        //            openvdb::tools::UNSIGNED_DISTANCE_FIELD);

        auto grid =
            openvdb::tools::meshToVolume<openvdb::FloatGrid>(f.mesh, tf, 1, 1);

        grid->setName("surface");

        to_write.push_back(grid);
    }

    openvdb::io::File file(outpath);
    file.write(to_write);
    file.close();
}

int import_xmf(fs::path path, fs::path outpath) {
    if (!fs::is_regular_file(path)) {
        std::cout << "Unable to open file " << path << std::endl;
        return EXIT_FAILURE;
    }

    auto in_stream = std::ifstream(path);

    pt::ptree tree;

    pt::xml_parser::read_xml(
        in_stream, tree, boost::property_tree::xml_parser::trim_whitespace);

    print_tree(tree, 0);

    for (auto const& node : tree) {
        std::cerr << node.first << std::endl;
        if (node.first == "Xdmf") {
            auto xmf = import_xmf_root(node.second, path);
            burn_to_vdb(xmf, outpath);
        }
    }

    return EXIT_SUCCESS;
}

int mesh_to_volume(Arguments& c) {
    auto fname   = fs::path(toml::find<std::string>(c.root, "input"));
    auto outname = fs::path(toml::find<std::string>(c.root, "output"));

    if (fname.extension() == ".xmf" or fname.extension() == ".xdmf") {
        return import_xmf(fname, outname);
    }

    std::cout << "unable to handle file type " << fname.extension()
              << std::endl;

    return EXIT_FAILURE;
}
