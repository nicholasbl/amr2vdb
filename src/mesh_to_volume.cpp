#include "mesh_to_volume.h"

#include "argparse.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <filesystem>

namespace fs = std::filesystem;
namespace pt = boost::property_tree;

using openvdb::Vec3I;
using openvdb::Vec3s;
using openvdb::Vec4I;

struct Mesh {
    std::vector<Vec3s> points;
    std::vector<Vec3I> triangles;
    std::vector<Vec4I> quads;
};

struct XMFFile {
    std::string name;
    fs::path    source_file;

    std::unordered_map<fs::path, std::vector<std::byte>> data_source_map;
};

struct DataItem { };

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

DataItem import_data_item(pt::ptree const& data_tree, XMFFile& xmf) { }

void import_xmf_topo(pt::ptree const& tree, XMFFile& xmf) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    auto type  = std::string {};
    auto count = int64_t {};

    if (auto n = get_attribs(tree); n) {
        type = n->get_optional<std::string>("TopologyType")
                   .get_value_or(std::string {});
        count = n->get_optional<int64_t>("NumberOfElements").get_value_or(0);
    }

    auto data = import_data_item(tree.get_child("DataItem"), xmf);
}

void import_xmf_geom(pt::ptree const& tree, XMFFile& xmf) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
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

void import_xmf_root(pt::ptree const& dtree, fs::path const& path) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    for (auto const& node : dtree) {
        std::cerr << node.first << std::endl;
        if (node.first == "Domain") { import_xmf_domain(node.second, path); }
    }
}

int import_xmf(fs::path path) {
    if (!fs::is_regular_file(path)) {
        std::cout << "Unable to open file " << path << std::endl;
        return EXIT_FAILURE;
    }

    auto in_stream = std::ifstream(path);

    pt::ptree tree;

    pt::xml_parser::read_xml(in_stream, tree);

    print_tree(tree, 0);

    for (auto const& node : tree) {
        std::cerr << node.first << std::endl;
        if (node.first == "Xdmf") { import_xmf_root(node.second, path); }
    }

    return EXIT_SUCCESS;
}

int mesh_to_volume(Arguments& c) {
    auto fname = fs::path(c.take_first_positional());

    if (fname.extension() == ".xmf" or fname.extension() == ".xdmf") {
        return import_xmf(fname);
    }

    std::cout << "unable to handle file type " << fname.extension()
              << std::endl;

    return EXIT_FAILURE;
}
