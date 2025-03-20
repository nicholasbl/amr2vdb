#include "argparse.h"
#include "export.h"
#include "mesh_to_volume.h"
#include "points.h"
#include "special.h"
#include "volume.h"
#include "volume2.h"

#include <cstdlib>
#include <string>

#include <openvdb/openvdb.h>

#include "spdlog/spdlog.h"

int main(int argc, char** argv) {
    if (argc <= 1) return EXIT_FAILURE;

    openvdb::initialize();

    std::deque<std::string> arguments = { argv + 1, argv + argc };

    if (arguments.empty()) return EXIT_FAILURE;

    auto args = Arguments::parse(std::move(arguments));

    auto cfg = toml::format(args.root);

    spdlog::info("Configuration: {}", cfg);

    if (args.root.contains("debug")) {
        spdlog::info("Enable debug output");
        spdlog::set_level(spdlog::level::debug);
    } else {
        spdlog::set_level(spdlog::level::info);
    }

    if (args.root.contains("amr")) { return amr_to_volume(args); }
    if (args.root.contains("flatten")) { return flatten_to_vdb(args); }
    if (args.root.contains("mesh")) { return mesh_to_volume(args); }
    if (args.root.contains("all_iso_merge")) { return all_iso_merge(args); }
    if (args.root.contains("points")) { return txt_to_points(args); }
    if (args.root.contains("export")) { return export_quantity(args); }

    spdlog::error("no command given");
}
