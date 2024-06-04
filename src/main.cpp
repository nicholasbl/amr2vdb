#include "argparse.h"
#include "mesh_to_volume.h"
#include "special.h"
#include "volume.h"
#include "volume2.h"

#include <cstdlib>
#include <string>

#include <openvdb/openvdb.h>

#include "use_sol.h"

#include "spdlog/spdlog.h"

int main(int argc, char** argv) {
    if (argc <= 1) return EXIT_FAILURE;

    openvdb::initialize();

    sol::state lua;

    lua.open_libraries(sol::lib::base, sol::lib::package, sol::lib::io);

    std::deque<std::string> arguments = { argv + 1, argv + argc };

    if (arguments.empty()) return EXIT_FAILURE;

    std::string script_name = arguments.at(0);

    arguments.pop_front();

    bool use_debug = false;

    auto debug_env = std::getenv("DEBUG");

    if (debug_env) { use_debug = std::atoi(debug_env); }

    if (use_debug) {
        spdlog::info("Enable debug output");
        spdlog::set_level(spdlog::level::debug);
    } else {
        spdlog::set_level(spdlog::level::info);
    }

    lua["args"] = arguments;

    lua.script_file(script_name);


    // if (args.root.contains("amr")) { return amr_to_volume(args); }
    // if (args.root.contains("flatten")) { return flatten_to_vdb(args); }
    // if (args.root.contains("mesh")) { return mesh_to_volume(args); }
    // if (args.root.contains("all_iso_merge")) { return all_iso_merge(args); }

    // spdlog::error("no command given");
}
