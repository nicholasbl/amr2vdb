#include "argparse.h"
#include "mesh_to_volume.h"
#include "volume.h"
#include "volume2.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc <= 1) return EXIT_FAILURE;

    std::deque<std::string> arguments = { argv + 1, argv + argc };

    if (arguments.empty()) return EXIT_FAILURE;

    auto args = Arguments::parse(std::move(arguments));

    std::cout << "Configuration: " << args.root << std::endl;

    if (args.root.contains("amr")) { return amr_to_volume(args); }
    if (args.root.contains("flatten")) { return flatten_to_vdb(args); }
    if (args.root.contains("mesh")) { return mesh_to_volume(args); }

    std::cerr << "no command given" << std::endl;
}
