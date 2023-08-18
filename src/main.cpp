#include "argparse.h"
#include "mesh_to_volume.h"
#include "volume.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <boost/algorithm/string/replace.hpp>


int main(int argc, char** argv) {
    if (argc <= 1) return EXIT_FAILURE;

    std::vector<std::string> arguments = { argv + 1, argv + argc };

    if (arguments.empty()) return EXIT_FAILURE;

    Arguments args;

    bool        need_flag = false;
    std::string last_key;
    for (auto& str : arguments) {
        if (need_flag) {
            assert(last_key.size());
            args.flags[last_key].push_back(std::move(str));
            need_flag = false;
            continue;
        }

        if (str.starts_with("--")) {
            need_flag = true;
            boost::algorithm::replace_all(str, "-", "");
            last_key = std::move(str);
            continue;
        }

        args.positional.push_back(std::move(str));
    }

    auto command = args.take_first_positional();

    if (command.starts_with("amr")) { return amr_to_volume(args); }
    if (command.starts_with("mesh")) { return mesh_to_volume(args); }
}
