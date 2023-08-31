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

    std::vector<std::string> arguments = { argv + 1, argv + argc };

    if (arguments.empty()) return EXIT_FAILURE;

    Arguments args;

    for (auto& str : arguments) {

        if (str.starts_with("--")) {
            std::string_view view = str;
            view.remove_prefix(2);

            auto split_at = view.find('=');

            if (split_at == view.npos) {
                args.flags[std::string(view)] = {};
            } else {
                auto value = view.substr(split_at + 1);
                auto key   = view.substr(0, split_at);
                args.flags[std::string(key)].push_back(std::string(value));
            }
            continue;
        }

        args.positional.push_back(std::move(str));
    }

    for (auto const& [k, v] : args.flags) {
        std::string value;

        if (v.size()) {
            value = v.at(0);
            for (size_t i = 1; i < v.size(); i++) {
                value += ", " + v[i];
            }
        }


        std::cout << "Option " << k << ": " << value << std::endl;
    }

    auto command = args.take_first_positional();

    if (command == "amr") { return amr_to_volume(args); }
    if (command == "amr2") { return flatten_to_vdb(args); }
    if (command == "mesh") { return mesh_to_volume(args); }

    std::cerr << "unknown command " << command << std::endl;
}
