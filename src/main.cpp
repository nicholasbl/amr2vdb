#include "amr.h"

#include "config.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    std::vector<std::string> arguments = { argv, argv + argc };

    std::string home = std::getenv("HOME");

    std::string target = home + "/Downloads/plt07400";
    if (arguments.size() > 1) { target = arguments.at(1); }

    auto amr = load_file(target, Config {});

    write_to_vdbs(*amr);
}
