#include "amr.h"

#include "config.h"

#include <string>
#include <vector>
#include <iostream>

int main(int argc, char** argv) {
    std::vector<std::string> arguments = { argv, argv + argc };

    std::string target = "/Volumes/LocalStore/plt07400";
    if (arguments.size() > 1) { target = arguments.at(1); }

    auto amr = load_file(target, Config {});
}
