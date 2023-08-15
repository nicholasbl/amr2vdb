#pragma once

#include <set>

struct Config {
    int                   max_level = -1;
    std::set<std::string> variables;
};
