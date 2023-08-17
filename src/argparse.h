#pragma once

#include <deque>
#include <string>
#include <unordered_map>

struct Arguments {
    std::deque<std::string>                                   positional;
    std::unordered_map<std::string, std::vector<std::string>> flags;

    std::string take_first_positional() {
        if (positional.empty()) return {};
        auto ret = positional.front();
        positional.pop_front();
        return ret;
    }
};
