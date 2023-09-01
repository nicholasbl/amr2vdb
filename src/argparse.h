#pragma once

#include <deque>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

#include <toml.hpp>

struct Arguments {
    toml::value root;

    static Arguments parse(std::deque<std::string>&& args);
};
