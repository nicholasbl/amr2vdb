#pragma once

#include <deque>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

#include <toml.hpp>

void merge_values(toml::value& old_value, toml::value& new_value);

struct Arguments {
    toml::value root;

    static Arguments parse(std::deque<std::string>&& args);
};
