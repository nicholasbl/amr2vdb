#pragma once

#include <deque>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

#include <toml.hpp>

/// merge toml nodes. in the case of scalar or array, just replace old values
/// with new ones. if a table, merge the keys, replacing old keys with new keys
void merge_values(toml::value& old_value, toml::value& new_value);

struct Arguments {
    toml::value root;

    /// Consume command line arguments. The first argument is the toml file to
    /// consume. from that point on, further arguments are overrides
    static Arguments parse(std::deque<std::string>&& args);
};
