#include "argparse.h"

#include <iostream>

std::string take_first(std::deque<std::string>& args) {
    if (args.empty()) return {};
    auto ret = std::move(args.front());
    args.pop_front();
    return ret;
}

void merge_values(toml::value& old_value, toml::value& new_value) {
    switch (old_value.type()) {
    case toml::value_t::empty:
    case toml::value_t::boolean:
    case toml::value_t::integer:
    case toml::value_t::floating:
    case toml::value_t::string:
    case toml::value_t::offset_datetime:
    case toml::value_t::local_datetime:
    case toml::value_t::local_date:
    case toml::value_t::local_time:
    case toml::value_t::array: old_value = std::move(new_value); break;
    case toml::value_t::table:
        if (new_value.is_table()) {
            auto& a_as_table = old_value.as_table();
            auto& b_as_table = new_value.as_table();
            for (auto& [k, v] : b_as_table) {
                auto iter = a_as_table.find(k);

                if (iter == a_as_table.end()) {
                    a_as_table[k] = v;
                } else {
                    merge_values(a_as_table[k], v);
                }
            }
        }
        break;
    }
}


Arguments Arguments::parse(std::deque<std::string>&& args) {
    Arguments ret;

    // take the first one, that should be a file

    auto filename = take_first(args);

    if (filename.empty()) { return ret; }

    std::cout << "Reading configuration from " << filename << std::endl;

    auto config = toml::parse(filename);

    while (true) {

        auto next = take_first(args);

        if (next.empty()) break;

        std::cout << "Override " << next << std::endl;

        std::istringstream stream(next);

        auto override = toml::parse(stream);

        std::cout << "PARSED: " << override << std::endl;

        merge_values(config, override);
    }

    std::cout << "Done" << std::endl;

    ret.root = config;

    return ret;
}
