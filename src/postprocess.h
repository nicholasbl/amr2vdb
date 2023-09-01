#pragma once

#include <openvdb/openvdb.h>

#include <toml.hpp>

#include <optional>


struct PostProcessOptions {
    struct PPVFrac {
        std::string name             = "v?frac";
        float       value            = .5;
        bool        keep_above_value = true;
    };

    struct PPVelName {
        std::string x, y, z;
    };

    std::optional<PPVFrac>   trim_vfrac;
    std::optional<PPVelName> merge_velocity;

    bool compute_mag_vort = false;

    static PostProcessOptions from_toml(toml::value const& root);
};

using GridMap = std::unordered_map<std::string, openvdb::GridBase::Ptr>;

void postprocess(PostProcessOptions, GridMap&);
