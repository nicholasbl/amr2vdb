#pragma once

#include <openvdb/openvdb.h>

#include <toml.hpp>

#include <optional>

using GridMap = std::unordered_map<std::string, openvdb::GridBase::Ptr>;

void postprocess(toml::value const& root, GridMap&);
