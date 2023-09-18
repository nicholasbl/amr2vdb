#pragma once

#include <openvdb/openvdb.h>

#include <toml.hpp>

using GridMap = std::unordered_map<std::string, openvdb::GridBase::Ptr>;

/// Run postprocessing on this grid collection
void postprocess(toml::value const& root, GridMap&);
