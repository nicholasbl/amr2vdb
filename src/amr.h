#pragma once

#include <filesystem>
#include <memory>

#include <openvdb/openvdb.h>
#include <openvdb/tools/MultiResGrid.h>

struct Config;
struct AMRState;


using FloatMultiGrid = openvdb::tools::MultiResGrid<openvdb::v10_0::FloatTree>;

std::shared_ptr<AMRState> load_file(std::filesystem::path, Config const& c);

void write_to_vdbs(AMRState const&);
