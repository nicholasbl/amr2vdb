#pragma once

#include <filesystem>
#include <memory>

#include <openvdb/openvdb.h>
#include <openvdb/tools/MultiResGrid.h>

struct Config;
struct AMRState;


using FloatMultiGrid    = openvdb::tools::MultiResGrid<openvdb::FloatTree>;
using FloatMultiGridPtr = openvdb::tools::MultiResGrid<openvdb::FloatTree>::Ptr;
using FloatGrid         = openvdb::FloatGrid;
using FloatGridPtr      = openvdb::FloatGrid::Ptr;

std::shared_ptr<AMRState> load_file(std::filesystem::path, Config const& c);


struct SampledGrid {
    std::string       name;
    FloatMultiGridPtr grid;
};

std::vector<SampledGrid> write_to_vdbs(AMRState const&);
