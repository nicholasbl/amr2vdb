#pragma once

#include <PltFileManager.H>

#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/MultiResGrid.h>

#include <memory>

using FloatMultiGrid    = openvdb::tools::MultiResGrid<openvdb::FloatTree>;
using FloatMultiGridPtr = openvdb::tools::MultiResGrid<openvdb::FloatTree>::Ptr;
using FloatGrid         = openvdb::FloatGrid;
using FloatGridPtr      = openvdb::FloatGrid::Ptr;


struct AMRState {
    amrex::AMReX* amr;

    AMRState() {
        // dummy params
        std::string          exe_name = "amr2vdb";
        std::array<char*, 1> argv_src = { exe_name.data() };
        int                  argc     = argv_src.size();
        char**               argv     = argv_src.data();

        amr = amrex::Initialize(argc, argv);
    }

    AMRState(AMRState const&)                  = delete;
    AMRState const& operator=(AMRState const&) = delete;
    AMRState(AMRState&&)                       = delete;
    AMRState const& operator=(AMRState&&)      = delete;

    ~AMRState() { amrex::Finalize(amr); }
};

struct LocalPlotFile : public pele::physics::pltfilemanager::PltFileManager {
public:
    LocalPlotFile(std::filesystem::path path)
        : pele::physics::pltfilemanager::PltFileManager(path) { }

    ~LocalPlotFile() = default;

    LocalPlotFile(LocalPlotFile const&)                  = delete;
    LocalPlotFile const& operator=(LocalPlotFile const&) = delete;
    LocalPlotFile(LocalPlotFile&&)                       = delete;
    LocalPlotFile const& operator=(LocalPlotFile&&)      = delete;

    auto const& all_data() const { return m_data; }

    auto const& get_data(int nlvl) { return m_data.at(nlvl); }

    auto const& get_dmap(int nlvl) { return m_dmaps.at(nlvl); }
};


struct VolumeConfig {
    int                   max_level = -1;
    std::set<std::string> variables;
};
