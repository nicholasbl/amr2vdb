#pragma once

#include <PltFileManager.H>

#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/MultiResGrid.h>

#include <filesystem>
#include <memory>
#include <optional>

using FloatMultiGrid    = openvdb::tools::MultiResGrid<openvdb::FloatTree>;
using FloatMultiGridPtr = openvdb::tools::MultiResGrid<openvdb::FloatTree>::Ptr;
using FloatGrid         = openvdb::FloatGrid;
using FloatGridPtr      = openvdb::FloatGrid::Ptr;

/// Manages an AMR state, automatically destroys the state when it goes out of
/// scope
struct AMRState {
    amrex::AMReX* amr;

    AMRState() {
        // dummy params. we do this to avoid using the ParamParse tools in
        // AMReX.
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

/// This class is an extension of PltFileManager, providing access to useful
/// members that are not exposed by default.
struct LocalPlotFile : public pele::physics::pltfilemanager::PltFileManager {
public:
    LocalPlotFile(std::filesystem::path path)
        : pele::physics::pltfilemanager::PltFileManager(path) { }

    ~LocalPlotFile() = default;

    LocalPlotFile(LocalPlotFile const&)                  = delete;
    LocalPlotFile const& operator=(LocalPlotFile const&) = delete;
    LocalPlotFile(LocalPlotFile&&)                       = delete;
    LocalPlotFile const& operator=(LocalPlotFile&&)      = delete;

    /// Obtain the list of multifabs, to be indexed by refinement level
    auto const& all_data() const { return m_data; }

    /// Obtain the multifab at a refinement level
    auto const& get_data(int nlvl) { return m_data.at(nlvl); }

    /// Obtain the distribution mapping at a given refinement level
    auto const& get_dmap(int nlvl) { return m_dmaps.at(nlvl); }
};

/// Configuration for a volume extraction
struct VolumeConfig {
    /// Max level of refinement to read
    int                        max_level = -1;
    /// Variable names to extract
    std::set<std::string>      variables;
    /// Optionally, save the AMR structure in the extraction with the given name
    std::optional<std::string> save_amr;
};
