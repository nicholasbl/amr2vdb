#pragma once

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

// Shamelessly stolen from Pele
// https://github.com/AMReX-Combustion/PelePhysics
class PltFileReader {
    // Plot metadata
    std::string                m_pltFile;     // pltfile name
    int                        m_nvars { 0 }; // number of variables
    amrex::Vector<std::string> m_vars; // list of variable in the plotfile
    int                        m_nlevels { 0 }; // number of levels
    amrex::Real                m_time { 0.0 };  // Simulation time
    bool m_dataLoaded { false }; // Flag to check is the data has been read in

    // Geometry, grid and data containers
    amrex::Vector<amrex::BoxArray>            m_grids;
    amrex::Vector<amrex::Geometry>            m_geoms;
    amrex::Vector<int>                        m_refRatio;
    amrex::Vector<amrex::DistributionMapping> m_dmaps;
    amrex::Vector<amrex::MultiFab>            m_data;

    void read_as_hdf5();

public:
    explicit PltFileReader(std::string a_pltFile);
    ~PltFileReader() = default;

    PltFileReader(PltFileReader const&)                  = delete;
    PltFileReader const& operator=(PltFileReader const&) = delete;
    PltFileReader(PltFileReader&&)                       = delete;
    PltFileReader const& operator=(PltFileReader&&)      = delete;

    void readGenericPlotfileHeader(const std::string& a_pltFileHeader);

    void readPlotFileMetaData();

    void readPlotFileData();

    void readLevelBoxArray(int a_lev, amrex::BoxArray& a_grid);

    amrex::Vector<std::string> getVariableList() { return m_vars; }
    int                        getNlev() const { return m_nlevels; }
    amrex::Real                getTime() const { return m_time; }
    const amrex::BoxArray&     getGrid(int a_lev) { return m_grids.at(a_lev); }
    const amrex::Geometry&     getGeom(int a_lev) { return m_geoms.at(a_lev); }


    /// Obtain the list of multifabs, to be indexed by refinement level
    auto const& all_data() const { return m_data; }

    /// Obtain the multifab at a refinement level
    auto const& get_data(int nlvl) { return m_data.at(nlvl); }

    /// Obtain the distribution mapping at a given refinement level
    auto const& get_dmap(int nlvl) { return m_dmaps.at(nlvl); }
};
