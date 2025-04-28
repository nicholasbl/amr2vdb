#include "pltfilereader.h"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#include <filesystem>
#include <spdlog/spdlog.h>
#include <utility>

#include <H5Cpp.h>
using namespace amrex;


namespace {
const std::string level_prefix { "Level_" };
}

void GotoNextLine(std::istream& is) {
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

void debug_type(H5::DataType dtype) {
    std::cout << "HDF5 data type class: ";
    switch (dtype.getClass()) {
    case H5T_INTEGER: std::cout << "INTEGER" << std::endl; break;
    case H5T_FLOAT: std::cout << "FLOAT" << std::endl; break;
    case H5T_STRING: std::cout << "STRING" << std::endl; break;
    case H5T_COMPOUND: std::cout << "COMPOUND" << std::endl; break;
    default: std::cout << "OTHER" << std::endl; break;
    }
    if (dtype.getClass() == H5T_COMPOUND) {
        H5::CompType compType(dtype.getId());
        int          n_members = compType.getNmembers();
        std::cout << "Compound type has " << n_members << " fields:\n";
        for (int i = 0; i < n_members; ++i) {
            std::string  name        = compType.getMemberName(i);
            H5::DataType member_type = compType.getMemberDataType(i);
            std::cout << " - " << name << " (size: " << member_type.getSize()
                      << ")\n";
        }
    }
}

struct BoxEntry3D {
    int lo_i;
    int lo_j;
    int lo_k;
    int hi_i;
    int hi_j;
    int hi_k;
};


void PltFileReader::read_as_hdf5() {
    using namespace H5;

    H5File file(m_pltFile, H5F_ACC_RDONLY);

    // === Read number of components ===
    Attribute attr_ncomp = file.openAttribute("num_components");
    int       ncomp      = 0;
    attr_ncomp.read(PredType::NATIVE_INT, &ncomp);
    m_nvars = ncomp;

    spdlog::info("Found {} components", ncomp);

    // === Read variable names ===
    m_vars.resize(ncomp);
    for (int i = 0; i < ncomp; ++i) {
        std::string name        = "component_" + std::to_string(i);
        Attribute   attr        = file.openAttribute(name);
        StrType     strdatatype = attr.getStrType();
        attr.read(strdatatype, m_vars[i]);
        spdlog::info("Adding component {}", m_vars[i]);
    }

    // === Read number of levels ===
    Attribute attr_levels = file.openAttribute("num_levels");
    attr_levels.read(PredType::NATIVE_INT, &m_nlevels);

    m_geoms.resize(m_nlevels);
    m_grids.resize(m_nlevels);
    m_dmaps.resize(m_nlevels);
    m_data.resize(m_nlevels);

    spdlog::info("Found {} levels", m_nlevels);

    H5::CompType boxType(sizeof(BoxEntry3D));
    boxType.insertMember(
        "lo_i", HOFFSET(BoxEntry3D, lo_i), H5::PredType::NATIVE_INT);
    boxType.insertMember(
        "lo_j", HOFFSET(BoxEntry3D, lo_j), H5::PredType::NATIVE_INT);
    boxType.insertMember(
        "lo_k", HOFFSET(BoxEntry3D, lo_k), H5::PredType::NATIVE_INT);
    boxType.insertMember(
        "hi_i", HOFFSET(BoxEntry3D, hi_i), H5::PredType::NATIVE_INT);
    boxType.insertMember(
        "hi_j", HOFFSET(BoxEntry3D, hi_j), H5::PredType::NATIVE_INT);
    boxType.insertMember(
        "hi_k", HOFFSET(BoxEntry3D, hi_k), H5::PredType::NATIVE_INT);

    static_assert(AMREX_SPACEDIM == 3);

    for (int lev = 0; lev < m_nlevels; ++lev) {
        spdlog::info("Reading level {}", lev);
        std::string level_name  = "level_" + std::to_string(lev);
        Group       level_group = file.openGroup(level_name);

        // === Read boxes ===
        DataSet   boxes_ds = level_group.openDataSet("boxes");
        DataSpace space    = boxes_ds.getSpace();
        hsize_t   dims[1];
        space.getSimpleExtentDims(dims);
        const int nboxes = dims[0];

        std::vector<BoxEntry3D> box_data(nboxes);
        // debug_type(boxes_ds.getDataType());
        boxes_ds.read(box_data.data(), boxType);

        amrex::BoxList bl;
        for (auto const& box : box_data) {
            amrex::IntVect lo {
                box.lo_i,
                box.lo_j,
                box.lo_k,
            };

            amrex::IntVect hi {
                box.hi_i,
                box.hi_j,
                box.hi_k,
            };

            // for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            //     lo[d] = box_data[i * 2 * AMREX_SPACEDIM + d];
            //     hi[d] = box_data[i * 2 * AMREX_SPACEDIM + AMREX_SPACEDIM +
            //     d];
            // }
            bl.push_back(amrex::Box(lo, hi));
        }

        amrex::BoxArray ba(std::move(bl));
        m_grids[lev] = ba;

        // === Read dx, prob_lo, prob_hi ===
        std::vector<double> dx(AMREX_SPACEDIM);
        std::vector<double> plo(AMREX_SPACEDIM);
        std::vector<double> phi(AMREX_SPACEDIM);
        level_group.openAttribute("Vec_dx").read(PredType::NATIVE_DOUBLE,
                                                 dx.data());
        level_group.openAttribute("prob_lo").read(PredType::NATIVE_DOUBLE,
                                                  plo.data());
        level_group.openAttribute("prob_hi").read(PredType::NATIVE_DOUBLE,
                                                  phi.data());

        amrex::RealBox  rb(plo.data(), phi.data());
        amrex::Geometry geom(m_grids[lev].minimalBox(), &rb);
        m_geoms[lev] = geom;

        // === Create DistributionMapping and MultiFab ===
        amrex::DistributionMapping dm(m_grids[lev]);
        m_dmaps[lev] = dm;

        int nghost = 0;
        m_data[lev].define(m_grids[lev], m_dmaps[lev], m_nvars, nghost);

        // === Read actual data ===
        // Assume single dataset: "data:datatype=0"
        DataSet data_ds = level_group.openDataSet("data:datatype=0");
        hsize_t totalSize;
        data_ds.getSpace().getSimpleExtentDims(&totalSize);

        std::vector<amrex::Real> flat_data(totalSize);
        data_ds.read(flat_data.data(),
                     sizeof(amrex::Real) == 4 ? PredType::NATIVE_FLOAT
                                              : PredType::NATIVE_DOUBLE);

        // === Fill MultiFab ===
        size_t offset = 0;
        for (amrex::MFIter mfi(m_data[lev]); mfi.isValid(); ++mfi) {
            auto&     fab  = m_data[lev][mfi];
            const int npts = fab.box().numPts();
            std::memcpy(fab.dataPtr(),
                        &flat_data[offset],
                        npts * m_nvars * sizeof(amrex::Real));
            offset += npts * m_nvars;
        }
    }

    m_dataLoaded = true;
}

PltFileReader::PltFileReader(std::string a_pltFile)
    : m_pltFile { std::move(a_pltFile) } {

    // this breaks
    // if (!hdf5::file::is_hdf5_file(a_pltFile)) {
    if (!std::filesystem::is_directory(a_pltFile)) {
        // try for HDF5
        read_as_hdf5();
        return;
    }

    // Get the plt metadata and resize part of the data vectors
    const std::string pltFileHeader(m_pltFile + "/Header");
    readGenericPlotfileHeader(pltFileHeader);

    // Resize the actual data container
    m_grids.resize(m_nlevels);
    m_dmaps.resize(m_nlevels);
    m_data.resize(m_nlevels);

    // Read the pltfile metadata only
    readPlotFileMetaData();
}

void PltFileReader::readGenericPlotfileHeader(
    const std::string& a_pltFileHeader) {
    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(a_pltFileHeader, fileCharPtr);
    std::string        fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // Title line
    std::getline(is, line);

    // Number of variables
    m_nvars = 0;
    is >> m_nvars;
    GotoNextLine(is);

    // Extract variables names
    m_vars.resize(m_nvars);
    for (int n = 0; n < m_nvars; n++) {
        is >> m_vars[n];
        GotoNextLine(is);
    }

    // Get and check space dimension
    int PLT_SPACEDIM = AMREX_SPACEDIM;
    is >> PLT_SPACEDIM;
    GotoNextLine(is);
    AMREX_ASSERT(PLT_SPACEDIM == AMREX_SPACEDIM);

    // Simulation time
    is >> m_time;
    GotoNextLine(is);

    // Number of levels
    is >> m_nlevels;
    GotoNextLine(is);
    m_nlevels += 1; // Finest is stored, need to add 1

    // Setup geometry data
    m_geoms.resize(m_nlevels);
    m_refRatio.resize(m_nlevels - 1);

    // Level 0 geometry
    Real prob_lo[AMREX_SPACEDIM];
    Real prob_hi[AMREX_SPACEDIM];
    // Low coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int                i = 0;
        while (lis >> word) {
            prob_lo[i++] = std::stod(word);
        }
    }

    // High coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int                i = 0;
        while (lis >> word) {
            prob_hi[i++] = std::stod(word);
        }
    }

    // Set up PltFile domain real box
    RealBox rb(prob_lo, prob_hi);

    std::getline(is, line);
    {
        std::istringstream lis(line);
        int                i = 0;
        while (lis >> word) {
            m_refRatio[i++] = std::stoi(word);
        }
    }

    // Get levels Domains
    Vector<Box> Domains(m_nlevels);
    for (int lev = 0; lev < m_nlevels; ++lev) {
        is >> Domains[lev];
    }
    GotoNextLine(is);
    GotoNextLine(is); // Skip nsteps line
    for (int lev = 0; lev < m_nlevels; ++lev) {
        GotoNextLine(is); // Skip dx line
    }

    // Coordinate system
    int coord_sys = 0;
    is >> coord_sys;
    GotoNextLine(is);

    // Populate the geometry vector, assume no periodicity
    Array<int, AMREX_SPACEDIM> perio({ AMREX_D_DECL(0, 0, 0) });
    m_geoms[0] = Geometry(Domains[0], rb, coord_sys, perio);
    for (int lev = 1; lev < m_nlevels; ++lev) {
        m_geoms[lev] = refine(m_geoms[lev - 1], m_refRatio[lev - 1]);
    }
}

void PltFileReader::readPlotFileMetaData() {
    // Set BoxArray, DistMap on each level
    for (int lev = 0; lev < m_nlevels; ++lev) {
        readLevelBoxArray(lev, m_grids[lev]);
        m_dmaps[lev] = DistributionMapping(m_grids[lev]);
    }
}

void PltFileReader::readPlotFileData() {
    if (m_dataLoaded) return;
    // Load the actual data
    // TODO: only load a subset of the pltfile variables
    for (int lev = 0; lev < m_nlevels; ++lev) {
        m_data[lev].define(m_grids[lev], m_dmaps[lev], m_nvars, 0);
        VisMF::Read(
            m_data[lev],
            MultiFabFileFullPrefix(lev, m_pltFile, level_prefix, "Cell"));
    }

    m_dataLoaded = true;
}

void PltFileReader::readLevelBoxArray(int a_lev, BoxArray& a_grid) {
    const std::string lvlHeader(m_pltFile + "/" + level_prefix +
                                std::to_string(a_lev) + "/Cell_H");

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(lvlHeader, fileCharPtr);
    std::string        fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line;

    std::getline(is, line); // Skip
    GotoNextLine(is);       // Skip
    GotoNextLine(is);       // Skip
    GotoNextLine(is);       // Skip
    a_grid.readFrom(is);
}
