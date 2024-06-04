#pragma once

#include <memory>

#include <AMReX_Vector.H>
#include <PltFileManager.H>

namespace sol {
class state;
}

namespace amrex {
class MultiFab;
class DistributionMapping;
} // namespace amrex

namespace amr_basic {

/// This class is an extension of PltFileManager, providing access to useful
/// members that are not exposed by default.
struct LocalPlotFile;
using LocalPlotFilePtr = std::shared_ptr<LocalPlotFile>;

amrex::Vector<amrex::MultiFab> const& get_multifabs(LocalPlotFilePtr const&);
amrex::MultiFab const& get_multifab_at(LocalPlotFilePtr const&, int nlvl);
amrex::DistributionMapping const& get_dmaps(LocalPlotFilePtr const& ptr,
                                            int                     nlvl);
std::vector<std::string> get_variable_list(LocalPlotFilePtr const& ptr);

pele::physics::pltfilemanager::PltFileManager&
get_manager(LocalPlotFilePtr& ptr);

void register_functions(sol::state&);

} // namespace amr_basic
