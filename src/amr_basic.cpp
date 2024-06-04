#include "amr_basic.h"

#include "use_sol.h"

#include "spdlog/spdlog.h"


namespace amr_basic {

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

amrex::Vector<amrex::MultiFab> const&
get_multifabs(LocalPlotFilePtr const& ptr) {
    return ptr->all_data();
}
amrex::MultiFab const& get_multifab_at(LocalPlotFilePtr const& ptr, int nlvl) {
    return ptr->get_data(nlvl);
}
amrex::DistributionMapping const& get_dmaps(LocalPlotFilePtr const& ptr,
                                            int                     nlvl) {
    return ptr->get_dmap(nlvl);
}

pele::physics::pltfilemanager::PltFileManager&
get_manager(LocalPlotFilePtr& ptr) {
    return *ptr;
}

void register_functions(sol::state& state) {
    state.set_function("load_plt", [](std::string path) {
        spdlog::info("Loading plt from {}", path);
        return std::make_shared<LocalPlotFile>(path);
    });
}

} // namespace amr_basic
