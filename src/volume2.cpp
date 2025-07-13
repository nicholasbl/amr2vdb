#include "volume2.h"

#include "amr_common.h"
#include "argparse.h"
#include "pltfilereader.h"

#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <spdlog/spdlog.h>

#include <optional>

namespace {
struct Result {
    // openvdb::BBoxd           bbox;
    openvdb::GridPtrVec grids;
};
} // namespace


static void loop_box(amrex::Box const&                       bx,
                     amrex::Array4<amrex::Real const> const& a,
                     openvdb::FloatGrid::Accessor&           accessor) {

    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    // std::cout << "BOX " << lo << hi << std::endl;

    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {

                float value = *a.ptr(i, j, k, 0);

                accessor.setValue({ i, j, k }, value);
            }
        }
    }
}

static FloatGridPtr extract_grid(PltFileReader& plt, int var_id) {
    auto geom = plt.getGeom(plt.getNlev() - 1);

    auto domain = geom.Domain();

    auto grid = amrex::BoxArray { domain };

    auto larray = domain.smallEnd().toArray();
    auto harray = domain.bigEnd().toArray();

    auto const dmap = amrex::DistributionMapping(grid);

    std::cout << "Domain L: " << larray[0] << " " << larray[1] << " "
              << larray[2] << std::endl;
    std::cout << "Domain H: " << harray[0] << " " << harray[1] << " "
              << harray[2] << std::endl;

    amrex::MultiFab output(grid, dmap, 1, 0);

    std::cout << "Interpolating..." << std::endl;
    // plt.fillPatchFromPlt(0, geom, var_id, 0, 1, output);
    spdlog::error("Functionality has been removed");

    std::cout << "Transform to VDB..." << std::endl;

    auto vdb_grid = FloatGrid::create();

    auto accessor = vdb_grid->getAccessor();

    auto iter = amrex::MFIter(output);

    for (; iter.isValid(); ++iter) {
        auto const& box = iter.validbox();

        amrex::FArrayBox const& fab = output[iter];

        auto const& a = fab.array();

        loop_box(box, a, accessor);
    }

    return vdb_grid;
}

static std::optional<Result> load_file(std::filesystem::path path,
                                       VolumeConfig          c) {
    auto state = AMRState();

    std::vector<std::string> var_names;
    std::vector<int>         var_ids;

    if (!std::filesystem::exists(path)) {
        std::cerr << "Path does not exist: " << path << std::endl;
        return std::nullopt;
    }

    auto plt_data = std::make_unique<PltFileReader>(path);
    auto plt_vars = plt_data->getVariableList();

    {
        std::set<std::string> extant_names(plt_vars.begin(), plt_vars.end());

        for (auto const& vname : c.variables) {
            if (!extant_names.contains(vname)) {
                std::cerr << "Variable " << vname
                          << " does not exist in pltfile.\n";
                return std::nullopt;
            }
        }
    }

    std::cout << "Loading from " << path << std::endl;

    for (long i = 0; i < plt_vars.size(); i++) {
        auto const& vname = plt_vars[i];
        if (c.variables.contains(vname)) {
            std::cout << "Variable " << vname << " at " << i << std::endl;
            var_names.push_back(vname);
            var_ids.push_back(i);
        }
    }

    if (c.max_level < 0) { c.max_level = plt_data->getNlev() - 1; }

    if (c.max_level >= plt_data->getNlev()) {
        std::cerr << "Asking for refinement level " << c.max_level
                  << " but pltfile only has " << plt_data->getNlev()
                  << std::endl;
        return std::nullopt;
    }

    plt_data->readPlotFileData();

    std::cout << "Done loading" << std::endl;

    Result res;

    for (int i = 0; i < var_names.size(); i++) {
        auto const& name = var_names.at(i);
        auto const& id   = var_ids.at(i);

        auto grid = extract_grid(*plt_data, id);

        grid->setName(name);

        res.grids.emplace_back(grid);
    }

    return res;
}

int flatten_to_vdb(Arguments const& c) {
    std::cout << "Using AMR flattening approach...\n";
    std::string  source_plt = toml::find<std::string>(c.root, "input");
    std::string  dest_vdb   = toml::find<std::string>(c.root, "output");
    VolumeConfig config;

    auto amr_config_node = toml::find(c.root, "flatten");

    auto source_vars = toml::find(amr_config_node, "variables").as_array();

    for (auto const& var : source_vars) {
        config.variables.insert(var.as_string());
    }

    auto amr = load_file(source_plt, config);

    if (!amr.has_value()) {
        std::cerr << "Unable to convert.\n";
        return EXIT_FAILURE;
    }

    openvdb::io::File file(dest_vdb);
    file.write(amr->grids);
    file.close();

    return EXIT_SUCCESS;
}
