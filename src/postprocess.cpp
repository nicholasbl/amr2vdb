#include "postprocess.h"

#include <boost/algorithm/string/case_conv.hpp>

#include <openvdb/tools/Activate.h>
#include <openvdb/tools/Clip.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/Mask.h>


#include <span>

using FloatGrid    = openvdb::FloatGrid;
using FloatGridPtr = openvdb::FloatGrid::Ptr;
using Grid         = openvdb::GridBase;
using GridPtr      = openvdb::GridBase::Ptr;

template <typename OpType>
void process_typed_grid(openvdb::GridBase::Ptr grid, OpType&& op) {
// this is from:
// https://www.openvdb.org/documentation/doxygen/codeExamples.html#sGenericProg
#define CALL_OP(GridType) op(openvdb::gridPtrCast<GridType>(grid))
    if (grid->isType<openvdb::BoolGrid>()) CALL_OP(openvdb::BoolGrid);
    else if (grid->isType<openvdb::FloatGrid>())
        CALL_OP(openvdb::FloatGrid);
    else if (grid->isType<openvdb::DoubleGrid>())
        CALL_OP(openvdb::DoubleGrid);
    else if (grid->isType<openvdb::Int32Grid>())
        CALL_OP(openvdb::Int32Grid);
    else if (grid->isType<openvdb::Int64Grid>())
        CALL_OP(openvdb::Int64Grid);
    else if (grid->isType<openvdb::Vec3IGrid>())
        CALL_OP(openvdb::Vec3IGrid);
    else if (grid->isType<openvdb::Vec3SGrid>())
        CALL_OP(openvdb::Vec3SGrid);
    else if (grid->isType<openvdb::Vec3DGrid>())
        CALL_OP(openvdb::Vec3DGrid);
#undef CALL_OP
}

static inline std::string_view advance(std::string_view s) {
    return s.substr(1);
}

inline bool special_compare(std::string_view a, std::string_view b) {
    if (a.empty() or b.empty()) return false;
    return a.at(0) == b.at(0);
}

static bool match_simple_impl(std::string_view pattern,
                              std::string_view string) {
    while (!pattern.empty()) {
        char p = pattern.at(0);

        if (p == '?') {
            if (string.empty()) return false;
        } else if (p == '*') {
            if (match_simple_impl(advance(pattern), string)) { return true; }

            if (!string.empty() and
                match_simple_impl(pattern, advance(string))) {
                return true;
            }

            return false;
        } else {
            if (!special_compare(pattern, string)) return false;
        }

        pattern = advance(pattern);
        string  = advance(string);
    }

    return string.empty() and pattern.empty();
}


bool match_simple(std::string_view pattern, std::string_view string) {
    if (pattern.empty()) return false;
    return match_simple_impl(pattern, string);
}

GridPtr find_grid_by_name(std::string_view pattern, GridMap const& grids) {

    for (auto const& [k, v] : grids) {
        auto name = k;

        name = boost::algorithm::to_lower_copy(name);

        if (match_simple(pattern, name)) { return v; }
    }

    return nullptr;
}

GridPtr take_grid_by_name(std::string_view pattern, GridMap& grids) {

    for (auto iter = grids.begin(); iter != grids.end(); iter++) {
        auto name = iter->first;

        name = boost::algorithm::to_lower_copy(name);

        if (match_simple(pattern, name)) {
            auto ret = iter->second;
            grids.erase(iter);
            return ret;
        }
    }

    return nullptr;
}

openvdb::BoolGrid::Ptr
extract_mask(FloatGridPtr                       vfrac_grid,
             PostProcessOptions::PPVFrac const& vfrac_info) {

    auto copy_vfrac = openvdb::deepCopyTypedGrid<FloatGrid>(*vfrac_grid);

    // we use a roundabout way to select the parts we want. we deactivate the
    // parts that not desired

    constexpr float tolerance = 2.0;

    float at = vfrac_info.value;

    if (vfrac_info.keep_above_value) {
        at += tolerance / 2;
    } else {
        at -= tolerance / 2;
    }

    openvdb::tools::deactivate(*copy_vfrac, at, tolerance);

    // create mask from active values

    return openvdb::tools::interiorMask(*copy_vfrac);
}

template <class GridType>
auto cast_to(GridPtr ptr) {
    return openvdb::gridPtrCast<GridType>(ptr);
}

inline FloatGridPtr cast_to_float(GridPtr ptr) {
    return cast_to<FloatGrid>(ptr);
}

void trim_all_vfrac(PostProcessOptions::PPVFrac vfrac_info, GridMap& grids) {
    std::cout << "Trimming...\n";
    auto grid = find_grid_by_name(vfrac_info.name, grids);

    if (!grid) {
        std::cerr << "unable to find a compatible volume fraction grid"
                  << std::endl;
        return;
    }

    // to compute the mask, we first take a copy of the vfrac grid

    auto mask = extract_mask(cast_to<FloatGrid>(grid), vfrac_info);

    for (auto& [k, v] : grids) {
        if (v == grid) continue;

        process_typed_grid(v, [&w = v, &mask](auto gptr) {
            w = openvdb::tools::clip(*gptr, *mask);
        });
    }
}

static openvdb::Vec3fGrid::Ptr scalar_to_vector(FloatGridPtr&& grid) {
    auto op = [](openvdb::FloatGrid::ValueOnCIter const& iter, auto& accessor) {
        auto value = openvdb::Vec3f(iter);

        if (iter.isVoxelValue()) {
            accessor.setValue(iter.getCoord(), value);
        } else {
            openvdb::CoordBBox bbox;
            iter.getBoundingBox(bbox);
            accessor.getTree()->fill(bbox, value);
        }
    };

    auto ret = openvdb::Vec3fGrid::create();

    openvdb::tools::transformValues(grid->cbeginValueOn(), *ret, op);

    return ret;
}

void merge_velocity(PostProcessOptions::PPVelName const& names,
                    GridMap&                             grids) {
    std::cout << "Merging velocity...\n";
    auto vx = cast_to_float(take_grid_by_name(names.x, grids));
    auto vy = cast_to_float(take_grid_by_name(names.y, grids));
    auto vz = cast_to_float(take_grid_by_name(names.z, grids));

    if (!vx or !vy or !vz) return;

    auto upgrade_x = scalar_to_vector(std::move(vx));
    auto upgrade_y = scalar_to_vector(std::move(vy));

    auto combiner = [](openvdb::Vec3f const& a,
                       openvdb::Vec3f const& b,
                       openvdb::Vec3f&       result) {
        result = openvdb::Vec3f(a.x(), b.y(), 0);
    };

    std::cout << "Adding x + y...\n";
    upgrade_x->tree().combine(upgrade_y->tree(), combiner);

    // make sure y is done
    upgrade_y = nullptr;

    // y should be empty, x should have partial results. now get z

    auto upgrade_z = scalar_to_vector(std::move(vz));

    auto final_combiner = [](openvdb::Vec3f const& a,
                             openvdb::Vec3f const& b,
                             openvdb::Vec3f&       result) {
        result = openvdb::Vec3f(a.x(), a.y(), b.z());
    };

    std::cout << "Adding xy + z...\n";
    upgrade_x->tree().combine(upgrade_z->tree(), final_combiner);

    auto new_grid_name = "velocity";

    upgrade_x->setName(new_grid_name);

    grids[new_grid_name] = upgrade_x;
}

void add_in_magvort(GridMap& grids) {
    std::cout << "Computing magvort...\n";
    auto vel =
        cast_to<openvdb::Vec3fGrid>(find_grid_by_name("velocity", grids));

    if (!vel) return;

    auto vort = openvdb::tools::curl(*vel);

    auto magvort = openvdb::v10_0::tools::magnitude(*vort);

    auto new_grid_name = "magvort";

    magvort->setName(new_grid_name);

    grids[new_grid_name] = magvort;
}

void postprocess(PostProcessOptions opts, GridMap& l) {
    std::cout << "Starting postprocess\n";
    // first trim by vfrac

    if (opts.trim_vfrac.has_value()) {
        trim_all_vfrac(opts.trim_vfrac.value(), l);
    }

    // merge velocity if we need to
    if (opts.merge_velocity.has_value()) {
        merge_velocity(opts.merge_velocity.value(), l);
    }

    if (opts.compute_mag_vort) { add_in_magvort(l); }
}


PostProcessOptions PostProcessOptions::from_toml(toml::value const& root) {
    PostProcessOptions ret;

    auto node = toml::find(root, "post");

    if (node.contains("trim")) {
        auto trim_node = toml::find(node, "trim");

        bool gt = trim_node.contains("greater");
        bool lt = trim_node.contains("less");

        bool keep = true;

        // not sure about these defaults
        if (gt and lt) {
            keep = true;
        } else if (gt and !lt) {
            keep = true;
        } else if (!gt and lt) {
            keep = false;
        } else {
            keep = true;
        }

        PPVFrac trim_info {
            .name             = toml::find_or(trim_node, "variable", "v?frac"),
            .value            = toml::find_or<float>(trim_node, "value", .5),
            .keep_above_value = keep,
        };

        ret.trim_vfrac = trim_info;
    }

    if (node.contains("merge_velocity")) {
        auto merge_node = toml::find(node, "merge_velocity");

        PPVelName info;

        info.x = toml::find<std::string>(merge_node, "x");
        info.y = toml::find<std::string>(merge_node, "y");
        info.z = toml::find<std::string>(merge_node, "z");

        ret.merge_velocity = info;
    }

    if (node.contains("compute_magvort")) { ret.compute_mag_vort = true; }

    return ret;
}
