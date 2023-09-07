#include "postprocess.h"

#include <boost/algorithm/string/case_conv.hpp>

#include <openvdb/tools/Activate.h>
#include <openvdb/tools/Clip.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/Mask.h>

#include "spdlog/spdlog.h"

#include <span>

using FloatGrid    = openvdb::FloatGrid;
using FloatGridPtr = openvdb::FloatGrid::Ptr;
using Grid         = openvdb::GridBase;
using GridPtr      = openvdb::GridBase::Ptr;

// =============================================================================

struct PPVFrac {
    std::string name             = "v?frac";
    float       value            = .5;
    bool        keep_above_value = true;
};

struct PPVelName {
    std::string x, y, z;
};

// =============================================================================

#define BIN_OP(OPER)                                                           \
    inline LVec operator OPER(LVec b) const {                                  \
        LVec ret;                                                              \
        for (int i = 0; i < this->size(); i++) {                               \
            ret[i] = (*this)[i] OPER b[i];                                     \
        }                                                                      \
        return ret;                                                            \
    }

#define BIN_OP_SCALAR(OPER)                                                    \
    inline LVec operator OPER(float b) const {                                 \
        LVec ret;                                                              \
        for (int i = 0; i < this->size(); i++) {                               \
            ret[i] = (*this)[i] OPER b;                                        \
        }                                                                      \
        return ret;                                                            \
    }

template <int N>
struct LVec : std::array<float, N> {
    BIN_OP(-)
    BIN_OP(+)
    BIN_OP_SCALAR(-)
    BIN_OP_SCALAR(+)
};

template <int N>
inline std::ostream& operator<<(std::ostream& s, LVec<N> const&) {
    s << "[]";
    return s;
}

using Vec6 = LVec<6>;
using Mat3 = LVec<9>;

namespace openvdb {
namespace OPENVDB_VERSION_NAME {
template <>
inline Vec6 zeroVal<Vec6>() {
    return {};
}

template <>
inline Mat3 zeroVal<Mat3>() {
    return {};
}

namespace math {
template <>
inline Vec6 negative(Vec6 const& val) {
    Vec6 ret;
    for (auto& i : ret) {
        i = -i;
    }
    return ret;
}
template <>
inline ::Mat3 negative(::Mat3 const& val) {
    ::Mat3 ret;
    for (auto& i : ret) {
        i = -i;
    }
    return ret;
}


} // namespace math

} // namespace OPENVDB_VERSION_NAME
} // namespace openvdb

inline Vec6 Abs(Vec6 v) {
    for (auto& i : v) {
        i = openvdb::math::Abs(i);
    }
    return v;
}

inline ::Mat3 Abs(::Mat3 v) {
    for (auto& i : v) {
        i = openvdb::math::Abs(i);
    }
    return v;
}

using Vec6Tree = openvdb::tree::Tree4<Vec6, 5, 4, 3>::Type;
using Mat3Tree = openvdb::tree::Tree4<Mat3, 5, 4, 3>::Type;

using Vec6Grid = openvdb::Grid<Vec6Tree>;
using Mat3Grid = openvdb::Grid<Mat3Tree>;

// =============================================================================

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

// =============================================================================

openvdb::BoolGrid::Ptr extract_mask(FloatGridPtr   vfrac_grid,
                                    PPVFrac const& vfrac_info) {
    spdlog::info("Extracting mask...");

    auto ret = openvdb::BoolGrid::create();

    auto xfrmr_keep_lt = [&](FloatGrid::ValueOnCIter const& iter,
                             openvdb::BoolGrid::Accessor&   accessor) {
        auto value = *iter;

        auto out_value = (value <= vfrac_info.value);

        if (iter.isVoxelValue()) {
            accessor.setValue(iter.getCoord(), out_value);
        } else {
            openvdb::CoordBBox bbox;
            iter.getBoundingBox(bbox);
            accessor.getTree()->fill(bbox, out_value);
        }
    };

    auto xfrmr_keep_gt = [&](FloatGrid::ValueOnCIter const& iter,
                             openvdb::BoolGrid::Accessor&   accessor) {
        auto value = *iter;

        auto out_value = (value >= vfrac_info.value);

        if (iter.isVoxelValue()) {
            accessor.setValue(iter.getCoord(), out_value);
        } else {
            openvdb::CoordBBox bbox;
            iter.getBoundingBox(bbox);
            accessor.getTree()->fill(bbox, out_value);
        }
    };

    if (vfrac_info.keep_above_value) {
        openvdb::tools::transformValues(
            vfrac_grid->cbeginValueOn(), *ret, xfrmr_keep_gt);
    } else {
        openvdb::tools::transformValues(
            vfrac_grid->cbeginValueOn(), *ret, xfrmr_keep_lt);
    }

    openvdb::tools::deactivate(*ret, false);

    openvdb::tools::pruneInactive(ret->tree());

    // create mask from active values

    return ret;
}

// =============================================================================

template <class GridType>
auto cast_to(GridPtr ptr) {
    return openvdb::gridPtrCast<GridType>(ptr);
}

inline FloatGridPtr cast_to_float(GridPtr ptr) {
    return cast_to<FloatGrid>(ptr);
}

// =============================================================================

void trim_all_vfrac(PPVFrac vfrac_info, GridMap& grids) {
    spdlog::info("Trimming...");
    auto grid = find_grid_by_name(vfrac_info.name, grids);

    if (!grid) {
        spdlog::error("unable to find a compatible volume fraction grid");
        return;
    }

    spdlog::info("Trim will use {}", grid->getName());

    // to compute the mask, we first take a copy of the vfrac grid

    auto mask = extract_mask(cast_to<FloatGrid>(grid), vfrac_info);

    for (auto& [k, v] : grids) {
        if (v == grid) continue;

        spdlog::info("Trimming {}", v->getName());

        process_typed_grid(v, [&w = v, &mask](auto gptr) {
            w = openvdb::tools::clip(*gptr, *mask);
        });
    }
}

// =============================================================================

template <class OutGrid, class InGrid, class Function>
auto convert(InGrid&& in, Function&& tf) {
    auto op = [&](auto const& iter, auto& accessor) {
        auto value = tf(*iter);

        if (iter.isVoxelValue()) {
            accessor.setValue(iter.getCoord(), value);
        } else {
            openvdb::CoordBBox bbox;
            iter.getBoundingBox(bbox);
            accessor.getTree()->fill(bbox, value);
        }
    };

    auto ret = OutGrid::create();

    openvdb::tools::transformValues(in->cbeginValueOn(), *ret, op);

    return ret;
}

using Vec2fGrid = openvdb::Grid<openvdb::Vec2STree>;

static Vec2fGrid::Ptr scalar_to_vector2(FloatGridPtr&& grid) {
    auto op = [](openvdb::FloatGrid::ValueOnCIter const& iter, auto& accessor) {
        auto value = openvdb::Vec2f(*iter);

        if (iter.isVoxelValue()) {
            accessor.setValue(iter.getCoord(), value);
        } else {
            openvdb::CoordBBox bbox;
            iter.getBoundingBox(bbox);
            accessor.getTree()->fill(bbox, value);
        }
    };

    auto ret = Vec2fGrid::create();

    openvdb::tools::transformValues(grid->cbeginValueOn(), *ret, op);

    return ret;
}

static openvdb::Vec3fGrid::Ptr scalar_to_vector3(FloatGridPtr&& grid) {
    auto op = [](openvdb::FloatGrid::ValueOnCIter const& iter, auto& accessor) {
        auto value = openvdb::Vec3f(*iter);

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

void merge_velocity(PPVelName const& names, GridMap& grids) {
    spdlog::info("Merging velocity...");
    auto vx = cast_to_float(take_grid_by_name(names.x, grids));
    auto vy = cast_to_float(take_grid_by_name(names.y, grids));
    auto vz = cast_to_float(take_grid_by_name(names.z, grids));

    if (!vx or !vy or !vz) return;

    auto upgrade_x = scalar_to_vector3(std::move(vx));
    auto upgrade_y = scalar_to_vector3(std::move(vy));

    auto combiner = [](openvdb::Vec3f const& a,
                       openvdb::Vec3f const& b,
                       openvdb::Vec3f&       result) {
        result = openvdb::Vec3f(a.x(), b.x(), 1);
    };

    spdlog::debug("Adding x + y...");
    upgrade_x->tree().combine(upgrade_y->tree(), combiner);

    // make sure y is done
    upgrade_y = nullptr;

    // y should be empty, x should have partial results. now get z

    auto upgrade_z = scalar_to_vector3(std::move(vz));

    auto final_combiner = [](openvdb::Vec3f const& xy,
                             openvdb::Vec3f const& z,
                             openvdb::Vec3f&       result) {
        result = openvdb::Vec3f(xy.x(), xy.y(), z.x());
    };

    spdlog::debug("Adding xy + z...");
    upgrade_x->tree().combine(upgrade_z->tree(), final_combiner);

    auto new_grid_name = "velocity";

    upgrade_x->setName(new_grid_name);
    upgrade_x->setVectorType(openvdb::VecType::VEC_CONTRAVARIANT_RELATIVE);

    grids[new_grid_name] = upgrade_x;
}

// =============================================================================

void add_in_magvort(GridMap& grids) {
    spdlog::info("Computing magvort...");
    auto vel =
        cast_to<openvdb::Vec3fGrid>(find_grid_by_name("velocity", grids));

    if (!vel) return;

    auto vort = openvdb::tools::curl(*vel);

    auto magvort = openvdb::v10_0::tools::magnitude(*vort);

    auto new_grid_name = "magvort";

    magvort->setName(new_grid_name);
    magvort->setGridClass(openvdb::GridClass::GRID_FOG_VOLUME);

    grids[new_grid_name] = magvort;
}

// =============================================================================

// Todo: generalize
void add_in_mag_grad_density(std::string density_name, GridMap& grids) {
    spdlog::info("Computing mag(grad({}))...", density_name);

    auto density =
        cast_to<openvdb::FloatGrid>(find_grid_by_name(density_name, grids));

    if (!density) {
        spdlog::error("Unable to find density grid");
        return;
    }

    auto g  = openvdb::tools::gradient(*density);
    density = nullptr; // clear out mem

    auto m = openvdb::tools::magnitude(*g);
    g      = nullptr; // clear out mem, again

    auto new_grid_name = "mag_grad_" + density_name;

    m->setName(new_grid_name);
    m->setGridClass(openvdb::GridClass::GRID_FOG_VOLUME);

    grids[new_grid_name] = m;
}

// =============================================================================

// template <typename Q>
// void ComputeQCriterionFromGradient(double* gradients, Q qCriterion)
//{
//     // see http://public.kitware.com/pipermail/paraview/2015-May/034233.html
//     for
//     // paper citation and formula on Q-criterion.
//     qCriterion[0] =
//         -(gradients[0] * gradients[0] + gradients[4] * gradients[4] +
//         gradients[8] * gradients[8]) /
//             2. -
//         (gradients[1] * gradients[3] + gradients[2] * gradients[6] +
//         gradients[5] * gradients[7]);
// }

void crunch_qcrit(Mat3Grid::Ptr jac, GridMap& grids) {
    auto fld = convert<FloatGrid>(jac, [](Mat3 const& m) {
        return -(m[0] * m[0] + m[4] * m[4] + m[8] * m[8]) / 2. -
               (m[1] * m[3] + m[2] * m[6] + m[5] * m[7]);
    });

    auto name = "q_criterion";

    fld->setName(name);
    fld->setGridClass(openvdb::GridClass::GRID_FOG_VOLUME);
    grids[name] = fld;
}

void qcrit_3d_vel(std::string name, GridMap& grids) { }
void qcrit_1d_vel(std::string vx,
                  std::string vy,
                  std::string vz,
                  GridMap&    grids) {
    // get all
    auto vx_grid = cast_to<openvdb::FloatGrid>(find_grid_by_name(vx, grids));
    auto vy_grid = cast_to<openvdb::FloatGrid>(find_grid_by_name(vy, grids));
    auto vz_grid = cast_to<openvdb::FloatGrid>(find_grid_by_name(vz, grids));

    // compute gradients of each
    auto gx = openvdb::tools::gradient(*vx_grid);
    auto gy = openvdb::tools::gradient(*vy_grid);


    auto grad_to_v6 = [](openvdb::Vec3fGrid::Ptr&& grid) {
        return convert<Vec6Grid>(std::move(grid), [](openvdb::Vec3f v) {
            return Vec6 { v.x(), v.y(), v.z(), 0, 0, 0 };
        });
    };

    auto grad_to_v9 = [](openvdb::Vec3fGrid::Ptr&& grid) {
        return convert<Mat3Grid>(std::move(grid), [](openvdb::Vec3f v) {
            return Mat3 {
                v.x(), v.y(), v.z(), //
                0,     0,     0,     //
                0,     0,     0,     //
            };
        });
    };

    auto v6_to_v9 = [](Vec6Grid::Ptr&& grid) {
        return convert<Mat3Grid>(std::move(grid), [](Vec6 const& v) {
            return Mat3 {
                v[0], v[1], v[2], //
                v[3], v[4], v[5], //
                0,    0,    0,    //
            };
        });
    };

    auto upgraded_gx = grad_to_v6(std::move(gx));
    auto upgraded_gy = grad_to_v6(std::move(gy));

    // now combine

    upgraded_gx->tree().combine(upgraded_gy->tree(),
                                [](Vec6 const& a, Vec6 const& b, Vec6& result) {
                                    result = Vec6 {
                                        a[0], a[1], a[2], b[0], b[1], b[2],
                                    };
                                });

    // clear out memory
    upgraded_gy = nullptr;

    // and now for z
    auto gz = openvdb::tools::gradient(*vz_grid);

    auto upgraded_gz = convert<Mat3Grid>(std::move(gz), [](openvdb::Vec3f v) {
        return Mat3 {
            0,     0,     0,     //
            0,     0,     0,     //
            v.x(), v.y(), v.z(), //
        };
    });

    // now we have to change the gx tree into a mat3 as well

    auto result = convert<Mat3Grid>(std::move(upgraded_gx), [](Vec6 const& in) {
        return Mat3 {
            in[0], in[1], in[2], //
            in[3], in[4], in[5], //
            0,     0,     0,     //
        };
    });

    // now zip and return

    result->tree().combine(
        upgraded_gz->tree(), [](Mat3 const& a, Mat3 const& b, Mat3& result) {
            result =
                Mat3 { a[0], a[1], a[2], a[3], a[4], a[5], b[6], b[7], b[8] };
        });

    crunch_qcrit(result, grids);
}

void add_in_qcriterion(std::vector<std::string> names, GridMap& grids) {
    spdlog::info("Computing qcriterion...");

    // is there a vec3 velocity? or is this multiple

    switch (names.size()) {
    case 1: return qcrit_3d_vel(names.at(0), grids); break;
    case 3:
        return qcrit_1d_vel(names.at(0), names.at(1), names.at(2), grids);
        break;
    }


    //    m->setName(new_grid_name);
    //    m->setGridClass(openvdb::GridClass::GRID_FOG_VOLUME);

    //    grids[new_grid_name] = m;
}

// =============================================================================

void postprocess(toml::value const& root, GridMap& l) {
    spdlog::info("Starting postprocess");

    auto node = toml::find(root, "post");

    if (node.contains("trim")) {
        auto trim_node = toml::find(node, "trim");

        bool keep_above = false;

        if (trim_node.contains("greater")) {
            auto val   = toml::find<bool>(trim_node, "greater");
            keep_above = val;
        } else if (trim_node.contains("less")) {
            auto val   = toml::find<bool>(trim_node, "less");
            keep_above = !val;
        }

        PPVFrac trim_info {
            .name             = toml::find_or(trim_node, "variable", "v*frac"),
            .value            = toml::find_or<float>(trim_node, "value", .5),
            .keep_above_value = keep_above,
        };

        trim_all_vfrac(trim_info, l);
    }

    // do this first before velocity work
    if (node.contains("compute_qcriterion")) {
        auto qcrit_node = toml::find(node, "compute_qcriterion");

        auto vars = toml::find(qcrit_node, "variables");

        std::vector<std::string> var_names;

        for (auto const& i : vars.as_array()) {
            var_names.push_back(i.as_string());
        }

        add_in_qcriterion(std::move(var_names), l);
    }

    if (node.contains("merge_velocity")) {
        auto merge_node = toml::find(node, "merge_velocity");

        PPVelName info;

        info.x = toml::find<std::string>(merge_node, "x");
        info.y = toml::find<std::string>(merge_node, "y");
        info.z = toml::find<std::string>(merge_node, "z");

        merge_velocity(info, l);
    }

    if (node.contains("compute_magvort")) { add_in_magvort(l); }

    if (node.contains("compute_mgd")) {
        auto mgd_node = toml::find(node, "compute_mgd");

        auto mgd = toml::find_or(mgd_node, "variable", "density");

        add_in_mag_grad_density(mgd, l);
    }
}
