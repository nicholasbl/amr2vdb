#pragma once

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>

// https://github.com/dfelinto/blender/blob/87a0770bb969ce37d9a41a04c1658ea09c63933a/source/blender/draw/engines/workbench/shaders/workbench_volume_frag.glsl#L23

struct Tricubic {
    const char* name() { return "Tricubic"; }
    static int  radius() { return 2; }
    static bool mipmap() { return true; }
    static bool consistent() { return false; }

    openvdb::tools::BoxSampler box_sampler;

    template <class TreeT>
    bool sample(TreeT const&               tree,
                openvdb::Vec3R const&      coord,
                typename TreeT::ValueType& val) const;
};


template <class TreeT>
bool Tricubic::sample(TreeT const&               tree,
                      openvdb::Vec3R const&      coord,
                      typename TreeT::ValueType& val) const {

    // auto tc = (coord - .5);
    // tc      = { std::floor(tc.x()), std::floor(tc.y()), std::floor(tc.z()) };
    // tc += openvdb::Vec3R(.5);

    openvdb::Vec3R tc = { std::floor(coord.x()),
                          std::floor(coord.y()),
                          std::floor(coord.z()) };

    auto f  = tc - coord;
    auto f2 = f * f;
    auto f3 = f2 * f;

    auto w3 = f3 / 6.0;
    auto w0 = -w3 + f2 * 0.5 - f * 0.5 + 1.0 / 6.0;
    auto w1 = f3 * 0.5 - f2 + 2.0 / 3.0;
    auto w2 = openvdb::Vec3R(1.0) - w0 - w1 - w3;

    auto s0 = w0 + w1;
    auto s1 = w2 + w3;

    auto f0 = w1 / (w0 + w1);
    auto f1 = w3 / (w2 + w3);


    openvdb::Vec4R final_co = {
        tc.x() - 1.0 + f0.x(),
        tc.y() - 1.0 + f0.y(),
        tc.x() + 1.0 + f1.x(),
        tc.y() + 1.0 + f1.y(),
    };

    openvdb::Vec2R final_z = {
        tc.z() - 1.0 + f0.z(),
        tc.z() + 1.0 + f1.z(),
    };

    // color = texture(ima, vec3(final_co.xy, final_z.x)) * s0.x * s0.y * s0.z;
    auto ret =
        box_sampler.sample(tree, { final_co.x(), final_co.y(), final_z.x() }) *
        s0.x() * s0.y() * s0.z();

    // color += texture(ima, vec3(final_co.zy, final_z.x)) * s1.x * s0.y * s0.z;
    ret +=
        box_sampler.sample(tree, { final_co.z(), final_co.y(), final_z.x() }) *
        s1.x() * s0.y() * s0.z();

    // color += texture(ima, vec3(final_co.xw, final_z.x)) * s0.x * s1.y * s0.z;
    ret +=
        box_sampler.sample(tree, { final_co.x(), final_co.w(), final_z.x() }) *
        s0.x() * s1.y() * s0.z();

    // color += texture(ima, vec3(final_co.zw, final_z.x)) * s1.x * s1.y * s0.z;
    ret +=
        box_sampler.sample(tree, { final_co.z(), final_co.w(), final_z.x() }) *
        s1.x() * s1.y() * s0.z();


    // color += texture(ima, vec3(final_co.xy, final_z.y)) * s0.x * s0.y * s1.z;
    ret +=
        box_sampler.sample(tree, { final_co.x(), final_co.y(), final_z.y() }) *
        s0.x() * s0.y() * s1.z();

    // color += texture(ima, vec3(final_co.zy, final_z.y)) * s1.x * s0.y * s1.z;
    ret +=
        box_sampler.sample(tree, { final_co.z(), final_co.y(), final_z.y() }) *
        s1.x() * s0.y() * s1.z();

    // color += texture(ima, vec3(final_co.xw, final_z.y)) * s0.x * s1.y * s1.z;
    ret +=
        box_sampler.sample(tree, { final_co.x(), final_co.w(), final_z.y() }) *
        s0.x() * s1.y() * s1.z();

    // color += texture(ima, vec3(final_co.zw, final_z.y)) * s1.x * s1.y * s1.z;
    ret +=
        box_sampler.sample(tree, { final_co.z(), final_co.w(), final_z.y() }) *
        s1.x() * s1.y() * s1.z();

    val = ret;

    return true;
}
