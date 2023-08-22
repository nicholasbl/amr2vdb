#pragma once

#include <openvdb/openvdb.h>


struct LZS3D {
    // we are using Lanczos 3 3D
    // paper suggests this has less artifacts than a = 2
    static constexpr int lz_a = 3;

    const char* name() { return "Lanczos 3D"; }
    static int  radius() { return lz_a; }
    static bool mipmap() { return true; }
    static bool consistent() { return false; }

    template <class TreeT>
    bool sample(TreeT const&               tree,
                openvdb::Vec3R const&      coord,
                typename TreeT::ValueType& val) const;
};

template <int LA, class T>
T L_func(T x) {
    if (x == 0) return 1;
    if (-LA < x and x < LA) {
        static constexpr double PISQ = M_PI * M_PI;
        return LA * std::sin(M_PI * x) * std::sin(M_PI * x / LA) /
               (PISQ * x * x);
    }
    return T(0);
}

template <class TreeT>
bool LZS3D::sample(TreeT const&               tree,
                   openvdb::Vec3R const&      coord,
                   typename TreeT::ValueType& val) const {

    int xi = std::floor(coord.x()) - lz_a + 1;
    int yi = std::floor(coord.y()) - lz_a + 1;
    int zi = std::floor(coord.z()) - lz_a + 1;

    int xf = std::floor(coord.x()) + lz_a;
    int yf = std::floor(coord.y()) + lz_a;
    int zf = std::floor(coord.z()) + lz_a;

    // static constexpr int size = lz_a * 2 - 1;
    static constexpr int size = lz_a * 2;

    int temp_x[size][size] = {};
    int temp_y[size]       = {};

    int m = 0;
    int n = 0;

    for (int k = zi; k <= zf; k++) {
        n = 0;
        for (int j = yi; j <= yf; j++) {
            for (int i = xi; i <= xf; i++) {
                // assert(m < size);
                // assert(n < size);
                temp_x[m][n] += tree.getValue(openvdb::Coord { i, j, k }) *
                                L_func<lz_a>(coord.x() - i);
            }
            n++;
        }
        m++;
    }

    m = 0;

    for (int k = zi; k <= zf; k++) {
        n = 0;
        for (int j = yi; j <= yf; j++) {
            // assert(m < size);
            // assert(n < size);
            temp_y[m] = temp_x[m][n] * L_func<lz_a>(coord.x() - j);
            n++;
        }
        m++;
    }

    m = 0;

    double value = 0;

    for (int k = zi; k <= zf; k++) {
        // assert(m < size);
        value += temp_y[m] * L_func<lz_a>(coord.z() - k);
        m++;
    }

    val = value;

    return true;
}
