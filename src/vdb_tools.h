#pragma once

#include <iostream>
#include <optional>

#include <openvdb/openvdb.h>

template <class T>
struct Pair {
    T first, second;
};

template <class T>
Pair<T> make_pair(T a, T b) {
    return Pair<T> { a, b };
}

template <class Reader, class IterA>
auto vdb_chunk(Reader const& a,
               Pair<size_t>  xs,
               Pair<size_t>  ys,
               Pair<IterA>   zs) {
    auto sub_grid = openvdb::FloatGrid::create();
    auto accessor = sub_grid->getAccessor();

    openvdb::Coord ijk;

    int& x = ijk[0];
    int& y = ijk[1];
    int& z = ijk[2];

    using RetType = std::invoke_result_t<Reader, size_t, size_t, size_t>;

    for (z = zs.first; z != zs.second; ++z) {
        for (y = ys.first; y < ys.second; ++y) {
            for (x = xs.first; x < xs.second; ++x) {
                float value = 0.0;
                bool  ok    = a(x, y, z, value);

                if (ok) { accessor.setValue(ijk, value); }
            }
        }

        std::cout << "P: " << z << "/" << zs.second - 1 << std::endl;
    }
    return sub_grid;
}

template <class Reader>
[[nodiscard]] auto build_open_vdb(std::array<size_t, 3> dims, Reader const& a) {
    std::cout << "Starting VDB build..." << std::endl;

    return vdb_chunk(
        a, { 0, dims[0] }, { 0, dims[1] }, make_pair(size_t { 0 }, dims[2]));
}
