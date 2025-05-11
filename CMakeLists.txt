cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

set(THIRD_PARTY_DIR ${CMAKE_SOURCE_DIR}/third_party/)

list(PREPEND CMAKE_PREFIX_PATH ${THIRD_PARTY_DIR})

if (UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

# Project ======================================================================

project(
    amr2vdb
    VERSION 0.1
    LANGUAGES C CXX
)

# Depends ======================================================================

include(cmake/CPM.cmake)

# IDE helpers
CPMAddPackage(
  NAME GroupSourcesByFolder.cmake
  GITHUB_REPOSITORY TheLartians/GroupSourcesByFolder.cmake
  VERSION 1.0
)

# Create BIN ===============================================================

add_executable(amr2vdb
    src/amr_common.h
    src/argparse.h
    src/lzs3d.h
    src/main.cpp
    src/mesh_to_volume.cpp
    src/mesh_to_volume.h
    src/tricubic.h
    src/volume.cpp
    src/volume.h
    src/volume2.cpp
    src/volume2.h
    src/postprocess.h src/postprocess.cpp
    src/argparse.cpp
    src/special.h src/special.cpp
    src/points.h src/points.cpp
    src/export.h src/export.cpp
    src/pltfilereader.h src/pltfilereader.cpp
    src/volume3.h src/volume3.cpp
)

set(BUILD_SHARED_LIBS OFF)

set_target_properties(amr2vdb PROPERTIES CXX_STANDARD 20)

if (CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
    set (X86 TRUE)
    target_compile_options(amr2vdb PUBLIC -mavx)
else ()
    set (X86 FALSE)
endif ()

if  (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(amr2vdb PUBLIC -fsanitize=address)
    target_link_options(amr2vdb PUBLIC  -fsanitize=address)
endif()

# Source Files =================================================================

set(CMAKE_INCLUDE_CURRENT_DIR ON)

target_include_directories(amr2vdb PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(amr2vdb PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/include/)

# Link Deps ====================================================================

find_package(AMReX REQUIRED)
#find_package(PelePhysics REQUIRED)

list(PREPEND CMAKE_MODULE_PATH "${THIRD_PARTY_DIR}/lib/cmake/OpenVDB/")
list(PREPEND CMAKE_MODULE_PATH "${THIRD_PARTY_DIR}/lib64/cmake/OpenVDB/")
set(Boost_USE_STATIC_LIBS ON)
find_package(TBB REQUIRED)

SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")

find_library(zstd_lib zstd ${THIRD_PARTY_DIR}/lib)
find_library(z_lib NAMES z.a z ${THIRD_PARTY_DIR}/lib)
find_library(blosc_lib blosc ${THIRD_PARTY_DIR}/lib)
find_library(openvdb_lib openvdb ${THIRD_PARTY_DIR}/lib)
find_library(boost_io boost_iostreams ${THIRD_PARTY_DIR}/lib)
find_library(hdf5_cpp hdf5_cpp PATHS ${THIRD_PARTY_DIR}/lib NO_DEFAULT_PATH REQUIRED)

find_package(spdlog)

target_link_libraries(amr2vdb PRIVATE
    AMReX::amrex_3d
    TBB::tbb ${openvdb_lib} ${boost_io}
    ${blosc_lib} spdlog::spdlog ${zstd_lib} ${z_lib} ${hdf5_cpp}
)

MESSAGE( ${hdf5_cpp})

# Install ======================================================================

include(GNUInstallDirs)

set(THIRD_PARTY_LIB_DIR "${THIRD_PARTY_DIR}/lib")
set(THIRD_PARTY_PLUGIN_DIR "${THIRD_PARTY_DIR}/plugin")

file(GLOB_RECURSE THIRD_PARTY_LIBS
    "${THIRD_PARTY_LIB_DIR}/*.so"
    "${THIRD_PARTY_LIB_DIR}/*.dylib"
)

file(GLOB_RECURSE THIRD_PARTY_PLUGINS
    "${THIRD_PARTY_PLUGIN_DIR}/*.so"
    "${THIRD_PARTY_PLUGIN_DIR}/*.dylib"
)

install(FILES ${THIRD_PARTY_LIBS} ${THIRD_PARTY_PLUGINS}
    DESTINATION ${CMAKE_INSTALL_BINDIR}
)

set_target_properties(amr2vdb PROPERTIES
    INSTALL_RPATH "\$ORIGIN"
    BUILD_WITH_INSTALL_RPATH TRUE
)

# Install the binary
install(TARGETS amr2vdb
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
