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

CPMAddPackage(
  NAME lua
  GIT_REPOSITORY https://github.com/lua/lua.git
  VERSION 5.4.6
  DOWNLOAD_ONLY YES
)

if(lua_ADDED)
  file(GLOB lua_sources ${lua_SOURCE_DIR}/*.c)
  list(REMOVE_ITEM lua_sources "${lua_SOURCE_DIR}/lua.c" "${lua_SOURCE_DIR}/luac.c")
  add_library(lua STATIC ${lua_sources})
  target_include_directories(lua SYSTEM PUBLIC $<BUILD_INTERFACE:${lua_SOURCE_DIR}>)
endif()

CPMAddPackage("gh:ThePhD/sol2@3.3.0")

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
    src/amr_basic.cpp src/amr_basic.h
    src/use_sol.h
)

set(BUILD_SHARED_LIBS OFF)

set_target_properties(amr2vdb PROPERTIES CXX_STANDARD 20)

if (CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
    set (X86 TRUE)
    target_compile_options(amr2vdb PUBLIC -mavx)
else ()
    set (X86 FALSE)
endif ()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(amr2vdb PUBLIC -fsanitize=address)
    target_link_options(amr2vdb PUBLIC  -fsanitize=address)
endif()

# Source Files =================================================================

set(CMAKE_INCLUDE_CURRENT_DIR ON)

target_include_directories(amr2vdb PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(amr2vdb PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/include/)

# Link Deps ====================================================================

find_package(AMReX REQUIRED)
find_package(PelePhysics REQUIRED)

list(PREPEND CMAKE_MODULE_PATH "${THIRD_PARTY_DIR}/lib/cmake/OpenVDB/")
list(PREPEND CMAKE_MODULE_PATH "${THIRD_PARTY_DIR}/lib64/cmake/OpenVDB/")
set(Boost_USE_STATIC_LIBS ON)
find_package(TBB REQUIRED)

find_library(zstd_lib zstd ${THIRD_PARTY_DIR}/lib)
find_library(z_lib z ${THIRD_PARTY_DIR}/lib)
find_library(blosc_lib blosc ${THIRD_PARTY_DIR}/lib)
find_library(openvdb_lib openvdb ${THIRD_PARTY_DIR}/lib)
find_library(boost_io boost_iostreams ${THIRD_PARTY_DIR}/lib)

find_package(spdlog)

target_link_libraries(amr2vdb PRIVATE
    AMReX::amrex_3d
    PelePhysics::PelePhysics-drm19-Fuego-Simple
    TBB::tbb ${openvdb_lib} ${boost_io}
    ${blosc_lib} spdlog::spdlog ${zstd_lib} ${z_lib}
    sol2 lua
)

# IDE Helpers ==================================================================

# Install ======================================================================
