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

# we need fortran for AMReX data services. Don't ask me why.

# hack for clangd
#set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
#if(CMAKE_EXPORT_COMPILE_COMMANDS)
#  set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES
#      ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
#endif()
# end hack for clangd

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
)

set(BUILD_SHARED_LIBS OFF)

set_target_properties(amr2vdb PROPERTIES CXX_STANDARD 20)

if (CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
    set (X86 TRUE)
    target_compile_options(amr2vdb PUBLIC -mavx)
else ()
    set (X86 FALSE)
endif ()

target_compile_options(amr2vdb PUBLIC -fsanitize=address)
target_link_options(amr2vdb PUBLIC  -fsanitize=address)

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
find_package(OpenVDB REQUIRED)

find_library(blosc_lib blosc ${THIRD_PARTY_DIR}/lib)

find_package(spdlog)

target_link_libraries(amr2vdb PRIVATE
    AMReX::amrex_3d
    PelePhysics::PelePhysics-drm19-Fuego-Simple
    TBB::tbb OpenVDB::openvdb
    ${blosc_lib} spdlog::spdlog
)



# IDE Helpers ==================================================================

# Install ======================================================================
