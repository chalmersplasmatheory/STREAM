cmake_minimum_required(VERSION 3.12..3.14)

if (${CMAKE_VERSION} VERSION_LESS 3.14)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
    cmake_policy(VERSION 3.14)
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Declare project
project(stream VERSION 1.0
               DESCRIPTION "Startup Runaway Electron Analysis Model"
               LANGUAGES CXX)

option(COLOR_TERMINAL "Allow colourful output" ON)
option(GIT_SUBMODULE "Check submodules during build" ON)

# Disable DREAM settings
set(DREAM_BUILD_TESTS OFF)
set(DREAM_BUILD_PYFACE OFF)

# Add CMake module path
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" "${PROJECT_SOURCE_DIR}/extern/DREAM/cmake" ${CMAKE_MODULE_PATH})

# Get git commit hash
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# Generate 'config.h'
configure_file(
    "${PROJECT_SOURCE_DIR}/include/stream.h.in"
    "${PROJECT_BINARY_DIR}/include/STREAM/stream.h"
)

# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if (EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory.")
endif()

# Find git
find_package(Git QUIET)
if (GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    if (GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            RESULT_VARIABLE GIT_SUBMOD_RESULT)

        if (NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

# Enable compilation warnings on all targets
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -Wno-ignored-qualifiers)
    if (CMAKE_BUILD_TYPE STREQUAL "Release")
        add_compile_options(-O3)
    endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-Wall -Wextra -Wpedantic -Wformat-security -Wno-misleading-indentation -Wno-ignored-qualifiers)
    if (CMAKE_BUILD_TYPE STREQUAL "Release")
        add_compile_options(-O3)
    endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    add_compile_options(-Wall -Wextra -Wpedantic -Wno-misleading-indentation -Wno-ignored-qualifiers)
    if (CMAKE_BUILD_TYPE STREQUAL "Release")
        add_compile_options(-O3)
    endif()
else ()
    message(WARNING "You are using an un-tested compiler.")
endif ()


# Build libDREAM
add_subdirectory(extern/DREAM)

# Sub-directories
add_subdirectory(src)
add_subdirectory(iface)

