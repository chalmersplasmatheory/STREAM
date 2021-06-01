
set(streami_files
    "${PROJECT_SOURCE_DIR}/iface/Main.cpp"
)

add_executable(streami ${streami_files})
target_link_libraries(streami PUBLIC stream)
set_target_properties(streami PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    include(CheckIPOSupported)
    # Clang does not support IPO (although Clang++ does), so we need
    # to explicitly check both C and C++ (since Clang is used for linking)
    check_ipo_supported(RESULT result LANGUAGES C CXX)
    if (result)
        set_target_properties(streami PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endif ()

