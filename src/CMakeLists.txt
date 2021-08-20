
set(stream
    "${PROJECT_SOURCE_DIR}/src/Grid/EllipticalRadialGridGenerator.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/DefineOptions.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/ElectricField.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/EquationSystem.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Grid.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Heat.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Ions.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Process.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/n_re.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/ConfinementTime.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/ElectronTransportDiffusion.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/ElectronHeatTransport.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/ElectronHeatTransportDiffusion.cpp"
#    "${PROJECT_SOURCE_DIR}/src/Equations/IonTransportDiffusion1D.cpp"
#    "${PROJECT_SOURCE_DIR}/src/Equations/IonHeatTransportDiffusion1D.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/IonTransport.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/IonHeatTransport.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Transport.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/SputteredRecycledCoefficient.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/NeutralInflux.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/RunawayElectronConfinementTime.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/MeanFreePathTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Transport.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Constants.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/PlasmaVolume.cpp" 
    "${PROJECT_SOURCE_DIR}/src/Settings/MeanFreePath.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/IonRateEquation.cpp"
    "${PROJECT_SOURCE_DIR}/src/EquationSystem.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/NeutralInflux.cpp" 
    "${PROJECT_SOURCE_DIR}/src/Equations/NeutralTransport.cpp" 
    "${PROJECT_SOURCE_DIR}/src/Equations/RadiatedPowerTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/ChargeExchangeTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/OtherQuantityHandler.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/OtherQuantityHandler.cpp"
)

set_property(
    SOURCE "${PROJECT_SOURCE_DIR}/src/OtherQuantityHandler.cpp"
    PROPERTY COMPILE_FLAGS -Wno-unused-parameter
)

add_library(stream STATIC ${stream})
target_include_directories(stream PUBLIC "${PROJECT_BINARY_DIR}/include" "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(stream PUBLIC fvm)
target_link_libraries(stream PUBLIC dream)

set_target_properties(stream PROPERTIES
    # Require C++17
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
)

if (POS_INDP_CODE)
    # Make code position-independent (required for linking to Python interface)
    set_property(TARGET stream PROPERTY POSITION_INDEPENDENT_CODE YES)
endif ()

######################
# MODULES
######################
# Math library
find_library(MATH_LIBRARY m)
if (MATH_LIBRARY)
    target_link_libraries(stream PUBLIC ${MATH_LIBRARY})
endif()

# Interprocedural optimizations
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    include(CheckIPOSupported)
    # Clang does not support IPO (although Clang++ does), so we need
    # to explicitly check both C and C++ (since Clang is used for linking)
    check_ipo_supported(RESULT result LANGUAGES C CXX)
    if (result)
        set_target_properties(stream PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endif ()

# GNU Scientific Library
find_package(GSL REQUIRED)
if (GSL_FOUND)
    if (GSL_VERSION VERSION_LESS 2.0)
        message(FATAL_ERROR "stream needs GSL >= 2.0")
    endif()

    target_include_directories(stream PUBLIC "${GSL_INCLUDE_DIRS}")
    target_link_libraries(stream PUBLIC "${GSL_LIBRARIES}")
    target_link_libraries(stream PUBLIC "${GSL_CBLAS_LIBRARY}")
endif()

find_package(HDF5 COMPONENTS CXX)
if (HDF5_FOUND)
    target_include_directories(stream PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(stream PUBLIC ${HDF5_LIBRARIES})
else (HDF5_FOUND)
    message(FATAL_ERROR "No HDF5 C++ library was found on this system")
endif (HDF5_FOUND)

# MPI (just as a dependency for PETSc)
#find_package(MPI COMPONENTS CXX)
#if (MPI_CXX_FOUND)
#    target_include_directories(dream PUBLIC ${MPI_CXX_INCLUDE_DIRS})
#	set_target_properties(dream PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_OPTIONS} ${MPI_CXX_COMPILE_DEFINITIONS}")
#	set_target_properties(dream PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
#	target_link_libraries(dream PUBLIC "${MPI_CXX_LIBRARIES}")
#endif (MPI_CXX_FOUND)

# PETSc
#find_package(PETSc COMPONENTS CXX REQUIRED)
#if (PETSC_FOUND)
#    target_include_directories(dream PUBLIC "${PETSC_INCLUDES}")
#    target_link_libraries(dream PUBLIC "${PETSC_LIBRARIES}")
#    add_definitions(${PETSC_DEFINITIONS})
#endif()

