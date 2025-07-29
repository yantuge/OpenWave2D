cmake_minimum_required(VERSION 3.16)
project(openwave2d CXX Fortran)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set Fortran compiler flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")

# Find required packages
find_package(PkgConfig REQUIRED)

# Add include directories
include_directories(include)

# Fortran library
add_library(seismic_kernels STATIC
    fortran/seismic_kernels.f90
)

# Set Fortran module directory
set_target_properties(seismic_kernels PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran_modules
)

# C++ source files organized by category
set(CXX_SOURCES
    # Core components
    src/seismic/core/Config.cpp
    src/seismic/core/Grid.cpp
    
    # Initializers
    src/seismic/initializers/PMLInitializer.cpp
    
    # IO components
    src/seismic/io/ModelLoader.cpp
    src/seismic/io/Recorder.cpp
    src/seismic/io/Writer.cpp
    
    # Main application
    src/seismic/main/Dispatcher.cpp
    src/seismic/main/main.cpp
    
    # Solvers
    src/seismic/solvers/RK4Integrator.cpp
    src/seismic/solvers/CPMLIntegrator.cpp
    src/seismic/solvers/TimeIntegratorFactory.cpp
    
    # Note: Sources are header-only implementations
)

# Main executable
add_executable(OpenWave2D ${CXX_SOURCES})

# Link libraries
target_link_libraries(OpenWave2D PRIVATE seismic_kernels)

# Include Fortran modules
target_include_directories(OpenWave2D PRIVATE ${CMAKE_BINARY_DIR}/fortran_modules)

# Try to find nlohmann/json
find_package(nlohmann_json QUIET)
if(nlohmann_json_FOUND)
    target_link_libraries(OpenWave2D PRIVATE nlohmann_json::nlohmann_json)
else()
    # If not found, assume it's header-only and available in system include paths
    message(STATUS "nlohmann/json not found via find_package, assuming header-only library")
endif()

# Install target
install(TARGETS OpenWave2D DESTINATION bin)
install(TARGETS seismic_kernels DESTINATION lib)

# Create output directories
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/output)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/seismograms)

# Copy example configuration
configure_file(config/elastic_homogeneous.json ${CMAKE_BINARY_DIR}/elastic_homogeneous.json COPYONLY)
