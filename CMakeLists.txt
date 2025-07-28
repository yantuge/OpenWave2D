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

# C++ source files
set(CXX_SOURCES
    src/Config.cpp
    src/Grid.cpp
    src/PMLInitializer.cpp
    src/ModelLoader.cpp
    src/Recorder.cpp
    src/Writer.cpp
    src/Dispatcher.cpp
    src/main.cpp
)

# Main executable
add_executable(openwave2d ${CXX_SOURCES})

# Link libraries
target_link_libraries(openwave2d PRIVATE seismic_kernels)

# Include Fortran modules
target_include_directories(SeismicHybrid PRIVATE ${CMAKE_BINARY_DIR}/fortran_modules)

# Try to find nlohmann/json
find_package(nlohmann_json QUIET)
if(nlohmann_json_FOUND)
    target_link_libraries(SeismicHybrid PRIVATE nlohmann_json::nlohmann_json)
else()
    # If not found, assume it's header-only and available in system include paths
    message(STATUS "nlohmann/json not found via find_package, assuming header-only library")
endif()

# Install target
install(TARGETS SeismicHybrid DESTINATION bin)
install(TARGETS seismic_kernels DESTINATION lib)

# Create output directories
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/output)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/seismograms)

# Copy example configuration
configure_file(config/elastic_homogeneous.json ${CMAKE_BINARY_DIR}/elastic_homogeneous.json COPYONLY)
