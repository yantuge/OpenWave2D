# OpenWave2D - 2D Elastic Wave Simulation (BUG...PROJECT IS DEVELOPING)

A high-performance 2D seismic wave simulation framework using C++/Fortran hybrid architecture with PML absorbing boundaries.

## Overview

SeismicHybrid is a modular seismic wave simulation system that combines:
- **C++ Control Layer**: Configuration management, I/O, memory management, and simulation orchestration
- **Fortran Computational Engine**: High-performance numerical kernels for wave propagation

## Features

- **Multiple PML Types**: CPML (Convolutional PML) and ADE-PML (Auxiliary Differential Equation PML)
- **Flexible Model Support**: Homogeneous, anisotropic, and binary file-based velocity models
- **Source Time Functions**: Ricker wavelet, Gaussian, and first derivative sources
- **Time Integration Schemes**: Explicit finite-difference and 4th-order Runge-Kutta
- **Output Formats**: Binary snapshots and ASCII/SU seismograms
- **JSON Configuration**: Flexible parameter specification via JSON files

## Architecture

```
SeismicHybrid/
├── build/                    # Build directory
├── config/                   # Configuration files
├── fortran/                  # Fortran computational kernels
├── include/seismic/          # C++ headers
│   ├── core/                 # Core data structures
│   ├── io/                   # I/O functionality  
│   ├── initializers/         # Initialization modules
│   ├── sources/              # Source time functions
│   └── main/                 # Main dispatcher
├── src/                      # C++ source files
├── CMakeLists.txt           # Build configuration
└── README.md                # This file
```

## Build Instructions

### Prerequisites

- CMake 3.16 or higher
- Modern C++ compiler (C++17 support)
- Fortran compiler (gfortran, ifort, etc.)
- nlohmann/json library

### Building

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

### Dependencies

For Ubuntu/Debian:
```bash
sudo apt-get install cmake gfortran nlohmann-json3-dev
```

For Windows with MSYS2:
```bash
pacman -S mingw-w64-x86_64-cmake mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-nlohmann-json
```

## Usage

### Basic Usage

```bash
./OpenWave2D ../config/elastic_homogeneous.json
```

### Configuration

Edit the JSON configuration file to customize simulation parameters:

```json
{
  "grid": {
    "nx": 401,
    "ny": 401,
    "deltax": 0.000625,
    "deltay": 0.000625
  },
  "model": {
    "type": "homogeneous",
    "vp": 2000.0,
    "vs": 1150.0,
    "density": 2000.0
  },
  "time": {
    "nstep": 3000,
    "dt": 5.0e-8,
    "scheme": "explicit"
  },
  "source": {
    "x": 0.125,
    "y": 0.125,
    "f0": 200000.0,
    "type": "ricker"
  },
  "pml": {
    "type": "cpml",
    "npoints_pml": 10
  }
}
```

## Model Types

### Homogeneous Model
```json
{
  "model": {
    "type": "homogeneous",
    "vp": 2000.0,
    "vs": 1150.0, 
    "density": 2000.0
  }
}
```

### Anisotropic Model
```json
{
  "model": {
    "type": "anisotropic",
    "c11": 4.0e10,
    "c12": 3.8e10,
    "c22": 20.0e10,
    "c33": 2.0e10,
    "density": 4000.0
  }
}
```

### Binary File Model
```json
{
  "model": {
    "type": "binary_file",
    "vp_file": "vp_model.bin",
    "vs_file": "vs_model.bin",
    "density_file": "density_model.bin"
  }
}
```

## PML Boundary Conditions

### CPML (Convolutional PML)
- Based on Roden and Gedney (2000)
- Second-order accurate in space and time
- Suitable for most applications

### ADE-PML (Auxiliary Differential Equation PML)
- Based on Roden and Gedney (2010)
- 8th-order spatial accuracy with Holberg coefficients
- 4th-order Runge-Kutta time integration
- Higher computational cost but better accuracy

## Output Files

### Snapshots
- Binary format: `snapshot_vx_XXXX.bin`, `snapshot_vy_XXXX.bin`
- Grid dimensions: nx × ny
- Data type: double precision

### Seismograms
- ASCII format: Time series with receiver columns
- SU format: Standard seismic Unix format

## Performance Notes

- Use Release build for production runs: `cmake -DCMAKE_BUILD_TYPE=Release ..`
- Compiler optimization flags are automatically applied
- Memory layout is optimized for Fortran column-major access
- PML profiles are pre-computed for efficiency

## References

1. Madariaga, R. (1976). Dynamics of an expanding circular fault. BSSA, 66(3), 639-666.
2. Virieux, J. (1986). P-SV wave propagation in heterogeneous media. Geophysics, 51(4), 889-901.
3. Roden, J. A., & Gedney, S. D. (2000). Convolutional PML (CPML). IEEE TAP, 48(3), 418-425.
4. Holberg, O. (1987). Computational aspects of the choice of operator and sampling interval for numerical differentiation. Geophysical Prospecting, 35(6), 629-655.

## License

This project is provided as-is for educational and research purposes.
