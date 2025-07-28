#include "seismic/io/Writer.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <algorithm>

namespace seismic {

void save_snapshot(const Grid& grid, const OutputConfig& output_config, 
                   Index it, const std::string& component) {
    // Create output directory if it doesn't exist
    std::filesystem::create_directories(output_config.output_dir);
    
    const std::string filename = output_config.output_dir + "/snapshot_" + 
                                component + "_" + std::to_string(it) + ".bin";
    
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error: Cannot create snapshot file " << filename << std::endl;
        return;
    }
    
    const Size grid_size = static_cast<Size>(grid.nx() * grid.ny());
    
    if (component == "vx") {
        file.write(reinterpret_cast<const char*>(grid.get_vx_ptr()), 
                  grid_size * sizeof(Real));
    } else if (component == "vy") {
        file.write(reinterpret_cast<const char*>(grid.get_vy_ptr()), 
                  grid_size * sizeof(Real));
    }
    
    file.close();
}

void save_snapshot_su(const Grid& grid, const OutputConfig& output_config,
                         Index it, const std::string& component) {
    // Create output directory if it doesn't exist
    std::filesystem::create_directories(output_config.output_dir);
    
    const std::string filename = output_config.output_dir + "/snapshot_" + 
                                component + "_" + std::to_string(it) + ".su";
    
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error: Cannot create SU snapshot file " << filename << std::endl;
        return;
    }
    
    // Simple SU header (240 bytes)
    char su_header[240] = {0};
    int* int_ptr = reinterpret_cast<int*>(su_header);
    
    // Basic SU header fields
    int_ptr[0] = 1;  // tracl (trace number within line)
    int_ptr[1] = 1;  // tracr (trace number within reel) 
    int_ptr[28] = static_cast<int>(grid.nx());  // ns (number of samples)
    int_ptr[29] = 1000;  // dt (sample interval in microseconds)
    
    file.write(su_header, 240);
    
    const Size grid_size = static_cast<Size>(grid.nx() * grid.ny());
    
    if (component == "vx") {
        file.write(reinterpret_cast<const char*>(grid.get_vx_ptr()), 
                  grid_size * sizeof(Real));
    } else if (component == "vy") {
        file.write(reinterpret_cast<const char*>(grid.get_vy_ptr()), 
                  grid_size * sizeof(Real));
    }
    
    file.close();
}

// SU format implementation would go in save_snapshot_su above

void write_seismograms_ascii(const std::vector<std::vector<Real>>& seismograms,
                            const std::string& filename, Real dt) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error: Cannot create seismogram file " << filename << std::endl;
        return;
    }
    
    if (seismograms.empty()) return;
    
    const Size nstep = seismograms[0].size();
    const Size nrec = seismograms.size();
    
    // Write header
    file << "# Time (s)";
    for (Size r = 0; r < nrec; ++r) {
        file << "  Receiver_" << r;
    }
    file << std::endl;
    
    // Write data
    file << std::scientific << std::setprecision(6);
    for (Size it = 0; it < nstep; ++it) {
        const Real time = it * dt;
        file << time;
        
        for (Size r = 0; r < nrec; ++r) {
            file << "  " << seismograms[r][it];
        }
        file << std::endl;
    }
    
    file.close();
}

void write_seismograms_su(const std::vector<std::vector<Real>>& seismograms,
                         const std::string& filename, Real dt) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error: Cannot create SU seismogram file " << filename << std::endl;
        return;
    }
    
    if (seismograms.empty()) return;
    
    const Size nstep = seismograms[0].size();
    const Size nrec = seismograms.size();
    
    // Simple SU header (240 bytes, all zeros except required fields)
    struct SUHeader {
        char data[240];
    };
    
    for (Size r = 0; r < nrec; ++r) {
        SUHeader header;
        std::fill_n(header.data, 240, 0);
        
        // Set number of samples (at byte offset 58)
        *reinterpret_cast<short*>(&header.data[58]) = static_cast<short>(nstep);
        
        // Set sample interval in microseconds (at byte offset 59)
        *reinterpret_cast<short*>(&header.data[59]) = static_cast<short>(dt * 1.0e6);
        
        // Write header
        file.write(header.data, 240);
        
        // Write trace data as float
        for (Size it = 0; it < nstep; ++it) {
            const float value = static_cast<float>(seismograms[r][it]);
            file.write(reinterpret_cast<const char*>(&value), sizeof(float));
        }
    }
    
    file.close();
}

} // namespace seismic
