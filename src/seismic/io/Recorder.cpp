#include "seismic/io/Recorder.hpp"
#include "seismic/io/Writer.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

namespace seismic {

SeismogramRecorder::SeismogramRecorder(const ReceiverConfig& receiver_config, Index nstep) 
    : m_nstep(nstep), m_recording_mode(RecordingMode::FOUR_POINT_AVG) {  // Default to 4-point average
    // 计算接收器数量（x和y位置数组中较小的那个）
    Size num_receivers = std::min(receiver_config.x_positions.size(), receiver_config.y_positions.size());
    
    m_receivers.reserve(num_receivers);
    m_seismograms_vx.resize(num_receivers);
    m_seismograms_vy.resize(num_receivers);
    
    for (Size i = 0; i < num_receivers; ++i) {
        m_receivers.push_back({
            0, 0,  // i, j 将在initialize方法中设置
            receiver_config.x_positions[i],
            receiver_config.y_positions[i]
        });
        
        m_seismograms_vx[i].resize(nstep, 0.0);
        m_seismograms_vy[i].resize(nstep, 0.0);
    }
}

void SeismogramRecorder::initialize(const GridConfig& grid_config) {
    // 将接收器位置转换为网格索引
    for (Size r = 0; r < m_receivers.size(); ++r) {
        // 找到最接近的网格点
        const Index best_i = static_cast<Index>(m_receivers[r].x / grid_config.deltax + 0.5);
        const Index best_j = static_cast<Index>(m_receivers[r].y / grid_config.deltay + 0.5);
        
        m_receivers[r].i = best_i;
        m_receivers[r].j = best_j;
        
        // 更新实际位置为网格点位置
        m_receivers[r].x = best_i * grid_config.deltax;
        m_receivers[r].y = best_j * grid_config.deltay;
    }
}

void SeismogramRecorder::set_recording_mode(const std::string& solver_type) {
    if (solver_type == "cpml") {
        m_recording_mode = RecordingMode::SINGLE_POINT;
    } else {
        m_recording_mode = RecordingMode::FOUR_POINT_AVG;  // Default for RK4/ADE-PML
    }
}

void SeismogramRecorder::record_time_step(const Grid& grid, Index it) {
    for (Size r = 0; r < m_receivers.size(); ++r) {
        const Index i = m_receivers[r].i;
        const Index j = m_receivers[r].j;
        
        if (it < m_nstep) {
            if (m_recording_mode == RecordingMode::FOUR_POINT_AVG) {
                // ADE-PML模式：使用原始Fortran代码的4点平均方法 (对vx应用)
                // vx使用4点平均: (vx(i,j) + vx(i+1,j) + vx(i,j+1) + vx(i+1,j+1))/4.0
                if (i + 1 < grid.nx() && j + 1 < grid.ny()) {
                    m_seismograms_vx[r][it] = (grid.vx(i, j) + grid.vx(i+1, j) + 
                                              grid.vx(i, j+1) + grid.vx(i+1, j+1)) / 4.0;
                } else {
                    m_seismograms_vx[r][it] = grid.vx(i, j);
                }
                
                // vy直接使用单点值
                m_seismograms_vy[r][it] = grid.vy(i, j);
            } else {
                // CPML模式：使用单点采样
                m_seismograms_vx[r][it] = grid.vx(i, j);
                m_seismograms_vy[r][it] = grid.vy(i, j);
            }
        }
    }
}

void SeismogramRecorder::save_seismograms(const OutputConfig& output_config, Real dt) const {
    if (!output_config.save_seismograms) return;
    
    if (output_config.seismogram_format == "fortran" || output_config.seismogram_format == "ascii_fortran") {
        save_seismograms_fortran_style(output_config, dt);
    } else {
        // 保持原有的输出格式
        const std::string vx_filename = output_config.output_dir + "/seismograms_vx." + output_config.seismogram_format;
        const std::string vy_filename = output_config.output_dir + "/seismograms_vy." + output_config.seismogram_format;
        
        if (output_config.seismogram_format == "ascii") {
            write_seismograms_ascii(m_seismograms_vx, vx_filename, dt);
            write_seismograms_ascii(m_seismograms_vy, vy_filename, dt);
        } else if (output_config.seismogram_format == "su") {
            write_seismograms_su(m_seismograms_vx, vx_filename, dt);
            write_seismograms_su(m_seismograms_vy, vy_filename, dt);
        }
    }
}

void SeismogramRecorder::save_seismograms_fortran_style(const OutputConfig& output_config, Real dt) const {
    // 按照原始Fortran代码的格式：每个接收器一个文件，每行为 "时间 值"
    
    // X分量 - 按照Fortran格式: Vx_file_001.dat, Vx_file_002.dat, ...
    for (Size r = 0; r < m_receivers.size(); ++r) {
        std::ostringstream filename;
        filename << output_config.output_dir << "/Vx_file_" 
                 << std::setfill('0') << std::setw(3) << (r + 1) << ".dat";
        
        std::ofstream file(filename.str());
        if (!file) {
            std::cerr << "Error: Cannot create seismogram file " << filename.str() << std::endl;
            continue;
        }
        
        for (Index it = 0; it < m_nstep; ++it) {
            const Real time = it * dt;
            // 使用与Fortran代码相同的格式：科学计数法 + 单精度
            file << std::scientific << std::setprecision(6) 
                 << time << " " << static_cast<float>(m_seismograms_vx[r][it]) << std::endl;
        }
        file.close();
    }
    
    // Y分量 - 按照Fortran格式: Vy_file_001.dat, Vy_file_002.dat, ...
    for (Size r = 0; r < m_receivers.size(); ++r) {
        std::ostringstream filename;
        filename << output_config.output_dir << "/Vy_file_" 
                 << std::setfill('0') << std::setw(3) << (r + 1) << ".dat";
        
        std::ofstream file(filename.str());
        if (!file) {
            std::cerr << "Error: Cannot create seismogram file " << filename.str() << std::endl;
            continue;
        }
        
        for (Index it = 0; it < m_nstep; ++it) {
            const Real time = it * dt;
            file << std::scientific << std::setprecision(6) 
                 << time << " " << static_cast<float>(m_seismograms_vy[r][it]) << std::endl;
        }
        file.close();
    }
    
    std::cout << "Seismograms saved in Fortran-style ASCII format (" 
              << m_receivers.size() << " receivers, " << m_nstep << " time steps)" << std::endl;
}

void SeismogramRecorder::write_seismograms_ascii(const std::vector<std::vector<Real>>& seismograms, 
                                                 const std::string& filename, Real dt) const {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error: Cannot create seismogram file " << filename << std::endl;
        return;
    }
    
    // 写入头部信息
    file << "# Seismograms - Time vs Amplitude for " << m_receivers.size() << " receivers" << std::endl;
    file << "# Time step: " << dt << " seconds" << std::endl;
    file << "# Number of time steps: " << m_nstep << std::endl;
    
    // 写入数据：每行为时间步，每列为接收器
    for (Index it = 0; it < m_nstep; ++it) {
        const Real time = it * dt;
        file << std::scientific << std::setprecision(6) << time;
        for (Size r = 0; r < m_receivers.size(); ++r) {
            file << " " << static_cast<float>(seismograms[r][it]);
        }
        file << std::endl;
    }
    file.close();
}

void SeismogramRecorder::write_seismograms_su(const std::vector<std::vector<Real>>& seismograms, 
                                              const std::string& filename, Real dt) const {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error: Cannot create SU seismogram file " << filename << std::endl;
        return;
    }
    
    // 创建SU格式头部 (240字节 = 120个short整数)
    std::vector<int16_t> header(120, 0);
    header[57] = static_cast<int16_t>(m_nstep);        // ns: 样本数
    header[58] = static_cast<int16_t>(dt * 1e6);       // dt: 时间间隔 (微秒)
    
    // 按照原始Fortran格式：每个接收器一个trace
    for (Size r = 0; r < m_receivers.size(); ++r) {
        // 写入头部
        file.write(reinterpret_cast<const char*>(header.data()), 240);
        
        // 写入数据：转换为float并写入
        for (Index it = 0; it < m_nstep; ++it) {
            const float sample = static_cast<float>(seismograms[r][it]);
            file.write(reinterpret_cast<const char*>(&sample), sizeof(float));
        }
    }
    
    file.close();
    std::cout << "Seismograms saved in SU binary format (" 
              << m_receivers.size() << " traces, " << m_nstep << " samples each)" << std::endl;
}

} // namespace seismic
