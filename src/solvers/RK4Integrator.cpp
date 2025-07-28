#include "seismic/solvers/RK4Integrator.hpp"
#include "seismic/sources/SourceTimeFunction.hpp"
#include <iostream>

namespace seismic {

// 静态成员定义
constexpr Real RK4Integrator::RK41_COEFFS[4];
constexpr Real RK4Integrator::RK42_COEFFS[4];

void RK4Integrator::step(Grid& grid, 
                        const SourceConfig& source_config,
                        const PMLConfig& pml_config, 
                        Real dt, 
                        Index it) {
    
    const Real current_time = it * dt;
    
    // 计算震源在网格中的位置索引
    const Index isource = static_cast<Index>(source_config.x / grid.deltax());
    const Index jsource = static_cast<Index>(source_config.y / grid.deltay());
    
    // 备份当前状态以进行RK4计算
    // 这是RK4算法的第一步：保存初始状态
    grid.backup_rk4_state();
    
    // 经典的4阶段RK4循环
    // 每个子步骤在不同的时间点计算导数，然后组合得到高精度的结果
    for (int rk_substep = 1; rk_substep <= 4; ++rk_substep) {
        // 计算当前RK子步骤的时间点
        const Real substep_time = current_time + RK41_COEFFS[rk_substep - 1] * dt;
        
        // 计算在该时间点的震源力
        Real force_x, force_y;
        calculate_source_forces_at_time(substep_time, source_config, force_x, force_y);
        
        // 调用ADE-PML RK4 Fortran内核
        // 注意：我们使用grid.fortran_nx()和grid.fortran_ny()来获取正确的网格尺寸
        adepml_rk4_step_c(
            grid.fortran_nx(), grid.fortran_ny(), 
            grid.deltax(), grid.deltay(), dt, rk_substep,
            isource, jsource, force_x, force_y,
            // 材料参数
            grid.get_lambda_ptr(), grid.get_mu_ptr(), grid.get_density_ptr(),
            // 波场变量
            grid.get_vx_ptr(), grid.get_vy_ptr(),
            grid.get_sigmaxx_ptr(), grid.get_sigmayy_ptr(), grid.get_sigmaxy_ptr(),
            // RK4中间变量
            grid.get_dvx_ptr(), grid.get_dvy_ptr(),
            grid.get_dsigmaxx_ptr(), grid.get_dsigmayy_ptr(), grid.get_dsigmaxy_ptr(),
            // PML参数
            grid.get_d_x_ptr(), grid.get_d_x_half_ptr(),
            grid.get_d_y_ptr(), grid.get_d_y_half_ptr(),
            grid.get_k_x_ptr(), grid.get_k_x_half_ptr(),
            grid.get_k_y_ptr(), grid.get_k_y_half_ptr(),
            grid.get_alpha_x_ptr(), grid.get_alpha_x_half_ptr(),
            grid.get_alpha_y_ptr(), grid.get_alpha_y_half_ptr(),
            // ADE-PML RK4特定系数
            grid.get_a_x_rk4_ptr(), grid.get_a_x_half_rk4_ptr(),
            grid.get_a_y_rk4_ptr(), grid.get_a_y_half_rk4_ptr(),
            grid.get_b_x_rk4_ptr(), grid.get_b_x_half_rk4_ptr(),
            grid.get_b_y_rk4_ptr(), grid.get_b_y_half_rk4_ptr(),
            // PML内存变量
            grid.get_memory_dvx_dx_ptr(), grid.get_memory_dvx_dy_ptr(),
            grid.get_memory_dvy_dx_ptr(), grid.get_memory_dvy_dy_ptr(),
            grid.get_memory_dsigmaxx_dx_ptr(), grid.get_memory_dsigmayy_dy_ptr(),
            grid.get_memory_dsigmaxy_dx_ptr(), grid.get_memory_dsigmaxy_dy_ptr()
        );
    }
    
    // RK4循环结束，此时grid中包含了下一时间步的解
}

void RK4Integrator::calculate_source_forces_at_time(Real time,
                                                   const SourceConfig& source_config,
                                                   Real& force_x, 
                                                   Real& force_y) const {
    // 计算指定时间的震源幅度
    const Real source_amplitude = calculate_source_amplitude(time, source_config);
    
    // 根据震源角度分解为x和y分量
    calculate_source_forces(source_amplitude, source_config, force_x, force_y);
}

} // namespace seismic
