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
    
    // 计算震源在网格中的位置索引 (转换为Fortran 1-based索引)
    const Index isource = static_cast<Index>(source_config.x / grid.deltax()) + 1;
    const Index jsource = static_cast<Index>(source_config.y / grid.deltay()) + 1;
    
    // 备份当前状态以进行RK4计算
    // 这是RK4算法的关键：保存初始状态
    grid.backup_rk4_state();
    
    // RK4积分：k1, k2, k3, k4四个子步骤
    // 每个子步骤计算在不同时间点的导数，最后进行加权平均
    for (int rk_substep = 1; rk_substep <= 4; ++rk_substep) {
        // 计算当前RK子步骤的时间点
        const Real substep_time = current_time + RK41_COEFFS[rk_substep - 1] * dt;
        
        // 计算在该时间点的震源力
        Real force_x, force_y;
        calculate_source_forces_at_time(substep_time, source_config, force_x, force_y);
        
        // 调用ADE-PML Fortran内核计算该子步骤的导数
        // 注意：Fortran内核现在只计算导数，不更新波场
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
        
        // 执行RK4积分：将计算出的导数累加到波场中
        // 使用RK42系数进行正确的加权
        perform_rk4_accumulation(grid, source_config, dt, rk_substep, force_x, force_y);
    }
    
    // RK4循环结束，此时grid中包含了正确的下一时间步解
}

void RK4Integrator::perform_rk4_accumulation(Grid& grid, 
                                            const SourceConfig& source_config,
                                            Real dt, 
                                            int rk_substep,
                                            Real force_x, 
                                            Real force_y) const {
    
    // 获取RK4系数
    const Real rk4_coeff = RK42_COEFFS[rk_substep - 1];
    const Real dt_rk = dt * rk4_coeff;
    
    // 计算震源位置
    const Index isource = static_cast<Index>(source_config.x / grid.deltax());
    const Index jsource = static_cast<Index>(source_config.y / grid.deltay());
    
    // 执行RK4积分：将导数累加到波场变量
    for (Index j = 0; j < grid.ny(); ++j) {
        for (Index i = 0; i < grid.nx(); ++i) {
            // 获取当前时间步开始时的备份状态 (stored in dvx(3,:,:), etc.)
            // 和当前子步骤计算出的导数 (stored in dvx(1,:,:), etc.)
            
            // 注意：我们需要通过raw pointer访问RK4中间数组
            // 数组布局：dvx[3 * (ny_extended * nx_extended)]
            // 其中：dvx[0 * size + index] = dvx(1,i,j) - 当前导数
            //      dvx[2 * size + index] = dvx(3,i,j) - 备份状态
            
            const Index nx_ext = grid.nx() + 8;  // extended grid: -4 to nx+4
            const Index ny_ext = grid.ny() + 8;
            const Index i_ext = i + 4;  // convert to extended grid index
            const Index j_ext = j + 4;
            const Index index = j_ext * nx_ext + i_ext;
            const Index grid_size = nx_ext * ny_ext;
            
            Real* dvx_ptr = grid.get_dvx_ptr();
            Real* dvy_ptr = grid.get_dvy_ptr();
            Real* dsigmaxx_ptr = grid.get_dsigmaxx_ptr();
            Real* dsigmayy_ptr = grid.get_dsigmayy_ptr();
            Real* dsigmaxy_ptr = grid.get_dsigmaxy_ptr();
            
            // RK4积分：从备份状态开始，累加加权导数
            // vx(i,j) = vx_backup + dt_rk * dvx/dt
            grid.vx(i, j) = dvx_ptr[2 * grid_size + index] + dt_rk * dvx_ptr[0 * grid_size + index];
            grid.vy(i, j) = dvy_ptr[2 * grid_size + index] + dt_rk * dvy_ptr[0 * grid_size + index];
            grid.sigmaxx(i, j) = dsigmaxx_ptr[2 * grid_size + index] + dt_rk * dsigmaxx_ptr[0 * grid_size + index];
            grid.sigmayy(i, j) = dsigmayy_ptr[2 * grid_size + index] + dt_rk * dsigmayy_ptr[0 * grid_size + index];
            grid.sigmaxy(i, j) = dsigmaxy_ptr[2 * grid_size + index] + dt_rk * dsigmaxy_ptr[0 * grid_size + index];
        }
    }
    
    // 添加震源项到速度场
    if (isource >= 0 && isource < grid.nx() && jsource >= 0 && jsource < grid.ny()) {
        const Real rho_source = grid.density(isource, jsource);
        grid.vx(isource, jsource) += force_x * dt_rk / rho_source;
        grid.vy(isource, jsource) += force_y * dt_rk / rho_source;
    }
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
