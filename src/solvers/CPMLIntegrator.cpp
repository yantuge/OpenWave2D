#include "seismic/solvers/CPMLIntegrator.hpp"
#include "seismic/sources/SourceTimeFunction.hpp"
#include <iostream>

namespace seismic {

void CPMLIntegrator::step(Grid& grid, 
                         const SourceConfig& source_config,
                         const PMLConfig& pml_config, 
                         Real dt, 
                         Index it) {
    
    const Real current_time = it * dt;
    
    // 计算震源在网格中的位置索引
    const Index isource = static_cast<Index>(source_config.x / grid.deltax());
    const Index jsource = static_cast<Index>(source_config.y / grid.deltay());
    
    // 计算当前时间步的震源力
    Real force_x, force_y;
    calculate_source_forces_at_time(current_time, source_config, force_x, force_y);
    
    // 调用CPML各向异性Fortran内核
    // CPML只需要一次内核调用，不像RK4需要4次
    cpml_anisotropic_step_c(
        grid.fortran_nx(), grid.fortran_ny(), 
        grid.deltax(), grid.deltay(), dt, it,
        isource, jsource, force_x, force_y,
        // 各向异性材料参数
        grid.get_c11_ptr(), grid.get_c12_ptr(), 
        grid.get_c22_ptr(), grid.get_c33_ptr(), grid.get_density_ptr(),
        // 波场变量
        grid.get_vx_ptr(), grid.get_vy_ptr(),
        grid.get_sigmaxx_ptr(), grid.get_sigmayy_ptr(), grid.get_sigmaxy_ptr(),
        // CPML参数
        grid.get_d_x_ptr(), grid.get_d_x_half_ptr(),
        grid.get_d_y_ptr(), grid.get_d_y_half_ptr(),
        grid.get_k_x_ptr(), grid.get_k_x_half_ptr(),
        grid.get_k_y_ptr(), grid.get_k_y_half_ptr(),
        grid.get_alpha_x_ptr(), grid.get_alpha_x_half_ptr(),
        grid.get_alpha_y_ptr(), grid.get_alpha_y_half_ptr(),
        // CPML系数 (1D数组)
        grid.get_a_x_ptr(), grid.get_a_x_half_ptr(),
        grid.get_a_y_ptr(), grid.get_a_y_half_ptr(),
        grid.get_b_x_ptr(), grid.get_b_x_half_ptr(),
        grid.get_b_y_ptr(), grid.get_b_y_half_ptr(),
        // PML内存变量
        grid.get_memory_dvx_dx_ptr(), grid.get_memory_dvx_dy_ptr(),
        grid.get_memory_dvy_dx_ptr(), grid.get_memory_dvy_dy_ptr(),
        grid.get_memory_dsigmaxx_dx_ptr(), grid.get_memory_dsigmayy_dy_ptr(),
        grid.get_memory_dsigmaxy_dx_ptr(), grid.get_memory_dsigmaxy_dy_ptr()
    );
}

void CPMLIntegrator::calculate_source_forces_at_time(Real time,
                                                    const SourceConfig& source_config,
                                                    Real& force_x, 
                                                    Real& force_y) const {
    // 计算指定时间的震源幅度
    const Real source_amplitude = calculate_source_amplitude(time, source_config);
    
    // 根据震源角度分解为x和y分量
    calculate_source_forces(source_amplitude, source_config, force_x, force_y);
}

} // namespace seismic
