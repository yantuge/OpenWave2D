#pragma once

#include "TimeIntegrator.hpp"

namespace seismic {

/**
 * @class RK4Integrator
 * @brief 经典的四阶龙格-库塔时间积分器
 * 
 * 这个类实现了ADE-PML的RK4时间积分方案。RK4是一种四阶精度的
 * 显式时间积分方法，特别适用于波动方程的求解。
 * 
 * RK4的特点：
 * - 四个子步骤，每个子步骤使用不同的时间权重
 * - 高精度（四阶）但计算开销较大
 * - 需要额外的存储空间来保存中间状态
 */
class RK4Integrator : public TimeIntegrator {
public:
    /**
     * @brief 构造函数
     */
    RK4Integrator() = default;
    
    /**
     * @brief 析构函数
     */
    ~RK4Integrator() override = default;

    /**
     * @brief 执行RK4时间步积分
     * 
     * 实现经典的四阶段RK4算法：
     * 1. 备份当前状态（k0）
     * 2. 计算四个RK子步骤，每个子步骤使用不同的时间权重
     * 3. 最终组合得到下一时间步的解
     */
    void step(Grid& grid, 
             const SourceConfig& source_config,
             const PMLConfig& pml_config, 
             Real dt, 
             Index it) override;

    /**
     * @brief 获取积分器名称
     */
    std::string get_name() const override { return "RK4"; }
    
    /**
     * @brief 获取RK4的稳定性系数
     * RK4通常比显式差分更稳定，但仍需要满足CFL条件
     */
    Real get_stability_factor() const override { return 2.8; }

private:
    // RK4系数现在是这个类的私有成员，不再硬编码于主循环
    // 这些系数定义了RK4的四个子步骤的时间权重
    static constexpr Real RK41_COEFFS[4] = {0.0, 0.5, 0.5, 1.0};
    
    // RK42系数用于最终的加权平均（如果需要的话）
    static constexpr Real RK42_COEFFS[4] = {1.0/6.0, 2.0/6.0, 2.0/6.0, 1.0/6.0};
    
    /**
     * @brief 计算震源在指定时间的贡献
     * @param time 当前时间
     * @param source_config 震源配置
     * @param force_x 输出：X方向力分量
     * @param force_y 输出：Y方向力分量
     */
    void calculate_source_forces_at_time(Real time,
                                        const SourceConfig& source_config,
                                        Real& force_x, 
                                        Real& force_y) const;
    
    /**
     * @brief 执行RK4积分的累加步骤
     * @param grid 网格对象
     * @param source_config 震源配置
     * @param dt 时间步长
     * @param rk_substep RK4子步骤编号(1-4)
     * @param force_x X方向震源力
     * @param force_y Y方向震源力
     */
    void perform_rk4_accumulation(Grid& grid, 
                                 const SourceConfig& source_config,
                                 Real dt, 
                                 int rk_substep,
                                 Real force_x, 
                                 Real force_y) const;
};

} // namespace seismic
