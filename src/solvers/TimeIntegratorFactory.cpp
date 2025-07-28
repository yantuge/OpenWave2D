#include "seismic/solvers/TimeIntegratorFactory.hpp"
#include "seismic/solvers/RK4Integrator.hpp"
#include "seismic/solvers/CPMLIntegrator.hpp"
#include <stdexcept>
#include <sstream>

namespace seismic {

std::unique_ptr<TimeIntegrator> TimeIntegratorFactory::create(const TimeConfig& time_config,
                                                             const PMLConfig& pml_config) {
    
    // 根据PML类型选择积分器
    // 这里的策略是：PML类型决定了积分器类型，而不是时间方案
    if (pml_config.type == "cpml") {
        // CPML使用显式积分
        return std::make_unique<CPMLIntegrator>();
    } 
    else if (pml_config.type == "adepml") {
        // ADE-PML目前只支持RK4积分
        if (time_config.scheme == "rk4") {
            return std::make_unique<RK4Integrator>();
        } else {
            throw std::runtime_error(
                "ADE-PML only supports RK4 time integration scheme. " 
                "Got: " + time_config.scheme + ". " +
                "Please set time.scheme to 'rk4' in your configuration."
            );
        }
    } 
    else {
        throw std::runtime_error(
            "Unknown PML type: " + pml_config.type + ". " +
            "Supported types: cpml, adepml. " +
            get_supported_types()
        );
    }
}

bool TimeIntegratorFactory::is_supported(const TimeConfig& time_config,
                                        const PMLConfig& pml_config) {
    try {
        // 尝试创建积分器，如果成功说明支持
        auto integrator = create(time_config, pml_config);
        return true;
    } catch (const std::runtime_error&) {
        return false;
    }
}

std::string TimeIntegratorFactory::get_supported_types() {
    std::ostringstream oss;
    oss << "Supported combinations:\n";
    oss << "  - PML type 'cpml' + any time scheme -> CPML Integrator\n";
    oss << "  - PML type 'adepml' + time scheme 'rk4' -> RK4 Integrator\n";
    oss << "\nNotes:\n";
    oss << "  - CPML is computationally efficient but less accurate\n";
    oss << "  - ADE-PML with RK4 is more accurate but requires more computation\n";
    return oss.str();
}

// 便利函数的实现
std::unique_ptr<TimeIntegrator> create_time_integrator(const TimeConfig& time_config,
                                                      const PMLConfig& pml_config) {
    return TimeIntegratorFactory::create(time_config, pml_config);
}

} // namespace seismic
