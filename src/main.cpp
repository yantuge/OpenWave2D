#include <iostream>
#include <fstream>
#include <stdexcept>
#include <nlohmann/json.hpp>
#include "seismic/main/Dispatcher.hpp"
#include "seismic/core/Config.hpp"

int main(int argc, char* argv[]) {
    try {
        // Parse command line arguments
        if (argc != 2) {
            std::cerr << "Usage: " << argv[0] << " <config.json>" << std::endl;
            return 1;
        }
        
        const std::string config_filename = argv[1];
        
        // Read JSON configuration file
        std::ifstream config_file(config_filename);
        if (!config_file) {
            std::cerr << "Error: Cannot open configuration file: " << config_filename << std::endl;
            return 1;
        }
        
        nlohmann::json j;
        config_file >> j;
        config_file.close();
        
        // Parse configuration
        auto config = seismic::SimulationConfig::from_json(j);
        
        // Run simulation
        seismic::run_simulation(config);
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred" << std::endl;
        return 1;
    }
}
