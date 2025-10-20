#include <iostream>
#include <string>

#include "Simulation.h"

int main(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " particle_file config_file output_file" << std::endl;
    }
    Simulation data = Simulation(argv[1]);
    data.useConfig(argv[2]);

    data.densityIterate();

    std::cout << "Density iteration complete" << std::endl;

    data.stepSimulation();

    std::cout << "Simulation Step Complete" << std::endl;

    // This argument should be verified prior to simulation.
    data.getSimData().toCSV(argv[3]);

    return 0;
}
