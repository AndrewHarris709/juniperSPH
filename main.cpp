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

    for (int step = 1; step < 100; step++) {
        std::cout << "Simulation Step " << step << "...";
        data.stepSimulation();
        std::cout << " Complete...";
        data.getSimData().toCSV("./dumps/step_" + std::to_string(step) + ".csv");
        std::cout << "Saved." << std::endl;
        data.densityIterate();
    }

    std::cout << "Simulation Steps Complete" << std::endl;

    // This argument should be verified prior to simulation.
    data.getSimData().toCSV(argv[3]);

    return 0;
}
