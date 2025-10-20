//
// Created by Andrew on 2025-09-16.
//

#include "SimData.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ranges>

SimData::SimData(const std::string& filename) {
    time = 0.0;
    m = 0.01;

    if (filename.empty()) {
        return;
    }

    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cout << std::filesystem::current_path().string() << std::endl;
        std::cerr << "Error opening file " << filename << std::endl;
        return;
    }

    std::vector<std::vector<std::string>> data;
    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string cell;

        while (std::getline(ss, cell, ',')) {
            // Linux systems will sometimes return a '/r' character from getline.
            std::erase(cell, '\r' );
            row.push_back(cell);
        }

        data.push_back(row);
    }

    for (const auto& row : data | std::views::drop(1)) {
        int i = 0;

        for (const auto& cell : row) {
            std::string column = data[0][i];
            if (std::ranges::find(posCols, column) != std::end(posCols)) {
                xyzh.push_back(std::stof(cell));
            }
            if (std::ranges::find(velCols, column) != std::end(velCols)) {
                if (column == "u") {
                    containsEnergy = true;
                }
                vxyzu.push_back(std::stof(cell));
            }
            if (std::ranges::find(varCols, column) != std::end(varCols)) {
                fxyz.push_back(std::stof(cell));
            }

            i++;
        }
    }

    file.close();
}

void SimData::toCSV(const std::string &filename) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cout << std::filesystem::current_path().string() << std::endl;
        std::cerr << "Error writing to file " << filename << std::endl;
    }

    file << "x,y,z,h,vx,vy,vz,u,fx,fy,fz" << std::endl;

    for (int i = 0; i < getParticleCount(); i++) {
        file << xyzh[4*i] << "," << xyzh[4*i+1] << "," << xyzh[4*i+2] << "," << xyzh[4*i+3] << ",";
        file << vxyzu[4*i] << "," << vxyzu[4*i+1] << "," << vxyzu[4*i+2] << "," << vxyzu[4*i+3] << ",";
        file << fxyz[3*i] << "," << fxyz[3*i+1] << "," << fxyz[3*i+2];
        file << std::endl;
    }

    file.close();
}

int SimData::getParticleCount() const {
    return std::floor(xyzh.size() / 4);
}

bool SimData::doesContainEnergy() const {
    return this->containsEnergy;
}

