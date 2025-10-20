//
// Created by Andrew on 2025-10-17.
//

#ifndef JUNIPEREXE_INTEGRATOR_H
#define JUNIPEREXE_INTEGRATOR_H
#include "SimData.h"


class Integrator {
public:
    virtual ~Integrator() = default;
    virtual void step(SimData& data, std::vector<float> accs, float timestep) = 0;
};


#endif //JUNIPEREXE_INTEGRATOR_H