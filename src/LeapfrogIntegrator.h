//
// Created by Andrew on 2025-10-17.
//

#ifndef JUNIPEREXE_LEAPFROGINTEGRATOR_H
#define JUNIPEREXE_LEAPFROGINTEGRATOR_H
#include "Integrator.h"


class LeapfrogIntegrator : public Integrator {
    void step(SimData& data, std::vector<float> accs, float timestep) override;
};


#endif //JUNIPEREXE_LEAPFROGINTEGRATOR_H