//
// Created by Andrew on 2025-10-17.
//

#ifndef JUNIPEREXE_LEAPFROGINTEGRATOR_H
#define JUNIPEREXE_LEAPFROGINTEGRATOR_H
#include "Integrator.h"
#include "Simulation.h"

class LeapfrogIntegrator : public Integrator {

    Simulation &sim;

    Point3f accForParticle(int particle, TreeNode& leaf);
    float energyChangeForParticle(int particle, TreeNode& leaf);
    void step(SimData& data, float timestep) override;

public:
    explicit LeapfrogIntegrator(Simulation& sim);
};


#endif //JUNIPEREXE_LEAPFROGINTEGRATOR_H