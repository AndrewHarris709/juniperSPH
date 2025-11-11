//
// Created by Andrew on 2025-10-17.
//

#ifndef JUNIPEREXE_LEAPFROGINTEGRATOR_H
#define JUNIPEREXE_LEAPFROGINTEGRATOR_H
#include "Integrator.h"
#include "Simulation.h"

class LeapfrogIntegrator : public Integrator {

    Simulation &sim;

    void periodicCorrection(int particle);

    void precompute();

    Point3f accForParticle(int particle, TreeNode& leaf);
    float energyChangeForParticle(int particle, TreeNode& leaf);

    void computeAccsAndEnergies();
    void step(SimData& data, float timestep) override;

    bool firstStep = true;

public:
    explicit LeapfrogIntegrator(Simulation& sim);
};


#endif //JUNIPEREXE_LEAPFROGINTEGRATOR_H