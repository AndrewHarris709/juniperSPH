//
// Created by Andrew on 2025-10-17.
//

#include "LeapfrogIntegrator.h"

void addToParticle(int particle, std::vector<float>& data1, int width1, float multiple, std::vector<float>& data2, int width2, float timestep) {
    // Width here is only for xyzh, to make sure that the particles are indexed correctly.
    // TODO: There should be a better way to link continguous data in vectors with width.
    for (int axis = 0; axis < 3; axis++) {
        data1[particle * width1 + axis] += multiple * timestep * data2[particle * width2 * axis];
    }
}

void LeapfrogIntegrator::step(SimData& data, std::vector<float> accs, float timestep) {
    for (int i = 0; i < data.getParticleCount(); i++) {

        addToParticle(i, data.vxyzv, 3, 0.5, accs, 3, timestep);
        addToParticle(i, data.xyzh, 4, 1, data.vxyzv, 3, timestep);

        accs[i*3 + 0] *= data.xyzh[i*3 + 0];
        accs[i*3 + 1] *= data.xyzh[i*3 + 1];
        accs[i*3 + 2] *= data.xyzh[i*3 + 2];

        addToParticle(i, data.vxyzv, 3, 0.5, accs, 3, timestep);
    }
}
