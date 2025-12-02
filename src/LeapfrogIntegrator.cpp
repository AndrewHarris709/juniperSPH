//
// Created by Andrew on 2025-10-17.
//

#include "LeapfrogIntegrator.h"

LeapfrogIntegrator::LeapfrogIntegrator(Simulation &sim) : sim(sim) { }

void addToParticle(int particle, std::vector<float>& data1, int width1, float multiple, std::vector<float>& data2, int width2, float timestep) {
    // Width here is for xyzh & vxyzu when u is present, to make sure that the particles are indexed correctly.
    // TODO: There should be a better way to link continguous data in vectors with width.
    for (int axis = 0; axis < 3; axis++) {
        data1[particle * width1 + axis] += multiple * timestep * data2[particle * width2 + axis];
    }
}

void LeapfrogIntegrator::periodicCorrection(int particle) {
    SimData& data = sim.getSimData();

    float x = data.xyzh[particle * 4 + 0], y = data.xyzh[particle * 4 + 1], z = data.xyzh[particle * 4 + 2];

    if (x > sim.xmax) {
        x -= (sim.xmax - sim.xmin);
    } else if (x < sim.xmin) {
        x += (sim.xmax - sim.xmin);
    }

    if (y > sim.ymax) {
        y -= (sim.ymax - sim.ymin);
    } else if (y < sim.ymin) {
        y += (sim.ymax - sim.ymin);
    }

    if (z > sim.zmax) {
        z -= (sim.zmax - sim.zmin);
    } else if (z < sim.zmin) {
        z += (sim.zmax - sim.zmin);
    }

    data.xyzh[particle * 4 + 0] = x;
    data.xyzh[particle * 4 + 1] = y;
    data.xyzh[particle * 4 + 2] = z;
}

void LeapfrogIntegrator::precompute() {
    SimData& simData = sim.getSimData();

    for (TreeNode* node : sim.getLeaves()) {
        for (int part: node->getParticleIndices()) {
            simData.density[part] = sim.densityAt(part, *node);
            simData.omega[part] = sim.omegaAt(part, *node);
        }
    }
}

Point3f LeapfrogIntegrator::accForParticle(int target, TreeNode& leaf) {
    float ax = 0, ay = 0, az = 0;
    SimData& data = sim.getSimData();
    Kernel& kernel = sim.getKernel();

    float targetDensity = data.density[target];
    float targetOmega = data.omega[target];

    std::vector<int> neighbours = sim.getNeighbours(target, leaf, true);
    for (int part: neighbours) {
        float targetQAB = sim.qabAt(target, part);
        float targetRatio = (sim.pressureAt(target) + targetQAB) / (targetDensity * targetDensity * targetOmega);

        float partDensity = data.density[part];
        float partOmega = data.omega[part];
        float partQAB = sim.qabAt(part, target);
        float partRatio = (sim.pressureAt(part) + partQAB) / (partDensity * partDensity * partOmega);

        Point3f disp = sim.displacementBetween(target, part);
        Point3f dispNorm = sim.norm(disp);
        float distance = sim.distBetween(target, part);
        float hTarget = data.xyzh[4 * target + 3], hPart = data.xyzh[4 * part + 3];

        float targetGradient = kernel.gradientAt(distance / hTarget) / (hTarget * hTarget * hTarget * hTarget);
        float partGradient = kernel.gradientAt(distance / hPart) / (hPart * hPart * hPart * hPart);

        ax -= data.m * (targetRatio * dispNorm.x * targetGradient + partRatio * dispNorm.x * partGradient);
        ay -= data.m * (targetRatio * dispNorm.y * targetGradient + partRatio * dispNorm.y * partGradient);
        az -= data.m * (targetRatio * dispNorm.z * targetGradient + partRatio * dispNorm.z * partGradient);
    }

    return Point3f{ax, ay, az};
}

float LeapfrogIntegrator::energyChangeForParticle(int target, TreeNode& leaf) {
    SimData& data = sim.getSimData();
    Kernel& kernel = sim.getKernel();

    float targetPressure = sim.pressureAt(target);
    float targetDensity = data.density[target];
    float targetOmega = data.omega[target];

    float targetRatio = targetPressure / (targetDensity * targetDensity * targetOmega);
    std::vector<int> neighbours = sim.getNeighbours(target, leaf, true);

    float uChange = 0;
    for (int part: neighbours) {
        Point3f vDiff = sim.velocityDiffBetween(target, part);
        Point3f disp = sim.displacementBetween(target, part);
        Point3f dispNorm = sim.norm(disp);
        float distance = sim.distBetween(target, part);

        float hTarget = data.xyzh[4 * target + 3];
        float gradient = kernel.gradientAt(distance / hTarget) / (hTarget * hTarget * hTarget * hTarget);

        float xComp = vDiff.x * dispNorm.x * gradient;
        float yComp = vDiff.y * dispNorm.y * gradient;
        float zComp = vDiff.z * dispNorm.z * gradient;

        uChange += data.m * (xComp + yComp + zComp);
    }

    return uChange * targetRatio;
}

void LeapfrogIntegrator::computeAccsAndEnergies() {
    precompute();
    SimData& data = sim.getSimData();
    std::vector<TreeNode*>& leaves = sim.getLeaves();

    for (TreeNode* leaf : leaves) {
        for (int i : leaf->getParticleIndices()) {
            Point3f iAccs = accForParticle(i, *leaf);
            data.accs[3 * i] = iAccs.x;
            data.accs[3 * i + 1] = iAccs.y;
            data.accs[3 * i + 2] = iAccs.z;
            data.energies[i] = energyChangeForParticle(i, *leaf);
        }
    }
}

void LeapfrogIntegrator::step(SimData& data, float timestep) {
    int vWidth = 3;

    if (firstStep) {
        computeAccsAndEnergies();
        firstStep = false;
    }

    // TODO: Huge code smell here, the cache should be replaced within Simulation
    // as soon as a particle's details changes, not externally like this.
    sim.resetNeighbourCache();

    for (int i = 0; i < sim.getParticleCount(); i++) {
        addToParticle(i, data.vxyzu, vWidth, 0.5, data.accs, 3, timestep);
        addToParticle(i, data.xyzh, 4, 1, data.vxyzu, vWidth, timestep);
        periodicCorrection(i);
        data.vxyzu[i * 4 + 3] += 0.5 * data.energies[i] * timestep;
    }

    sim.getSimData().toCSV("first_kick.csv");

    sim.buildTree();
    computeAccsAndEnergies();

    for (int i = 0; i < sim.getParticleCount(); i++) {
        addToParticle(i, data.vxyzu, vWidth, 0.5, data.accs, 3, timestep);
        data.vxyzu[i * 4 + 3] += 0.5 * data.energies[i] * timestep;
    }

    // TODO: See above.
    sim.resetNeighbourCache();
}
