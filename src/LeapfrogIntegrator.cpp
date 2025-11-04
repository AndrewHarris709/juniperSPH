//
// Created by Andrew on 2025-10-17.
//

#include "LeapfrogIntegrator.h"

LeapfrogIntegrator::LeapfrogIntegrator(Simulation &sim) : sim(sim) { }

void addToParticle(int particle, std::vector<float>& data1, int width1, float multiple, std::vector<float>& data2, int width2, float timestep) {
    // Width here is for xyzh & vxyzu when u is present, to make sure that the particles are indexed correctly.
    // TODO: There should be a better way to link continguous data in vectors with width.
    for (int axis = 0; axis < 3; axis++) {
        data1[particle * width1 + axis] += multiple * timestep * data2[particle * width2 * axis];
    }
}

Point3f LeapfrogIntegrator::accForParticle(int target, TreeNode& leaf) {
    float ax = 0, ay = 0, az = 0;
    SimData& data = sim.getSimData();
    Kernel& kernel = sim.getKernel();

    std::vector<int> neighbours = sim.getNeighbours(target, leaf, true);
    for (int part: leaf.getParticleIndices()) {
        float targetDensity = sim.densityAt(target, leaf);
        float targetOmega = sim.omegaAt(target, leaf);
        float targetQAB = sim.qabAt(target, part, leaf);
        float targetRatio = (sim.pressureAt(target, leaf) + targetQAB) / (targetDensity * targetDensity * targetOmega);

        float partDensity = sim.densityAt(part, leaf);
        float partOmega = sim.omegaAt(target, leaf);
        float partQAB = sim.qabAt(part, target, leaf);
        float partRatio = (sim.pressureAt(part, leaf) + partQAB) / (partDensity * partDensity * partOmega);

        Point3f disp = sim.displacementBetween(target, part);
        float hTarget = data.xyzh[4 * target + 3], hPart = data.xyzh[4 * part + 3];

        ax += data.m * (targetRatio * kernel.gradientAt(disp.x) / hTarget + partRatio * kernel.gradientAt(disp.x) / hPart);
        ay += data.m * (targetRatio * kernel.gradientAt(disp.y) / hTarget + partRatio * kernel.gradientAt(disp.y) / hPart);
        az += data.m * (targetRatio * kernel.gradientAt(disp.z) / hTarget + partRatio * kernel.gradientAt(disp.z) / hPart);
    }

    return Point3f{ax, ay, az};
}

float LeapfrogIntegrator::energyChangeForParticle(int target, TreeNode& leaf) {
    SimData& data = sim.getSimData();
    Kernel& kernel = sim.getKernel();

    float targetPressure = sim.pressureAt(target, leaf);
    float targetDensity = sim.densityAt(target, leaf);
    float targetOmega = sim.omegaAt(target, leaf);

    float targetRatio = targetPressure / (targetDensity * targetDensity * targetOmega);
    std::vector<int> neighbours = sim.getNeighbours(target, leaf, true);

    float uChange = 0;
    for (int part: leaf.getParticleIndices()) {
        Point3f vDiff = sim.velocityDiffBetween(target, part);
        Point3f disp = sim.displacementBetween(target, part);
        float hTarget = data.xyzh[4 * target + 3];

        float xComp = vDiff.x * (kernel.gradientAt(disp.x) / hTarget);
        float yComp = vDiff.y * (kernel.gradientAt(disp.y) / hTarget);
        float zComp = vDiff.z * (kernel.gradientAt(disp.z) / hTarget);

        uChange += data.m * (xComp + yComp + zComp);
    }

    return uChange * targetRatio;
}

void LeapfrogIntegrator::step(SimData& data, float timestep) {
    std::vector<float> accs, energies;
    int vWidth = 3;
    std::vector<TreeNode*>& leaves = sim.getLeaves();

    if (data.doesContainEnergy()) {
        for (TreeNode* leaf : leaves) {
            for (int i : leaf->getParticleIndices()) {
                Point3f iAccs = accForParticle(i, *leaf);
                accs.push_back(iAccs.x);
                accs.push_back(iAccs.y);
                accs.push_back(iAccs.z);
                energies.push_back(energyChangeForParticle(i, *leaf));
            }
        }
        vWidth = 4;
    } else {
        accs.assign(sim.getParticleCount() * 3, 0);
        energies.assign(sim.getParticleCount(), 0);
    }

    // TODO: Huge code smell here, the cache should be replaced within Simulation
    // as soon as a particle's details changes, not externally like this.
    sim.resetNeighbourCache();

    for (TreeNode* leaf : leaves) {
        for (int i : leaf->getParticleIndices()) {
            addToParticle(i, data.vxyzu, vWidth, 0.5, accs, 3, timestep);
            addToParticle(i, data.xyzh, 4, 1, data.vxyzu, vWidth, timestep);
            data.vxyzu[i * 4 + 3] += 0.5 * energies[i] * timestep;
        }
    }

    sim.buildTree();

    for (TreeNode* leaf : leaves) {
        for (int i: leaf->getParticleIndices()) {
            Point3f newIAccs = accForParticle(i, *leaf);
            accs[i*3 + 0] = newIAccs.x;
            accs[i*3 + 1] = newIAccs.y;
            accs[i*3 + 2] = newIAccs.z;
            energies[i] = energyChangeForParticle(i, *leaf);

            addToParticle(i, data.vxyzu, vWidth, 0.5, accs, 3, timestep);
            data.vxyzu[i * 4 + 3] += 0.5 * energies[i] * timestep;
        }
    }

    // TODO: See above.
    sim.resetNeighbourCache();
}
