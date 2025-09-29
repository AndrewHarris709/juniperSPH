#include <random>
#include <gtest/gtest.h>
#include "Simulation.h"
#include "kernel.h"

TEST(DensityTest, SimpleDensityTest) {
    // Note: This test currently depends on the default m=0.01.
    Simulation sim = Simulation();
    SimData& data = sim.getSimData();
    Kernel kernel = Kernel();

    int parts1d = 2;
    float L = 1.0f;
    float increment = L / parts1d;

    for (int x = 0; x < parts1d; x++) {
        for (int y = 0; y < parts1d; y++) {
            for (int z = 0; z < parts1d; z++) {
                data.xyzh.push_back(x * increment);
                data.xyzh.push_back(y * increment);
                data.xyzh.push_back(z * increment);
                // A large h value, so the distance between neighbouring particles is effectively zero.
                data.xyzh.push_back(100000000000000);
            }
        }
    }

    sim.setLimits(0, L, 0, L, 0, L);

    for (int part = 0; part < parts1d * parts1d * parts1d; part++) {
        ASSERT_FLOAT_EQ(sim.densityAt(part, kernel), 0.02546479089);
    }
}

TEST(DensityTest, PhantomDensityTest) {
    Simulation sim = Simulation("files/hydro32_00020.csv");
    SimData& data = sim.getSimData();
    Kernel kernel = Kernel();

    data.m = 3.0517578125e-05;

    // Even if we randomly perturb the h-values of some particles, we should still return to the same state.
    std::vector<float> oldxyzh(data.xyzh);
    std::default_random_engine el(15);
    std::uniform_real_distribution<float> distribution(0.9, 1.1);
    for (int i = 0; i < data.getParticleCount(); i++) {
        if (i % 4 == 0) {
            data.xyzh[4 * i + 3] = data.xyzh[4 * i + 3] * distribution(el);
        }
    }

    sim.densityIterate(kernel);
    for (int i = 0; i < data.getParticleCount(); i++) {
        EXPECT_NEAR(data.xyzh[4 * i + 3], oldxyzh[4 * i + 3], 0.0001);
    }
}

TEST(DensityTest, PeriodicBoundariesTest) {
    // Note: This test currently depends on the default m=0.01.
    Simulation sim = Simulation();
    SimData& data = sim.getSimData();
    Kernel kernel = Kernel();

    int parts = 2;
    float L = 1.0f;
    float increment = L / (parts - 1);

    for (int x = 0; x < parts; x++) {
        data.xyzh.push_back(x * increment);
        data.xyzh.push_back(0);
        data.xyzh.push_back(0);
        data.xyzh.push_back(0.25);
    }
    sim.setLimits(0, 1, 0, 0, 0, 0);

    ASSERT_NE(sim.densityAt(0, kernel), kernel.valueAt(0));
}
