#include <filesystem>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "Simulation.h"

void assertTestSim(Simulation& sim) {
    ASSERT_THAT(sim.getSimData().xyzh, testing::ElementsAre(3, 5, 2, 0.5, 1, 5, 3, 2, 2, 20, 1, 0.7, 5, 10, 2, 0.1));
    ASSERT_THAT(sim.getSimData().vxyzu, testing::ElementsAre(1.5, 2.5, 0.25, 0.3, 4, 1, -0.5, 0.3, 0.5, 0.5, 0.5, 0.3, 0.5, 0.5, 0.5, 0.3));
    ASSERT_THAT(sim.getSimData().fxyz, testing::ElementsAre(1, 2, 3, 3, 2.3, 7.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5));
    ASSERT_EQ(sim.xmin, 1);
    ASSERT_EQ(sim.xmax, 5);
    ASSERT_EQ(sim.ymin, 5);
    ASSERT_EQ(sim.ymax, 20);
    ASSERT_EQ(sim.zmin, 1);
    ASSERT_EQ(sim.zmax, 3);
    ASSERT_TRUE(sim.getSimData().doesContainEnergy());
}

TEST(IOTest, BasicReaderTest) {
    Simulation data = Simulation(".//files//test.csv");

    assertTestSim(data);
}

TEST(IOTest, BasicWriterTest) {
    Simulation sim = Simulation(".//files//test.csv");
    sim.getSimData().toCSV(".//files//temp_writing_test.csv");
    Simulation writtenSim = Simulation(".//files//temp_writing_test.csv");
    try {
        std::filesystem::remove(".//files//temp_writing_test.csv");
    } catch (const std::filesystem::filesystem_error& err) {
        std::cout << "Error deleting temporary file!: " << err.what() << "\n";
    }

    assertTestSim(writtenSim);
}

