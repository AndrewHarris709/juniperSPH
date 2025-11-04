#include <gtest/gtest.h>
#include "Simulation.h"

void assertAtDepth(TreeNode* node, int depth, int partCountAtDepth) {
    if (depth > 0) {
        assertAtDepth(node->getLeftChild().get(), depth - 1, partCountAtDepth);
        assertAtDepth(node->getRightChild().get(), depth - 1, partCountAtDepth);
    } else {
        ASSERT_EQ(node->getParticleCount(), partCountAtDepth);
        ASSERT_TRUE(node->isLeaf());
    }
}


TEST(NeighboursTest, TreeBuildTest) {
    Simulation sim = Simulation(".//files//kd_test.csv");
    sim.buildTree();

    assertAtDepth(&sim.getBaseNode(), 6, 5);
}
