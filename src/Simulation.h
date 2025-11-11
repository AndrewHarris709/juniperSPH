#pragma once
#include <vector>
#include <string>
#include <unordered_map>

#include "Integrator.h"
#include "kernel.h"
#include "SimData.h"
#include "TreeNode.h"

#ifndef SIMULATION_H
#define SIMULATION_H

class Simulation {
    // The order of these arguments determines the initialize order in the
    // constructor. So don't change this order!
    SimData simData;
    ParticleSet globalSet;
    TreeNode baseNode;
    Kernel kernel;
    Integrator* integrator;
    std::vector<TreeNode*> leaves;
    std::unordered_map<int, std::vector<int>> neighbourCache;

    float distBetween(float x1, float x2, float y1, float y2, float z1, float z2) const;
    float densityIterationForParticle(int particle, TreeNode& node);

public:
    explicit Simulation(const std::string& filename);
    explicit Simulation() : Simulation("") {};
    ~Simulation();

    float xmin, xmax, ymin, ymax, zmin, zmax;

    void useConfig(const std::string& filename);
    void stepSimulation();

    float distBetween(int part1, int part2) const;
    float distBetweenNodes(TreeNode& node1, TreeNode& node2) const;
    void densityIterate();
    std::vector<int> getNeighbours(int part);

    Point3f displacementBetween(float x1, float x2, float y1, float y2, float z1, float z2) const;

    std::vector<int> getNeighbours(int target, TreeNode& node, bool isCache);
    void resetNeighbourCache();

    float densityAt(int part, TreeNode& node);
    float pressureAt(int part);
    float omegaAt(int part, TreeNode& node);
    Point3f velocityDiffBetween(int target, int part);
    Point3f displacementBetween(int target, int part);

    void setLimits();
    void setLimits(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax);
    int getParticleCount() const;

    SimData& getSimData();
    TreeNode& getBaseNode();
    void buildTree();
    Kernel& getKernel();
    std::vector<TreeNode*>& getLeaves();

    // for artificial viscosity
    float qabAt(int partA, int partB);
    static Point3f norm(Point3f p);
    static float dot(Point3f p1, Point3f p2);
};

#endif //SIMULATION_H
