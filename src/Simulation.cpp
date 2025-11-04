//
// Created by Andrew on 2025-07-01.
//

#include "Simulation.h"
#include "toml.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <ranges>
#include <limits>
#include <stack>

#include "LeapfrogIntegrator.h"

constexpr int NODE_PARTICLE_MIN = 10;
constexpr int MAX_DENSITY_ITERATIONS = 400;
constexpr float GAMMA = 5.0 / 3;

Simulation::Simulation(const std::string& filename) : simData(filename), globalSet(simData),
                                                      baseNode(nullptr, globalSet) {
    integrator = new LeapfrogIntegrator(*this);

    // Do not set limits if data is not supplied, like in unit tests.
    if (filename.empty()) {
        return;
    }
    this->setLimits();
}

Simulation::~Simulation() {
    delete integrator;
}

void Simulation::useConfig(const std::string& filename) {
    toml::table tbl;
    try
    {
        tbl = toml::parse_file(filename);
        this->simData.m = tbl["simconfig"]["mass"].value<float>().value();
    }
    catch (const toml::parse_error& err)
    {
        std::cerr << "Parsing failed:\n" << err << "\n";
    }
}

void Simulation::buildTree() {
    leaves.clear();
    TreeNode newNode(nullptr, this->globalSet);
    baseNode = newNode;

    std::stack<TreeNode*> nodeStack;
    nodeStack.push(&this->baseNode);

    while (!nodeStack.empty()) {
        TreeNode* node = nodeStack.top();
        nodeStack.pop();
        if (node->getParticleCount() >= NODE_PARTICLE_MIN) {
            node->splitLeaf();
            nodeStack.push(node->getLeftChild().get());
            nodeStack.push(node->getRightChild().get());
        } else {
            leaves.push_back(node);
        }
    }
}

std::vector<int> Simulation::getNeighbours(int target, TreeNode& targetNode, bool isCache) {
    if (neighbourCache.contains(target) && isCache) {
        return neighbourCache[target];
    }

    std::vector<int> neighbours;
    std::stack<TreeNode*> nodeStack;
    nodeStack.push(&this->baseNode);

    while (!nodeStack.empty()) {
        TreeNode* nextNode = nodeStack.top();
        nodeStack.pop();

        float distance = distBetweenNodes(*nextNode, targetNode);
        float targetBounds = nextNode->size + targetNode.size + (kernel.getRadius() * targetNode.hmax);

        if (distance * distance < targetBounds * targetBounds) {
            if (nextNode->isLeaf()) {
                for (int candidate : nextNode->getParticleIndices()) {
                    neighbours.push_back(candidate);
                }
            } else {
                nodeStack.push(nextNode->getLeftChild().get());
                nodeStack.push(nextNode->getRightChild().get());
            }
        }
    }

    if (isCache)
        for (int part: targetNode.getParticleIndices())
            neighbourCache.insert(std::make_pair(part, neighbours));

    return neighbours;
}

void Simulation::resetNeighbourCache() {
    neighbourCache.clear();
}

std::vector<int> Simulation::getNeighbours(int part) {
    std::stack<TreeNode*> nodeStack;
    nodeStack.push(&this->baseNode);

    while (!nodeStack.empty()) {
        TreeNode* nextNode = nodeStack.top();
        nodeStack.pop();
        if (!nextNode->isLeaf()) {
            nodeStack.push(nextNode->getLeftChild().get());
            nodeStack.push(nextNode->getRightChild().get());
            continue;
        }

        std::vector<int> indices = nextNode->getParticleIndices();
        if (std::find(indices.begin(), indices.end(), part) != indices.end()) {
            return getNeighbours(part, *nextNode, false);
        }
    }

    std::cout << "Particle " << part << " not found!" << std::endl;
    return {};
}

float Simulation::distBetween(float x1, float x2, float y1, float y2, float z1, float z2) const {
    float dx = std::abs(x1 - x2);
    float dy = std::abs(y1 - y2);
    float dz = std::abs(z1 - z2);

    if (dx > (this->xmax - this->xmin) / 2) {
        dx = (this->xmax - this->xmin) - dx;
    }
    if (dy > (this->ymax - this->ymin) / 2) {
        dy = (this->ymax - this->ymin) - dy;
    }
    if (dz > (this->zmax - this->zmin) / 2) {
        dz = (this->zmax - this->zmin) - dz;
    }

    return std::sqrt(dx * dx + dy * dy + dz * dz);
}


float Simulation::distBetween(int part1, int part2) const {
    float x1 = simData.xyzh[4 * part1], y1 = simData.xyzh[4 * part1 + 1], z1 = simData.xyzh[4 * part1 + 2];
    float x2 = simData.xyzh[4 * part2], y2 = simData.xyzh[4 * part2 + 1], z2 = simData.xyzh[4 * part2 + 2];

    return distBetween(x1, x2, y1, y2, z1, z2);
}

float Simulation::distBetweenNodes(TreeNode& node1, TreeNode& node2) const {
    float x1 = node1.x, y1 = node1.y, z1 = node1.z;
    float x2 = node2.x, y2 = node2.y, z2 = node2.z;

    return distBetween(x1, x2, y1, y2, z1, z2);
}

float Simulation::densityAt(int part, TreeNode& node) {
    std::vector<int> neighbours = getNeighbours(part, node, true);

    float density = 0.0;
    for (int i : neighbours) {
        float dist = distBetween(part, i);
        float tarH = simData.xyzh[4 * part+3];
        density += kernel.valueAt(dist / tarH) * simData.m;
    }

    return density;
}

float Simulation::pressureAt(int part, TreeNode& node) {
    if (!simData.doesContainEnergy()) {
        std::cout << ("No energy found, returning 0!") << std::endl;
        return 0;
    }

    return (GAMMA - 1) * densityAt(part, node) * simData.vxyzu[4 * part + 3];
}

float Simulation::omegaAt(int part, TreeNode& node) {
    float omega = 0;
    float newH = simData.xyzh[part * 4 + 3];
    float grad = -3 * (newH / densityAt(part, node));

    std::vector<int> neighbours = getNeighbours(part, node, true);
    for (int neighbour : neighbours) {
        omega += simData.m * kernel.dWdhAt(distBetween(part, neighbour) / newH);
    }

    return 1 - grad * omega / (newH * newH * newH * newH);
}

float Simulation::densityIterationForParticle(int particle, TreeNode& node) {
    float oldH = std::numeric_limits<float>::max();
    float newH = simData.xyzh[particle * 4 + 3];
    int iterationCount = 0;

    while (std::abs(newH - oldH) / simData.xyzh[particle * 4 + 3] > 10e-4) {
        std::vector<int> neighbours = getNeighbours(particle, node, false);

        float hfact = 1.2;
        float density = simData.m * (hfact / newH) * (hfact / newH) * (hfact / newH);
        float grad = -3 * (newH / density);

        float density_sum = 0;
        float omega = 0;
        for (int neighbour : neighbours) {
            density_sum += simData.m * kernel.valueAt(distBetween(particle, neighbour) / newH) / (newH * newH * newH);
            omega += simData.m * kernel.dWdhAt(distBetween(particle, neighbour) / newH);
        }
        omega = 1 - grad * omega / (newH * newH * newH * newH);

        oldH = newH;
        newH = newH - (density_sum - density) / ((-3 * density * omega) / newH);
        iterationCount++;

        if (newH > 1.4 * oldH) {
            newH = 1.4 * oldH;
        }
        else if (newH < 0.7 * oldH) {
            newH = 0.7 * oldH;
        }

        if (iterationCount > MAX_DENSITY_ITERATIONS) {
            break;
        }
    }
    return newH;
}

Point3f Simulation::velocityDiffBetween(int target, int part) {
    return {simData.vxyzu[4 * target + 0] - simData.vxyzu[4 * part + 0],
                simData.vxyzu[4 * target + 1] - simData.vxyzu[4 * part + 1],
                simData.vxyzu[4 * target + 2] - simData.vxyzu[4 * part + 2]};
}

Point3f Simulation::displacementBetween(int target, int part) {
    return {simData.xyzh[4 * target + 0] - simData.xyzh[4 * part + 0],
                simData.xyzh[4 * target + 1] - simData.xyzh[4 * part + 1],
                simData.xyzh[4 * target + 2] - simData.xyzh[4 * part + 2]};
}

void Simulation::densityIterate() {
    buildTree();
    std::vector<int> neighbours;

    #pragma omp parallel for if(leaves.size() > 1)
    for (std::size_t leafIdx = 0; leafIdx < leaves.size(); leafIdx++) {
        TreeNode* leaf = leaves[leafIdx];
        auto indices = leaf->getParticleIndices();
        for (int i: indices) {
            simData.xyzh[i * 4 + 3] = densityIterationForParticle(i, *leaf);
        }
    }
}

void Simulation::stepSimulation() {
    integrator->step(this->simData, 2e-3);
}

void Simulation::setLimits() {
    xmin = xmax = simData.xyzh[0];
    ymin = ymax = simData.xyzh[1];
    zmin = zmax = simData.xyzh[2];

    for (int i = 0; i < this->getParticleCount(); i++) {
        float x = simData.xyzh[4 * i], y = simData.xyzh[4 * i + 1], z = simData.xyzh[4 * i + 2];

        if (xmin > x) { xmin = x; }
        if (ymin > y) { ymin = y; }
        if (zmin > z) { zmin = z; }
        if (xmax < x) { xmax = x; }
        if (ymax < y) { ymax = y; }
        if (zmax < z) { zmax = z; }
    }
}

void Simulation::setLimits(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax) {
    this->xmin = xmin;
    this->xmax = xmax;
    this->ymin = ymin;
    this->ymax = ymax;
    this->zmin = zmin;
    this->zmax = zmax;
}

int Simulation::getParticleCount() const {
    return this->simData.getParticleCount();
}

SimData& Simulation::getSimData() {
    return this->simData;
}

TreeNode& Simulation::getBaseNode() {
    return this->baseNode;
}

Kernel& Simulation::getKernel() {
    return this->kernel;
}

std::vector<TreeNode*>& Simulation::getLeaves() {
    return this->leaves;
}

float Simulation::qabAt(const int partA, const int partB, TreeNode& nodeA) {
    const Point3f velocityDiff = velocityDiffBetween(partA, partB);
    const Point3f displacement = displacementBetween(partA, partB);
    const Point3f dispNorm = norm(displacement);

    const float losDot = dot(velocityDiff, dispNorm);
    if (losDot >= 0) {
        return 0;
    }

    const float densityA = densityAt(partA, nodeA);
    const float pressureA = pressureAt(partA, nodeA);
    const float soundSpeed = sqrt(GAMMA * pressureA / densityA);
    const float signalSpeed = 1 * soundSpeed + 2 * abs(losDot);

    return -0.5 * pressureA * signalSpeed * losDot;
}

Point3f Simulation::norm(const Point3f p) {
    const float mag = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);

    if (mag <= 0) {
        return Point3f(0, 0, 0);
    }
    return {p.x / mag, p.y / mag, p.z / mag};
}

float Simulation::dot(const Point3f p1, const Point3f p2) {
    return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}
