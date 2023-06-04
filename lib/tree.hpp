#ifndef TREE_HPP
#define TREE_HPP

#include "node.hpp"

#include <vector>
#include <fstream>

class Tree {
public:
    // Constructor
    Tree() = delete;
    Tree(double mass, double* positions, double* forces_, double softening, unsigned nParticles, double size, Vector center, unsigned* nodeDepth);
    // Destructor
    ~Tree();
    // Update force on a particle
    void updateForce(unsigned particle);
    // Update tree
    void update();
    // get center of mass
    Vector getCOM() const;
    // get maximum depth
    unsigned getMaxDepth() const;
    // save tree to file
    void save2file(std::ofstream& file) const;

private:
    Node* root_;

    double mass_;
    double* positions_;
    double* forces_;
    unsigned nParticles_;
    double softening_;
    unsigned* nodeDepths_;

    // half of the size of the domain
    double size_;
    Vector center_;

    // allParticles contains all indices, so [0, 1, 2, ..., nParticles_-1]. It should never be changed. 
    List allParticles_;
};

#endif /* TREE_HPP */