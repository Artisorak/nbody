#include "tree.hpp"
#include "node.hpp"

#include <vector>

Tree::Tree(double mass, double* positions, double* forces, double softening, unsigned nParticles, double size, Vector center, unsigned* nodeDepths) {
    mass_ = mass;
    positions_ = positions;
    forces_ = forces;
    nParticles_ = nParticles;
    softening_ = softening;
    size_ = size;
    nodeDepths_ = nodeDepths;

    center_ = Vector{0.0, 0.0, 0.0};
    for (unsigned i=0; i < 3; ++i) center_[i] = center[i];

    allParticles_ = {};
    for (unsigned i=0; i < nParticles_; ++i) allParticles_.push_back(i);

    root_ = new Node(mass_, positions_, nParticles_, softening_, 0, size_, center_, allParticles_, nParticles_, nodeDepths_);
}

Tree::~Tree() {
    delete root_;
}

void Tree::updateForce(unsigned particle) {
    Vector force = root_->calculateForce(particle);
    for (unsigned i=0; i < 3; ++i) forces_[3*particle + i] = force[i];
}

void Tree::update() {
    root_->update(allParticles_, nParticles_);
}

Vector Tree::getCOM() const {
    return root_->getCOM();
}

unsigned Tree::getMaxDepth() const {
    return root_->getMaxDepth();
}

void Tree::save2file(std::ofstream& file) const {
    root_->save2file(file);
}