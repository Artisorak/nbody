#include "node.hpp"

#include <vector>
#include <cmath>
#include <iostream>
// #include <omp.h>

Node::Node(double mass, double* positions, unsigned nParticles, double softening, unsigned depth, 
        double size, const Vector& center, const List& localParticles, unsigned nLocalParticles, unsigned* nodeDepths) {
    // this only updates everything that is constant for the node
    mass_ = mass;
    positions_ = positions;
    nParticles_ = nParticles;
    softening_ = softening;
    nodeDepths_ = nodeDepths;

    depth_ = depth;
    size_ = size;
    center_ = center;

    for (unsigned i=0; i<8; ++i) children_[i] = nullptr;

    // update everything that is not constant for the node, and the child nodes
    update(localParticles, nLocalParticles);
}

Node::~Node() {
    for (unsigned i=0; i<8; ++i) {
        delete children_[i];
        children_[i] = nullptr;
    }
}

Node& Node::operator=(const Node& other) {
    positions_ = other.positions_;
    nParticles_ = other.nParticles_;
    depth_ = other.depth_;
    size_ = other.size_;
    totalMass_ = other.totalMass_;
    for (unsigned i=0; i<3; ++i) center_[i] = other.center_[i];
    for (unsigned i=0; i<8; ++i) {
        childParticles_[i] = other.childParticles_[i];
        children_[i] = other.children_[i];
    }
    return *this;
}

void Node::update(const List& localParticles, unsigned nLocalParticles) {
    nLocalParticles_ = nLocalParticles;
    // for each local particle, determine which octant it is in
    // get the mass in this node while we're at it
    for (unsigned i=0; i<8; ++i) {
        childParticles_[i].clear();
    }
    for (unsigned i : localParticles) {
        childParticles_[getOctant(i)].push_back(i);
        nodeDepths_[i] = depth_;
    }

    updateMassAndCOM();
    calculateQ();

    isLeaf_ = true;

    if (depth_ == maxDepth_) {
        // std::cout << "WARNING: Max depth reached" << std::endl;
        for (unsigned i=0; i<8; ++i) {
            delete children_[i];
            children_[i] = nullptr;
        }
        return;
    }

    #pragma omp parallel for num_threads(8) if(depth_ == 0)
    for (unsigned i=0; i<8; ++i) {
        if (children_[i] != nullptr) {
            // if child is not a nullptr and particleDivision[i] contains too few particles, delete the child
            // otherwise update the child
            if (childParticles_[i].size() <= minNParticles_) {
                delete children_[i];
                children_[i] = nullptr;
            } else {
                children_[i]->update(childParticles_[i], childParticles_[i].size());
                isLeaf_ = false;
            }
        } else {
            // if child wouldn't contain enough particles, don't create a new child
            if (childParticles_[i].size() <= minNParticles_) continue;
            
            // if child is a nullptr and particleDivision[i] contains enough particles, create a new child
            double childSize = size_/2;
            Vector childCenter = {0, 0, 0};
            for (unsigned j=0; j<3; ++j) {
                childCenter[j] = center_[j] + (i & (1 << j) ? childSize : -childSize);
            }
            children_[i] = new Node(mass_, positions_, nParticles_, softening_, depth_+1, 
                    childSize, childCenter, childParticles_[i], childParticles_[i].size(), nodeDepths_);
            isLeaf_ = false;
        }
    }
}

Vector Node::calculateForce(unsigned particle) const {
    Vector force = {0, 0, 0};
    // calculate theta
    const double px = positions_[3*particle + 0];
    const double py = positions_[3*particle + 1];
    const double pz = positions_[3*particle + 2];
    const double dx = centerOfMass_[0] - px;
    const double dy = centerOfMass_[1] - py;
    const double dz = centerOfMass_[2] - pz;
    const Vector y = {dx, dy, dz};
    const double ynorm = std::sqrt(dx*dx + dy*dy + dz*dz);
    const double theta = 2*size_/ynorm;
    
    // If theta small enough, calculate force on particle from this node with multipole expansion
    if (theta < theta0_)
        return multipole(y, ynorm);

    // If theta too large: 
    //      For each child that exists, calculate force on particle from that child.
    //      If child doesn't exist, calculate force on particle from particles in that octant.

    // get force for each child that exists and add them together
    for (unsigned c=0; c<8; ++c) {
        // if the child exists, get the force from the child
        if (children_[c] != nullptr) {
            Vector childForce = children_[c]->calculateForce(particle);
            for (unsigned j=0; j<3; ++j) {
                force[j] += childForce[j];
            }
            continue;
        }

        // if the child is a nullptr, calculate the force from the particles in the child (octant) directly
        Vector childForce = {0, 0, 0};
        for (unsigned p : childParticles_[c]) {
            // don't calculate force on particle from itself
            if (p == particle) continue;

            double pdx = positions_[3*p + 0] - px;
            double pdy = positions_[3*p + 1] - py;
            double pdz = positions_[3*p + 2] - pz;

            double r2 = pdx*pdx + pdy*pdy + pdz*pdz;
            double ir = 1.0/std::sqrt(r2);

            double factor = mass_*mass_ / (r2 + softening_*softening_) * ir;

            // add the force from each particle in the child
            childForce[0] += factor * pdx;
            childForce[1] += factor * pdy;
            childForce[2] += factor * pdz;
        }
        for (unsigned j=0; j<3; ++j) {
            force[j] += childForce[j];
        }
    }

    return force;
}

void Node::calculateQ() {
    for (unsigned i=0; i<3; ++i) {
        for (unsigned j=0; j<3; ++j) {
            // calculate the entry in Qij
            double qij = 0;
            // go through every local particle
            for (unsigned c=0; c<8; ++c) {
                for (unsigned k : childParticles_[c]) {
                    // calculate using the formula from page 17 of the pdf for lecture 7
                    qij += 3*(centerOfMass_[i] - positions_[3*k+i])*(centerOfMass_[j] - positions_[3*k+j]);
                    if (i!=j) continue;
                    for (unsigned l=0; l<3; ++l) {
                        double dkl = centerOfMass_[l] - positions_[3*k+l];
                        qij -= dkl*dkl;
                    }

                }
            }
            Q_[3*i+j] = qij;
        }
    }
}

Vector Node::monopole(const Vector& y, double ynorm) const {
    // force = G * mass_ * totalMass_/ynorm^2 * y/ynorm
    Vector result = {0, 0, 0};
    const double ynorm2 = ynorm*ynorm;
    for (unsigned i=0; i<3; ++i) {
        result[i] += mass_*totalMass_/ynorm2 * y[i]/ynorm;
    }
    return result;
}

Vector Node::multipole(const Vector& y, double ynorm) const {
    // force = G * mass_ * (totalMass_/ynorm^2 * y/ynorm - Q * y/ynorm^5 + 5/2 * y^T * Q * y/ynorm^7 * y)
    Vector result = {0, 0, 0};
    const double ynorm2 = ynorm*ynorm;
    const double ynorm5 = ynorm2*ynorm2*ynorm;
    const double ynorm7 = ynorm5*ynorm2;
    
    for (unsigned i=0; i<3; ++i) {
        result[i] += totalMass_/ynorm2 * y[i]/ynorm;
    }
    
    for (unsigned i=0; i<3; ++i) {
        for (unsigned j=0; j<3; ++j) {
            result[i] -= Q_[3*i+j]*y[j] / ynorm5;
        }
    }
    double yQy = 0;
    for (unsigned i=0; i<3; ++i) {
        for (unsigned j=0; j<3; ++j) {
            yQy += y[i]*Q_[3*i+j]*y[j];
        }
    }
    for (unsigned i=0; i<3; ++i) {
        result[i] += 2.5 * yQy * y[i] / ynorm7;
    }
    
    for (unsigned i=0; i<3; ++i) {
        result[i] *= mass_;
    }
    return result;
}

void Node::updateMassAndCOM() {
    totalMass_ = 0;
    centerOfMass_ = {0, 0, 0};
    for (unsigned i=0; i<8; ++i) {
        for (unsigned p : childParticles_[i]) {
            totalMass_ += mass_;
            for (unsigned j=0; j<3; ++j) {
                centerOfMass_[j] += mass_ * positions_[3*p+j];
            }
        }
    }
    for (unsigned j=0; j<3; ++j) {
        centerOfMass_[j] /= totalMass_;
    }
}

Vector Node::getCOM() const {
    return centerOfMass_;
}

unsigned Node::getMaxDepth() const {
    return maxDepth_;
}

unsigned Node::getOctant(unsigned particle) const {
    unsigned result = 0;
    for (unsigned i=0; i<3; ++i) {
        result |= (positions_[3*particle + i] > center_[i]) << i;
    }
    return result;
}

void Node::save2file(std::ofstream& file) const {
    file << center_[0] << " " << center_[1] << " " << center_[2] << " " << size_ << " " << nLocalParticles_ << std::endl;
    for (unsigned i=0; i<8; ++i) {
        if (children_[i] != nullptr) {
            children_[i]->save2file(file);
        }
    }
}