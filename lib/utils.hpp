#ifndef UTILS_HPP
#define UTILS_HPP

#include "tree.hpp"

class Utilities {
public:
    // Constructor
    Utilities();
    Utilities(unsigned nParticles, double* positions, double* velocities, double* forces, double mass, double softening, Tree* tree);
    
    // Loops over all pairs of particles and calculates forces, checks if forces are correct
    double verifyForces() const;
    // Loop over all pairs of particles and calculate mean interparticle distance
    double calculateMeanInterparticleDistance() const;
    // Loop over all force vectors and calculate mean force magnitude
    double calculateMeanForceMagnitude() const;
    // Loop over every particle and calculate the mean distance from the center
    double calculateMeanCenterDistance() const;
    // Loop over all particles and calculate the kinetic energy
    double calculateKineticEnergy() const;
    // Loop over all particles and calculate the potential energy
    double calculatePotentialEnergy() const;
    // Get center of mass
    Vector getCenterOfMass() const;

private:
    unsigned nParticles_;
    double* positions_;
    double* velocities_;
    double* forces_;
    double mass_;
    double softening_;
    Tree* tree_;
};

#endif /* UTILS_HPP */