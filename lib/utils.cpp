#include "utils.hpp"

#include <cmath>
#include <omp.h>
#include <fstream>
#include <iomanip>

Utilities::Utilities() {
    nParticles_ = 0;
    positions_ = nullptr;
    velocities_ = nullptr;
    forces_ = nullptr;
    mass_ = 0.0;
    tree_ = nullptr;
}

Utilities::Utilities(unsigned nParticles, double* positions, double* velocities, double* forces, double mass, double softening, Tree* tree) {
    nParticles_ = nParticles;
    positions_ = positions;
    velocities_ = velocities;
    forces_ = forces;
    mass_ = mass;
    softening_ = softening;
    tree_ = tree;
}

double Utilities::verifyForces() const {
    std::ofstream file("../out/forces.txt");
    double avgdev = 0.0;
    double maxdev = 0.0;
    // loop over each pair of particles and calculate the force between them
    #pragma omp parallel for reduction(+:avgdev) reduction(max:maxdev)
    for (unsigned i = 0; i < nParticles_; ++i) {
        double px = positions_[3*i];
        double py = positions_[3*i+1];
        double pz = positions_[3*i+2];
        double fx = 0.0;
        double fy = 0.0;
        double fz = 0.0;
        for (unsigned j = 0; j < nParticles_; ++j) {
            if (j == i) continue;
            // calculate distance between particles
            double dx = positions_[3*j+0] - px;
            double dy = positions_[3*j+1] - py;
            double dz = positions_[3*j+2] - pz;

            double r2 = dx*dx + dy*dy + dz*dz;
            double ir = 1.0/std::sqrt(r2);

            // calculate forces
            double factor = mass_*mass_ / (r2 + softening_*softening_) * ir;
            double fxij = factor * dx;
            double fyij = factor * dy;
            double fzij = factor * dz;
            fx += fxij;
            fy += fyij;
            fz += fzij;
        }
        double f = std::sqrt(fx*fx + fy*fy + fz*fz);
        double af = std::sqrt(forces_[3*i+0]*forces_[3*i+0] + forces_[3*i+1]*forces_[3*i+1] + forces_[3*i+2]*forces_[3*i+2]);
        double dx = fx - forces_[3*i+0];
        double dy = fy - forces_[3*i+1];
        double dz = fz - forces_[3*i+2];
        double df = std::sqrt(dx*dx + dy*dy + dz*dz);
        avgdev += df/f;
        maxdev = std::max(maxdev, df/f);
        #pragma omp critical
        file << std::setw(20) << i << " " << std::setw(20) << f << " " << std::setw(20) << af << " " << std::setw(20) << df/f << std::endl;
    }
    avgdev /= nParticles_;
    file.close();
    return avgdev;
}

double  Utilities::calculateMeanInterparticleDistance() const {
    double meanInterparticleDistance = 0.0;
    // loop over each pair of particles and calculate the distance between them
    #pragma omp parallel for reduction(+:meanInterparticleDistance) // schedule(guided, 64)
    for (unsigned i = 0; i < nParticles_; i++) {
        double px = positions_[3*i];
        double py = positions_[3*i+1];
        double pz = positions_[3*i+2];
        double localMeanInterparticleDistance = 0.0;
        for (unsigned j = i+1; j < nParticles_; j++) {
            // calculate distance between particles
            double dx = positions_[3*j]   - px;
            double dy = positions_[3*j+1] - py;
            double dz = positions_[3*j+2] - pz;
            double r2 = dx*dx + dy*dy + dz*dz;
            double r = std::sqrt(r2);
            localMeanInterparticleDistance += r;
        }
        meanInterparticleDistance += localMeanInterparticleDistance;
    }
    meanInterparticleDistance /= (nParticles_*(nParticles_-1)/2);
    return meanInterparticleDistance;
}

double  Utilities::calculateMeanForceMagnitude() const {
    double forceMagnitude = 0.0;
    #pragma omp parallel for reduction(+:forceMagnitude)
    for (unsigned i = 0; i < nParticles_; i++) {
        double fx = forces_[3*i];
        double fy = forces_[3*i+1];
        double fz = forces_[3*i+2];
        double f2 = fx*fx + fy*fy + fz*fz;
        double f = std::sqrt(f2);
        forceMagnitude += f;
    }
    forceMagnitude /= nParticles_;
    return forceMagnitude;
}

double  Utilities::calculateMeanCenterDistance() const {
    double avgDist = 0.0;
    #pragma omp parallel for reduction(+:avgDist)
    for (unsigned i = 0; i < nParticles_; i++) {
        double x = positions_[3*i];
        double y = positions_[3*i+1];
        double z = positions_[3*i+2];
        double d2 = x*x + y*y + z*z;
        double d = std::sqrt(d2);
        avgDist += d;
    }
    avgDist /= nParticles_;
    return avgDist;
}

double  Utilities::calculateKineticEnergy() const {
    double kineticEnergy = 0.0;
    #pragma omp parallel for reduction(+:kineticEnergy)
    for (unsigned i = 0; i < nParticles_; i++) {
        double vx = velocities_[3*i];
        double vy = velocities_[3*i+1];
        double vz = velocities_[3*i+2];
        double v2 = vx*vx + vy*vy + vz*vz;
        kineticEnergy += 0.5 * mass_ * v2;
    }
    return kineticEnergy;
}

double  Utilities::calculatePotentialEnergy() const {
    double potentialEnergy = 0.0;
    #pragma omp parallel for reduction(+:potentialEnergy)
    for (unsigned i = 0; i < nParticles_; i++) {
        double px = positions_[3*i];
        double py = positions_[3*i+1];
        double pz = positions_[3*i+2];
        double localPotentialEnergy = 0.0;
        for (unsigned j = i+1; j < nParticles_; ++j) {
            // calculate distance between particles
            double dx = positions_[3*j]   - px;
            double dy = positions_[3*j+1] - py;
            double dz = positions_[3*j+2] - pz;
            double r2 = dx*dx + dy*dy + dz*dz;
            double r = std::sqrt(r2);
            potentialEnergy -= 2 * mass_*mass_ / r;
        }
        potentialEnergy += localPotentialEnergy;
    }
    return potentialEnergy;
}

Vector  Utilities::getCenterOfMass() const {
    return tree_->getCOM(); 
}
