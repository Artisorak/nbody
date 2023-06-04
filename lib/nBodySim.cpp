#include "nBodySim.hpp"
#include <omp.h>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>

#define TIME(function) \
    do { \
        auto start = std::chrono::high_resolution_clock::now(); \
        function; \
        auto stop = std::chrono::high_resolution_clock::now(); \
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start); \
        std::cout << std::right << "  TIME" << std::setw(8) << duration.count()/1000.0 << " s  for  " << #function << std::endl; \
    } while (0)

nBodySim::nBodySim(std::string datafile) {
    std::string dummy;
    std::ifstream file(datafile);
    file >> dummy >> dummy >> nParticles_;

    for (unsigned i=0; i < nParticles_; ++i) {
        file >> mass_;
    }

    double maxCoord = 0;
    positions_ = new double[3*nParticles_];
    for (unsigned i=0; i < 3; ++i) {
        for (unsigned j=0; j < nParticles_; ++j) {
            file >> positions_[i+3*j];
            if (std::abs(positions_[i+3*j]) > maxCoord) maxCoord = std::abs(positions_[i+3*j]);
        }
    }
    
    velocities_ = new double[3*nParticles_];
    for (unsigned i=0; i < 3; ++i) {
        for (unsigned j=0; j < nParticles_; ++j) {
            file >> velocities_[i+3*j];
        }
    }

    for (unsigned i=0; i < nParticles_; ++i) {
        file >> softening_;
    }
    
    potential_ = new double[nParticles_];
    for (unsigned i=0; i < nParticles_; ++i) {
        file >> potential_[i];
    }

    forces_ = new double[3*nParticles_];
    updateForces();
    
    file.close();

    nodeDepths_ = new unsigned[nParticles_];

    softening_ = 0.01;

    Vector center{0.0, 0.0, 0.0};
    tree_ = new Tree(mass_, positions_, forces_, softening_, nParticles_, maxCoord, center, nodeDepths_);

    maxDepth_ = tree_->getMaxDepth();

    utils_ = Utilities(nParticles_, positions_, velocities_, forces_, mass_, softening_, tree_);
}

nBodySim::~nBodySim() {
    delete tree_;
    delete[] positions_;
    delete[] velocities_;
    delete[] potential_;
    delete[] forces_;
}

void nBodySim::runSimulation(double dt, unsigned nSteps) {
    dt_ = dt;
    // how often to write positions to file
    const unsigned skip = 10;

    // files that store positions and energy
    std::ofstream positionsFile("../out/positions.dat");
    std::ofstream energyFile("../out/energy.dat");

    // write header to positions file
    positionsFile << nParticles_ << " " << skip << " " << dt << std::endl;

    // calculate initial energy
    const double energy0 = utils_.calculateKineticEnergy() + utils_.calculatePotentialEnergy();

    const auto start = std::chrono::high_resolution_clock::now();

    write2file(positions_, positionsFile, 3, 10);

    for (unsigned step=0; step < nSteps; ++step) {

        std::cout << std::right << "Step "    << std::setw(5)  << step+1 << " / " << nSteps << " (" << std::setw(3) << (int)((step+1.0)*100.0/nSteps) << "%)" << std::endl;

        const auto t0 = std::chrono::high_resolution_clock::now();

        // Calculate energy and write to file
        if (step % skip == 0) {
            double ekin, epot;
            TIME(ekin = utils_.calculateKineticEnergy());
            TIME(epot = utils_.calculatePotentialEnergy());
            energyFile << std::right << std::setw(10) << ekin/energy0 << std::setw(10) << epot/energy0 << std::endl;
        }

        // these particles will be updated 
        std::vector<unsigned> particles;
        particles.reserve(nParticles_);
        // mask looks like this: 1, 11, 1, 111, 1, 11, 1, 1111, ...
        // this means that particles at max depth will be updated every step, particles at max depth-1 will be updated every 2nd step, 
        // particles at max depth-2 will be updated every 4th step, etc.
        // if adaptive is false, then mask is 1111...1111, so every particle is updated every step
        unsigned mask = ((step+1) ^ step) | (adaptive_ ? 0 : -1);
        for (unsigned j=0; j < nParticles_; ++j)
            if (mask & (1 << (maxDepth_ - nodeDepths_[j])))
                particles.push_back(j);

        std::cout << std::right << "  Updating " << std::setw(5) << particles.size() << " particles" << std::endl;

        if (particles.size() != 0) {
            // Integrate with leapfrog
            TIME(kick(particles));
            TIME(drift(particles));
            TIME(treeUpdateForces(particles));
            // Verify forces and print average deviation
            if (step % 100 == 0) {
                double dev;
                TIME(dev = utils_.verifyForces());
                std::cout << "  DEV = " << std::setw(12) << dev << std::endl;
            }
            TIME(kick(particles));
            // update entire tree
            TIME(updateTree());
        }

        // Write positions to file
        if (step % 10 == 0) write2file(positions_, positionsFile, 3, 10);
        
        const auto t1 = std::chrono::high_resolution_clock::now();
        const auto stepDuration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);

        const auto now = std::chrono::high_resolution_clock::now();
        const auto totalDuration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start);

        std::cout << "  Duration: "   << std::setw(8)  << stepDuration.count()/1000.0 << " s / Total Duration: " << std::setw(8)  << totalDuration.count()/60000.0 << " minutes" << std::endl << std::endl;
    }
    
    positionsFile.close();
    energyFile.close();
}

void nBodySim::updateForces() {
    #pragma omp parallel for
    for (unsigned i=0; i < 3*nParticles_; ++i) {
        forces_[i] = 0.0;
    }
    // loop over each pair of particles and calculate the force between them
    #pragma omp parallel for // schedule(guided)
    for (unsigned i = 0; i < nParticles_; i++) {
        double px = positions_[3*i];
        double py = positions_[3*i+1];
        double pz = positions_[3*i+2];
        double fx = 0.0;
        double fy = 0.0;
        double fz = 0.0;
        for (unsigned j = 0; j < nParticles_; j++) {
            if (j == i) continue;
            // calculate distance between particles
            double dx = positions_[3*j+0] - px;
            double dy = positions_[3*j+1] - py;
            double dz = positions_[3*j+2] - pz;
            double r2 = dx*dx + dy*dy + dz*dz;
            double r = std::sqrt(r2);

            // calculate forces
            double factor = mass_*mass_ / (r2 + softening_*softening_);
            double fxij = factor * dx / r;
            double fyij = factor * dy / r;
            double fzij = factor * dz / r;
            fx += fxij;
            fy += fyij;
            fz += fzij;
        }
        forces_[3*i]   += fx;
        forces_[3*i+1] += fy;
        forces_[3*i+2] += fz;
    }
}

void nBodySim::treeUpdateForces(const std::vector<unsigned>& particles) {
    #pragma omp parallel for
    for (unsigned i=0; i<particles.size(); ++i) {
        tree_->updateForce(particles[i]);
    }
}

void nBodySim::kick(const std::vector<unsigned>& particles) {
    #pragma omp parallel for
    for (unsigned i=0; i<particles.size(); ++i) {
        unsigned p = particles[i];
        double dt = 0.5*dt_;
        // if adaptive, dt is halved for each level of the tree
        if (adaptive_)
            for (unsigned j=maxDepth_; j>nodeDepths_[p]; --j)
                dt *= 2;
        for (unsigned j=0; j < 3; ++j) {
            velocities_[3*p+j] += forces_[3*p+j] * dt / mass_;
        }
    }
}

void nBodySim::drift(const std::vector<unsigned>& particles) {
    #pragma omp parallel for
    for (unsigned i=0; i<particles.size(); ++i) {
        unsigned p = particles[i];
        double dt = dt_;
        // if adaptive, dt is halved for each level of the tree
        if (adaptive_)
            for (unsigned j=maxDepth_; j>nodeDepths_[p]; --j)
                dt *= 2;
        for (unsigned j=0; j < 3; ++j) {
            positions_[3*p+j] += velocities_[3*p+j] * dt;
        }
    }
}

void nBodySim::updateTree() {
    tree_->update();
}

void nBodySim::write2file(const double* array, std::ofstream& file, unsigned dim, unsigned skip) const {
    // writes array to file, with N rows and dim columns
    for (unsigned i = 0; i < nParticles_; i+=skip) {
        for (unsigned j = 0; j < dim; j++) {
            file << std::right << std::setw(16) << array[dim*i+j] << " ";
        }
        file << std::endl;
    }
}

void nBodySim::saveState2file(unsigned step, std::ofstream& file) const {
    // save step number
    file << "step: " << step << std::endl;
    // save number of particles
    file << nParticles_ << std::endl;
    // save masses
    for (unsigned i=0; i < nParticles_; ++i) {
        file << mass_ << std::endl;
    }
    // save positions
    for (unsigned i=0; i < 3*nParticles_; ++i) {
        file << positions_[i] << std::endl;
    }
    // save velocities
    for (unsigned i=0; i < 3*nParticles_; ++i) {
        file << velocities_[i] << std::endl;
    }
    // save softening
    for (unsigned i=0; i < nParticles_; ++i) {
        file << softening_ << std::endl;
    }
    // save potential
    for (unsigned i=0; i < nParticles_; ++i) {
        file << potential_[i] << std::endl;
    }
}

void nBodySim::saveTree2file(std::string filename) const {
    std::ofstream file(filename);
    tree_->save2file(file);
}