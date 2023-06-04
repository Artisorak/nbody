#ifndef NBODYSIM_HPP
#define NBODYSIM_HPP

#include "tree.hpp"
#include "utils.hpp"
#include <string>

class nBodySim {
public:
    // Constructor reads in data from datafile
    nBodySim(std::string datafile);
    // Destructor deletes arrays 
    ~nBodySim();

    // Run simulation for nSteps time steps of size dt, integration with leapfrog algorithm
    void runSimulation(double dt, unsigned nSteps);

    // Loops over all pairs of particles and calculates forces 
    void updateForces();
    // Use treecode and multipole expansion to calculate forces
    void treeUpdateForces(const std::vector<unsigned>& particles);

    // Leapfrog: update velocities
    void kick(const std::vector<unsigned>& particles);
    // Leapfrog: update positions
    void drift(const std::vector<unsigned>& particles);
    // Update tree
    void updateTree();
    
    // Writes data to file, every row is a particle
    void write2file(const double* array, std::ofstream& file, unsigned dim, unsigned skip) const;
    // Write current state of simulation to file
    void saveState2file(unsigned step, std::ofstream& file) const;
    // save location and size of all nodes in tree
    void saveTree2file(std::string filename) const;

    // Functions to get private variables
    double getMass() { return mass_; }
    double* getPositions() { return positions_; }
    double* getForces() { return forces_; }
    unsigned getNParticles() const { return nParticles_; }

private:
    Utilities utils_;

    Tree* tree_;
    unsigned maxDepth_;
    // whether to use adaptive timestep
    const bool adaptive_ = true;
    double dt_;

    unsigned nParticles_;
    // 3d vectors are stored as a array of length 3*nParticles_ in the format x1, y1, z1, x2, y2, z2, ...
    double mass_;
    double* positions_;
    double* velocities_;
    double softening_;
    double* potential_;
    double* forces_;
    unsigned* nodeDepths_;
};

#endif /* NBODYSIM_HPP */