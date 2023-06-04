#include <iostream>
#include <fstream>
#include <chrono>
#include <omp.h>
#include "lib/nBodySim.hpp"

int main(int argc, char** argv) {
    std::string datafile;

    const double dt = 1e-7;
    const unsigned nSteps = 14400;

    switch (argc) {
        case 1:
            datafile = "../data/data.ascii";
            break;
        case 2:
            datafile = argv[1];
            break;
        default:
            std::cerr << "Usage: " << argv[0] << " [datafile]" << std::endl;
            return 1;
    }

    std::cout << "datafile: " << datafile << std::endl;
    std::cout << "Timestep = " << dt << std::endl;
    std::cout << "#Steps   = " << nSteps << std::endl;
    std::cout << "simulated time is " << dt*nSteps * 7.76e6 << " Kyr" << std::endl;
    std::cout << "#Threads = " << omp_get_max_threads() << std::endl;
    std::cout << std::endl;

    std::cout << "Initializing simulation..." << std::endl;
    nBodySim sim(datafile);
    std::cout << "Done." << std::endl;
    
    unsigned N = sim.getNParticles();
    std::cout << "N = " << N << std::endl;
    std::cout << std::endl;
    
    // sim.saveTree2file("../out/tree.txt");

    try {
        std::cout << "Running simulation... " << std::endl;

        auto start = std::chrono::high_resolution_clock::now();
        sim.runSimulation(dt, nSteps);
        auto stop = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout << std::endl << "Computation Time: " << duration.count()/60000.0 << " minutes (" << duration.count()/1000.0 << " seconds)" << std::endl;
    } catch (std::exception& e) {
        // save state of simulation
        std::ofstream errfile("../out/errstate.txt");
        sim.saveState2file(nSteps, errfile);
        errfile.close();
        
        std::cerr << e.what() << std::endl;
        return 1;
    }
    
    // save state of simulation
    std::ofstream endfile("../out/endstate.txt");
    sim.saveState2file(nSteps, endfile);
    endfile.close();

    std::cout << "Simulation finished." << std::endl << std::endl;

    return 0;
}