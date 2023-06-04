#ifndef NODE_HPP
#define NODE_HPP

#include <array>
#include <list>
#include <fstream>

using Vector = std::array<double, 3>;
using List = std::list<unsigned>;

class Node {
public:
    // Default constructor creates a leaf node
    Node() = default;
    // Whoever creates the node is responsible for deciding which particles are in which octant and for allocating the localParticles array 
    Node(double mass, double* positions, unsigned nParticles, double softening, unsigned depth, 
        double size, const Vector& center, const List& localParticles, unsigned nLocalParticles, unsigned* nodeDepths);
    // Destuctor
    ~Node();
    // Copy constructor (default)
    Node(const Node&) = default;
    // Copy assignment
    Node& operator=(const Node&);
    
    // Update node: calculate center of mass, get local particles
    void update(const List& particles, unsigned nParticles);
    // Calculate force on a particle recursively, can be called for several particles at the same time by multiple threads 
    Vector calculateForce(unsigned particle) const;
    // Calculate matrix Q for multipole expansion
    void calculateQ();
    // Calculate gravitational force of a node with monopole expansion
    Vector monopole(const Vector& y, double ynorm) const;
    // Calculate gravitational force of a node with multipole expansion 
    Vector multipole(const Vector& y, double ynorm) const;
    // Update center of mass and total mass in the node
    void updateMassAndCOM();
    // get center of mass
    Vector getCOM() const;
    // get maximum depth
    unsigned getMaxDepth() const;
    // Function that takes a particle index in the global particles_ array 
    // and returns the octant it is in relative to the center fo the current node 
    // this function doesn't check if the particle is actually in the node
    unsigned getOctant(unsigned particle) const;
    // save center, size, and number of particles of node to file
    void save2file(std::ofstream& file) const;

private:
    // Node* root_;
    // Node* parent_;
    Node* children_[8];

    double mass_;
    double* positions_;
    unsigned nParticles_;
    double softening_;
    unsigned* nodeDepths_;

    unsigned depth_;
    double size_;
    Vector center_;
    Vector centerOfMass_;

    // Contains the indices of the particles that are in the children of this node for the particles_ array
    List childParticles_[8];
    unsigned nLocalParticles_;

    // Matrix Q for multipole expansion stored in row-major order
    double Q_[9];

    double totalMass_;
    const double theta0_ = 0.4;
    const double minNParticles_ = 8;
    // maximum depth should be below 32, otherwise bad things will happen since we're doing some bit trickery with uints for the adaptive timestepping
    const unsigned maxDepth_ = 19;
    bool isLeaf_;
};

#endif /* NODE_HPP */