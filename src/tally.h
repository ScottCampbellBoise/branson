/*
 * Test verified with branson's photon.h
 */

#ifndef tally_h
#define tally_h

#include <cmath>

#include "photon.h"

using namespace std;

class Tally {

public:
    
    Tally(double radius ): radius(radius), number_hits(0) {}
    
    ~Tally() {}
    
    // Assumes the photon has started from the origin
    bool hitTally(const Photon phtn){
        const double* cur_pos = phtn.get_position();
        double dist = sqrt( pow(cur_pos[0],2) + pow(cur_pos[1],2) + pow(cur_pos[2],2) );
	if(dist >= radius) {
          number_hits++;
	  return true;
	}
	return false;
    }
    inline int getHits() { return number_hits; }
    inline double getEnergy() { return total_energy; }
    inline void resetHits() { number_hits = 0; }
    inline void resetEnergy() { total_energy = 0; }
    double getExpFlux() { return number_hits / (4 * pi * pow(radius, 2)); }
    double getTheorFlux() { cout << "\t\t ERROR: Theor. Flux Unknown ..." << endl; return -1; }    
    double getFluxRelError() { cout << "\t\tERROR: Rel. Error of flux is unknown ..." << endl; return -1; }

private:
    double radius;
    int number_hits;
    double total_energy;
    const double pi = 3.1415926535897932384626433832795; //!< Pi
};

#endif
