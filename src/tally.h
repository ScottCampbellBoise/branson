/*
 * Test verified with branson's photon.h
 */

#ifndef tally_h
#define tally_h

#include <math.h>

#include "photon.h"

//This class creates a Tally object/surface that is a sphere centered at the origin
class Tally {
    
public:
    
    Tally(double radius ): radius(radius), number_hits(0) {}
    
    ~Tally() {}
    
    // Assumes the photon has started from the origin
    bool hitTally(const Photon* phtn){
        const double* cur_pos = phtn->get_position();
        double dist = sqrt( pow(cur_pos[0],2) + pow(cur_pos[1],2) + pow(cur_pos[2],2) );
	if(dist >= radius) {
          numberHits++;
	  return true;
	}
	return false;
    }
    inline int getHits() { return number_hits; }
    inline double getEnergy() { return total_energy; }
    inline void resetHits() { number_hits = 0; }
    inline void resetEnergy() { total_energy = 0; }
    
private:
    double radius;
    int number_hits;
    double total_energy;
};

#endif
