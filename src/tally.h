/*
 * Test verified with branson's photon.h
 */

#ifndef tally_h
#define tally_h

#include "photon.h"

//This class creates a Tally object/surface that is a plane oriented along the yz-axis
class Tally {
    
public:
    
    Tally( double x_pos, double y_min, double y_max, double z_min, double z_max ): x_pos(x_pos), y_min(y_min), y_max(y_max), z_min(z_min), z_max(z_max), number_hits(0) {}
    
    ~Tally() {}
    
    // Assumes the photon has started from the origin
    bool hitTally(const Photon* phtn){
        const double* cur_pos = phtn->get_position();
        if( (cur_pos[0] >= x_pos) && (cur_pos[1] >= y_min && cur_pos[1] <= y_max) && (cur_pos[2] >= z_min && cur_pos[2] <= z_max)) {
            number_hits++;
            //total_energy += phtn->get_E();
            return true;
        }
        return false;
    }
    inline int getHits() { return number_hits; }
    inline double getEnergy() { return total_energy; }
    inline void resetHits() { number_hits = 0; }
    inline void resetEnergy() { total_energy = 0; }
    
private:
    double x_pos, y_min, y_max, z_min, z_max;
    int number_hits;
    double total_energy;
};

#endif
