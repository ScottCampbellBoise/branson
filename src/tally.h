/*
 * @file: tally.h
 * @author: Scott Campbell
 * @version: 14 June 2019
 * Test verified with branson's photon.h
 */

#ifndef tally_h
#define tally_h

#include <cmath>
#include "constants.h"	// Use for standard math constants
#include "photon.h"	// Use for particle information: pos, angle, weight
#include "mesh.h"	

using namespace std;
using Constants::pi;

class Tally {

public:
    
    Tally(double radius, double x1, double y1, double z1, Mesh& mesh) : 
	  radius(radius), x1(x1), y1(y1), z1(z1), mesh(mesh), 
	  n_hits(0), total_energy(0), filter_state(NO_FILTER) {}
    
    ~Tally() {}
    
    bool hit_tally(Photon phtn) {
        if(passed_spherical_tally(phtn)) {
	    n_hits++;
	    total_energy += phtn.get_E();
	    return true;
	}
	return false;
    }

    inline int get_hits() { return n_hits; }
    inline double get_E() { return total_energy; }
    inline void reset_hits() { n_hits = 0; }
    inline void reset_E() { total_energy = 0; }
    
    inline double get_x1() { return x1; }
    inline double get_y1() { return y1; }
    inline double get_z1() { return z1; }
    inline double get_radius() { return radius; }
    inline int get_motion_filter() { return filter_state; }
    
    inline void set_tally_pos(double new_r, double new_x1, double new_y1, double new_z1) {
	radius = new_r;
	x1 = new_x1;
	y1 = new_y1;
	z1 = new_z1;
    }
    
    bool set_motion_filter(int new_filter) {
	if(new_filter < -1 || new_filter > 1)
	    return false;
	filter_state = new_filter;
	return true;
    }

    inline void add_weight(double weight) { total_energy += weight; }

    void print_tally_info(ostream& outfile, double time_step) {
	outfile << time_step << " , " << total_energy << endl;
    }

    bool is_inside_tally(Photon& phtn) {
	const double* pos = phtn.get_position();
	double dist = sqrt(pow(pos[0]-x1,2) + pow(pos[1]-y1,2) + pow(pos[2]-z1,2));
	return dist <= radius;
    }

    // Find the distance of the photon to tally surface
    double get_dist_to_tally(Photon& phtn) {
	const double* pos = phtn.get_position();
	const double* ang = phtn.get_angle();
	const double* prev_pos = phtn.get_prev_position();

	// Check that the photon intersects the sphere
	// equation: a*t^2 + b*t + c = 0
	double a = pow(prev_pos[0] - pos[0], 2) + 
		   pow(prev_pos[1] - pos[1], 2) +
		   pow(prev_pos[2] - pos[2], 2);
	double b = -2*((prev_pos[0] - pos[0]) * (x1 - pos[0]) + 
		       (prev_pos[1] - pos[1]) * (y1 - pos[1]) +
		       (prev_pos[2] - pos[2]) * (z1 - pos[2]));
	double c = pow(x1 - pos[0], 2) + pow(y1 - pos[1], 2) + 
		   pow(z1 - pos[2], 2) - pow(radius, 2);
	
	if((pow(b, 2) - 4*a*c) < 0) { return 1e9; } // No intersection

	// Solve the quadratic equation for t
	double t_a = (-b + sqrt(pow(b, 2) - 4*a*c)) / (2*a);
	double t_b = (-b - sqrt(pow(b, 2) - 4*a*c)) / (2*a);	

	double pt_a[3] = { (prev_pos[0] - pos[0])*t_a + pos[0], 
			   (prev_pos[1] - pos[1])*t_a + pos[1],
			   (prev_pos[2] - pos[2])*t_a + pos[2] };
	double pt_b[3] = { (prev_pos[0] - pos[0])*t_b + pos[0], 
			   (prev_pos[1] - pos[1])*t_b + pos[1],
			   (prev_pos[2] - pos[2])*t_b + pos[2] };
	
	double dist_a = sqrt(pow(pos[0]-pt_a[0],2) + pow(pos[1]-pt_a[1],2) + pow(pos[2]-pt_a[2],2));
	double dist_b = sqrt(pow(pos[0]-pt_b[0],2) + pow(pos[1]-pt_b[1],2) + pow(pos[2]-pt_b[2],2));

	return min(dist_a, dist_b);
    }

    // Variables regarding the motion filter settings
    const static int NO_FILTER = 0;
    const static int POSITIVE_ONLY_FILTER = 1; // Only particles moving in (+) dir / into sphere
    const static int NEGATIVE_ONLY_FILTER = -1; // " " in (-) dir / out of sphere

private:
    
    // ----------------------------------------------------------------------
    // Private helper methods
    // ----------------------------------------------------------------------
   
    bool passed_spherical_tally(Photon& phtn) {
	const double* phtn_pos = phtn.get_position();
	const double* phtn_prev_pos = phtn.get_prev_position();
	
	//Calculate the distance the phtn is from the tally center
	double cur_dist = sqrt(pow(phtn_pos[0] - x1,2) + pow(phtn_pos[1] - y1,2) + pow(phtn_pos[2] - z1,2));
	double prev_dist = sqrt(pow(phtn_prev_pos[0] - x1,2) + pow(phtn_prev_pos[1] - y1,2) + pow(phtn_prev_pos[2] - z1,2));
	
	if(filter_state == POSITIVE_ONLY_FILTER) {
	    // Count only particles moving into the sphere
	    // total_energy += 
	    return (cur_dist <= radius && prev_dist > radius);
	} else if(filter_state == NEGATIVE_ONLY_FILTER) {
	    // Count only the particles moving out of the sphere
	    // total_energy -=
	    return (cur_dist > radius && prev_dist < radius);
	} else {
	    return (cur_dist <= radius && prev_dist > radius) ||
		   (cur_dist > radius && prev_dist < radius);
	}
    }

    // ----------------------------------------------------------------------
    // Private Variables/Constants
    // ----------------------------------------------------------------------
    Mesh& mesh;   
    
    double radius, x1, y1, z1;
   
    int n_hits;
    double total_energy;
    int filter_state;
};

#endif
