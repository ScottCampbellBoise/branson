
/*
 * @file: plane_tally.h
 * @author: Scott Campbell
 * @version: 14 June 2019
 * Test verified with branson's photon.h
 */

#ifndef plane_tally_h
#define plane_tally_h

#include <cmath>
#include "constants.h"  // Use for standard math constants
#include "photon.h"     // Use for particle information: pos, angle, weight
#include "mesh.h"

using namespace std;
using Constants::pi;

class Plane_Tally {
    
public:
    
    Plane_Tally(double* v1, double* v2, double* v3, double* v4, Mesh& mesh) :
    v1(v1), v2(v2), v3(v3), v4(v4), mesh(mesh),
    total_regular_E(0), total_response_E(0), n_response_hits(0),
    n_regular_hits(0) {}
    
    ~Plane_Tally() {}
    
    inline double get_regular_E() { return total_regular_E; }
    inline double get_response_E() { return total_response_E; }
    inline void reset_regular_E() { total_regular_E = 0; }
    inline void reset_response_E() { total_response_E = 0; }
    inline void add_response_weight(double weight) { total_response_E += weight; }
    inline void add_regular_weight(double weight) { total_regular_E += weight; }
    
    inline uint32_t get_regular_hits() { return n_regular_hits; }
    inline uint32_t get_response_hits() { return n_response_hits; }
    inline void reset_regular_hits() { n_regular_hits = 0; }
    inline void reset_response_hits() { n_response_hits = 0; }
    inline void add_regular_hit() { n_regular_hits++; }
    inline void add_response_hit() { n_response_hits++; }
   
    inline double get_x1() { return v1[0]; }
 
    // Find the distance of the photon to tally surface
    double get_dist_to_tally(Photon& phtn) {
        const double* pos = phtn.get_position();
        const double* ang = phtn.get_angle();

        double x_min = min(v1[0], min(v2[0], min(v3[0], v4[0])));
        double x_max = max(v1[0], max(v2[0], max(v3[0], v4[0])));
        double y_min = min(v1[1], min(v2[1], min(v3[1], v4[1])));
        double y_max = max(v1[1], max(v2[1], max(v3[1], v4[1])));
        double z_min = min(v1[2], min(v2[2], min(v3[2], v4[2])));
        double z_max = max(v1[2], max(v2[2], max(v3[2], v4[2])));   

        double vec_12[3] = {v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]};
        double vec_13[3] = {v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2]};
        
	// Equation of the plane: Ax + By + Cz = D
	    
        double A = vec_12[1]*vec_13[2] - vec_12[2]*vec_13[1];
        double B = vec_12[2]*vec_13[0] - vec_12[0]*vec_13[2];
        double C = vec_12[0]*vec_13[1] - vec_12[1]*vec_13[0];
        double D = A*v1[0] + B*v1[1] + C*v1[2];
        
        if(ang[0] < 0 || pos[0] > x_max || (A*ang[0] + B*ang[1] + C*ang[2] == 0)) {return 1e100;}
        
        double d = ((v1[0]-pos[0])*A + (v1[1]-pos[1])*B + (v1[2]-pos[2])*C) /
                   (ang[0]*A + ang[1]*B + ang[2]*C);
        
        double x_int = d*ang[0] + pos[0];
        double y_int = d*ang[1] + pos[1];
        double z_int = d*ang[2] + pos[2];
        
        double dist = sqrt(pow(x_int-pos[0],2) + pow(y_int-pos[1],2) + pow(z_int-pos[2],2));

	if(y_int <= y_max && y_int >= y_min 
	   && z_int <= z_max && z_int >= z_min) {return dist;} 
        
        return 1e100;
    }
    
    void get_point_on_plane(double*& pos) {
        double x_min = min(v1[0], min(v2[0], min(v3[0], v4[0])));
        double x_max = max(v1[0], max(v2[0], max(v3[0], v4[0])));
        double y_min = min(v1[1], min(v2[1], min(v3[1], v4[1])));
        double y_max = max(v1[1], max(v2[1], max(v3[1], v4[1])));
        double z_min = min(v1[2], min(v2[2], min(v3[2], v4[2])));
        double z_max = max(v1[2], max(v2[2], max(v3[2], v4[2])));
        
        pos[0] = x_min + ((double)rand() / ((double)RAND_MAX)) * (x_max - x_min);
        pos[1] = y_min + ((double)rand() / ((double)RAND_MAX)) * (y_max - y_min);
        pos[2] = z_min + ((double)rand() / ((double)RAND_MAX)) * (z_max - z_min);
    }
    
private:
    // ----------------------------------------------------------------------
    // Private Variables/Constants
    // ----------------------------------------------------------------------
    Mesh& mesh;
    
    double* v1; // corner 1
    double* v2; // corner 2
    double* v3; // corner 3
    double* v4; // corner 4
    
    int n_response_hits;
    int n_regular_hits;
    
    double total_regular_E;
    double total_response_E;
    
};

#endif
