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
#include "input.h"	// Contains information about the mesh dimesions

using namespace std;
using Constants::pi;

// This class creates a 'bar' tally
// Test verified that the tally will accurately count
//     particles only if they travel through the tally
//     within their last time/space step. Also verified 
//     the motion filter settings.
// Test verified to properly get the x and y dimensions
//     of the mesh object through an Input object

class Tally {

public:
    
    // NOTE: The Mesh class jhas been updated to have ref. to Input. Could
    // 	alternatively pass a mesh as an object, and extract the input
    // 	out of the mesh object.
    
    Tally(double x1, double y1, double x2, double y2, const Input& input) : x1(x1), 
	  y1(y1), x2(x2), y2(y2), n_hits(0), filter_state(NO_FILTER) {
	// Define the line equation to connect the two tally points
        tally_m = (y1 - y2) / (x1 - x2);
        tally_b = (-tally_m * x1) + y1;
   	
	// Set up the mesh boundary information
	mesh_start_x = input.get_x_start(0);
	mesh_start_y = input.get_y_start(0);
	mesh_start_z = input.get_z_start(0);
	mesh_end_x = input.get_x_end(input.get_n_x_divisions() - 1);
	mesh_end_y = input.get_y_end(input.get_n_y_divisions() - 1);
	mesh_end_z = input.get_z_end(input.get_n_z_divisions() - 1);     

	// Check that the tally line is fully inside/on the mesh surface
	if(x1 >= mesh_start_x && x1 <= mesh_end_x &&
	   x2 >= mesh_start_x && x2 <= mesh_end_x &&
	   y1 >= mesh_start_y && y1 <= mesh_end_y &&
	   y2 >= mesh_start_y && y2 <= mesh_end_y) {
	    std::cout << "\n\nERROR: The tally line MUST be inside or on" <<
		" the mesh surface!\n\n";
	    return;
	}	
    }
    
    ~Tally() {}
    
    bool hit_tally(const Photon phtn){
	if(passed_tally(phtn)) {
	    n_hits++;
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
    inline double get_x2() { return x2; }
    inline double get_y2() { return y2; }
    inline double get_tally_slope() { return tally_m; }
    inline double get_tally_intercept() { return tally_b; }
    inline int get_motion_filter() { return filter_state; }

    inline double get_mesh_start_x() { return mesh_start_x; }
    inline double get_mesh_end_x() { return mesh_end_x; }
    inline double get_mesh_start_y() { return mesh_start_y; }
    inline double get_mesh_end_y() { return mesh_end_y; }
    inline double get_mesh_start_z() { return mesh_start_z; }
    inline double get_mesh_end_z() { return mesh_end_z; }
    
    inline void set_tally_pos(double x1, double y1, double x2, double y2) {
	x1 = x1;
	y1 = y1;
	x2 = x2;
	y2 = x2;
	tally_m = (y1 - y2) / (x1 - x1);
	tally_b = (-tally_m * x1) + y1;
    }
    bool set_motion_filter(int new_filter) {
	if(new_filter < -1 || new_filter > 1) {
	    return false;
	}
	filter_state = new_filter;
	return true;
    }

    // Variables regarding the motion filter settings
    const static int NO_FILTER = 0;
    const static int POSITIVE_ONLY_FILTER = 1; // Only particles moving in (+) dir
    const static int NEGATIVE_ONLY_FILTER = -1; // " " in (-) dir

private:
    
    // ----------------------------------------------------------------------
    // Private helper methods/structs
    // ----------------------------------------------------------------------
    
    struct Point {
	double x;
	double y;
    };

    // determine if the photon has passed through the tally in the last time step
    bool passed_tally(Photon phtn) { 
	const double* phtn_pos = phtn.get_position();
	const double* phtn_prev_pos = phtn.get_prev_position();	

	// if both the prev point and the cur point are above/below tally,
	// then there is no way that it can be passing through the tally
	if(((phtn_pos[1] > tally_m * phtn_pos[0] + tally_b) &&
	    (phtn_prev_pos[1] > tally_m * phtn_prev_pos[0] + tally_b)) ||
	   ((phtn_pos[1] < tally_m * phtn_pos[0] + tally_b) &&
	    (phtn_prev_pos[1] < tally_m * phtn_prev_pos[0]))) {
	    return false;
	}

	//Check if there is a motion filter on
	if(filter_state == POSITIVE_ONLY_FILTER) {
	    if(phtn_pos[1] < tally_m * phtn_pos[0] + tally_b) {
		return false;
	    }
	} else if(filter_state == NEGATIVE_ONLY_FILTER) {
	    if(phtn_pos[1] > tally_m * phtn_pos[0] + tally_b) {
		return false;
	    }
	}

        Point p1 = {phtn_pos[0], phtn_pos[1]};
	Point q1 = {phtn_prev_pos[0], phtn_prev_pos[1]};
	
	Point p2 = {x1, y1};
	Point q2 = {x2, y2};

	// get the 4 orientations possible
	double o1 = find_orient(p1, q1, p2); 
	double o2 = find_orient(p1, q1, q2);
	double o3 = find_orient(p2, q2, p1);
	double o4 = find_orient(p2, q2, q1);

	// Chek if the photons line seg intersects the tally line
	return (o1 != o2 && o3 != o4) ||
	     (o1 == 0 && on_seg(p1, p2, q1)) ||
	     (o2 == 0 && on_seg(p1, q2, q1)) ||
	     (o3 == 0 && on_seg(p2, p1, q2)) ||
	     (o4 == 0 && on_seg(p2, q1, q2));
    }

    //find the orientation of the triplet (p, q, r)
    // 0 -> p, q & r are colinear
    // 1 -> Clockwise
    // 2 -> Counterclockwise
    int find_orient(Point p, Point q, Point r) {
	double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
	if(val == 0) { return 0; }
	return (val > 0)? 1 : 2;
    }

    // For colinear points p, q, & r check if q lies on seg 'pr'
    inline bool on_seg(Point p, Point q, Point r) {
	return (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) && 
		q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y));
    }

    // ----------------------------------------------------------------------
    // Private Variables/Constants
    // ----------------------------------------------------------------------
    
    double x1, y1, x2, y2;
    double tally_m, tally_b; // Define the line between the tally y = m*x + b
    int n_hits;
    double total_energy;
    int filter_state;
    
    // Variables to store information about the mesh dimensions
    double mesh_start_x, mesh_start_y, mesh_start_z;
    double mesh_end_x, mesh_end_y, mesh_end_z;

};

#endif
