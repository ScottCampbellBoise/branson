/*
 * @file: tally.h
 * @author: Scott Campbell
 * @version: 14 June 2019
 * Test verified with branson's photon.h
 */

#ifndef tally_h
#define tally_h

#include <cmath>
#include "constants.h"
#include "photon.h"

using namespace std;
using Constants::pi;

// This class creates a 'bar' tally
// Test verified that the tally will accurately count
//     particles only if they travel through the tally
//     within their last time/space step. Also verified 
//     the motion filter settings.

class Tally {

public:
    
    Tally(double x_1, double y_1, double x_2, double y_2 ) : x_1(x_1), 
	  y_1(y_1), x_2(x_2), y_2(y_2), n_hits(0), filter_state(NO_FILTER) {
	// Define the line equation to connect the two tally points
        tally_m = (y_1 - y_2) / (x_1 - x_2);
        tally_b = (-tally_m * x_1) + y_1;
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
    inline double get_x1() { return x_1; }
    inline double get_y1() { return y_1; }
    inline double get_x2() { return x_2; }
    inline double get_y2() { return y_2; }
    inline double get_tally_slope() { return tally_m; }
    inline double get_tally_intcpt() { return tally_b; }
    inline int get_motion_filter() { return filter_state; }

    inline set_tally_pos(double x1, double y1, double x2, double y2) {
	x_1 = x1;
	y_1 = y1;
	x_2 = x2;
	y_2 = x2;
	tally_m = (y_1 - y_2) / (x_1 - x_1);
	tally_b = (-tally_m * x_1) + y_1;
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
	
	Point p2 = {x_1, y_1};
	Point q2 = {x_2, y_2};

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

    double x_1, y_1, x_2, y_2;
    double tally_m, tally_b; // Define the line between the tally y = m*x + b
    int n_hits;
    double total_energy;
    int filter_state;

};

#endif
