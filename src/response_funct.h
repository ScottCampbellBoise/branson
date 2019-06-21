/*
 * @file:    response_funct.h
 * @author:  Scott Campbell
 * @version: 18 June 2019
 * @desc:    This class stores/holds the response value for a each
 *           cell in a grid system.
 */

#ifndef response_funct_h
#define response_funct_h

#include <cmath>
#include <algorithm>

#include "tally.h"
#include "mesh.h"

using namespace std;

class Response_Funct {

// Testing Results:
// 	The method to check if a line (y = mx + b) passes
// 	through a cell (box) is test verified to work as 
// 	expected.
// To Do:
// 	for the simple response function table function,
// 	define lines outward from the tally. NOTE: for
// 	simplcity, the tally MUST be placed on a cell 
// 	boundary. V1 will assume its vertical on the 
// 	rightmost wall. ALSO: each cell must have a
// 	uniform sigma_a value
//
//	** FINISH THE CONTRIBUTION METHOD **
//
//	Make a way to find the sigma_a value of a cell ...
public:

    struct Point {
	double x;
	double y;
    };

    struct Ray{
	double m; // Slope of the ray
	double b; // Y-Intercept of the ray
    };
	
    // Constructor: store the info about the mesh and tally
    Response_Funct(Tally& tally, Mesh& mesh) : tally(tally),
					       mesh(mesh) {
	mesh_start_x = mesh.get_start_x();
	mesh_end_x = mesh.get_end_x();
	mesh_start_y = mesh.get_start_y();
	mesh_end_y = mesh.get_end_y();
	mesh_dx = mesh.get_dx();
	mesh_dy = mesh.get_dy();
	mesh_dz = mesh.get_dz();

	tally_x1 = tally.get_x1();
	tally_y1 = tally.get_y1();
	tally_x2 = tally.get_x2();
	tally_y2 = tally.get_y2();
    }   

    ~Response_Funct() {}

    // This uses each cell as a region to accumulate over
    // Assumes the tally is on the far right wall ...
    double** get_simple_response(int n_points) {
	double** resp_table = get_blank_table(mesh.get_n_x_cells(),
					      mesh.get_n_y_cells());
	
	Ray ray;
	double slopes[3] = {0, -1, 1};
	
	for(double point_scale = 0; point_scale <= 1.0; point_scale += (1/n_points)) {
	    for(int m = 0; m < 3; m++) {
		// Generate a Ray object that intersects the tally at each point
		double x_coor = (tally.get_x1() + tally.get_x2()) * point_scale;
		double y_coor = (tally.get_y1() + tally.get_y2()) * point_scale;
		
		// Normalize to the standard grid
		double adj_y_int = slopes[m] * x_coor - y_coor;
		ray = {slopes[m] , adj_y_int};
		
		// Check all the cells to see if the ray goes through
		adjust_cells(resp_table, ray);	
	    }

	}
	return resp_table;
    }

    // Iterate through all cells and adjust their table value if the ray goes through
    void adjust_cells(double**& table, Ray ray) {
 	// Do a simple tracking of the photon ...
    }

    double** get_blank_table(int n_rows, int n_cols) {
	double** table = new double*[n_rows];
	for(int row = 0; row < n_rows; row++) {
	    table[row] = new double[n_cols];
	    for(int col = 0; col < n_cols; col++) {
	   	table[row][col] = 0;
	    }
	}
	return table;
    }
    
    double get_contribution(double dist, double sigma_a) {
	// FINISH THIS SECTION
	return 0;
    }

    Point make_point(double x, double y) {
   	Point p = {x, y};
	return p;
    }

    Ray make_ray(double m, double b) {
	Ray r = {m, b};
	return r;
    }

    // Check if a Ray intersects the cell defined by its
    // lower-left point and its upper right point, and if 
    // so, return the length of intersection. else -1
    double intersect_cell(Point l_left, Point u_right, Ray ray) {
	// Check if the vector intersects the cell
	Point l_right = {u_right.x, l_left.y};
	Point u_left = {l_left.x, u_right.y};
	double l_left_resid = l_left.y - ray.m * l_left.x - ray.b;
	double l_right_resid = l_right.y - ray.m * l_right.x - ray.b;
	double u_left_resid = u_left.y - ray.m * u_left.x - ray.b;
	double u_right_resid = u_right.y - ray.m * u_right.x - ray.b;

	double min_resid = min(l_left_resid, min(l_right_resid, min(u_left_resid, u_right_resid)));
	double max_resid = max(l_left_resid, max(l_right_resid, max(u_left_resid, u_right_resid)));

	if(min_resid <= 0 && max_resid >= 0) { // ray intersects the cell
	    // Calculate the distance the ray goes through the cell
	    
	    // If the ray is on the edge or goes horiz. through
	    if(ray.m == 0) { 
		return (l_right.x - l_left.x);
	    }
	    // if the ray is on the edge or eff. vertical through 
	    else if(ray.m > 1e6) {
		return (u_left.y - l_left.y);
	    }
	    double resid_array[4] = {l_left_resid,l_right_resid,u_right_resid,u_left_resid};
	    sort(resid_array, resid_array + 4);
            
	    //Case 1: three of the residuals are either positive or negative
	    if(resid_array[0] < 0 && resid_array[1] > 0) {
		if(resid_array[0] == l_left_resid) {
		    return sqrt(pow(l_left.x - (l_left.y-ray.b)/ray.m,2) +
			   pow((ray.m*l_left.x + ray.b) - l_left.y,2));
		} else if(resid_array[0] == l_right_resid) {
		    return sqrt(pow(l_right.x - (l_right.y-ray.b)/ray.m,2) +
			   pow((ray.m*l_right.x + ray.b) - l_right.y,2));
		} else if(resid_array[0] == u_left_resid) {
		    return sqrt(pow(u_left.x - (u_left.y-ray.b)/ray.m,2) +
			   pow((ray.m*u_left.x + ray.b) - u_left.y,2));
		} else  {
		    return sqrt(pow(u_right.x - (u_right.y-ray.b)/ray.m,2) +
			   pow((ray.m*u_right.x + ray.b) - u_right.y,2));
		}
	    } else if (resid_array[2] < 0 && resid_array[3] > 0) {
		if(resid_array[3] == l_left_resid) {
		    return sqrt(pow(l_left.x - (l_left.y-ray.b)/ray.m,2) +
			   pow((ray.m*l_left.x + ray.b) - l_left.y,2));
		} else if(resid_array[3] == l_right_resid) {
		    return sqrt(pow(l_right.x - (l_right.y-ray.b)/ray.m,2) +
			   pow((ray.m*l_right.x + ray.b) - l_right.y,2));
		} else if(resid_array[3] == u_left_resid) {
		    return sqrt(pow(u_left.x - (u_left.y-ray.b)/ray.m,2) +
			   pow((ray.m*u_left.x + ray.b) - u_left.y,2));
		} else  {
		    return sqrt(pow(u_right.x - (u_right.y-ray.b)/ray.m,2) +
			   pow((ray.m*u_right.x + ray.b) - u_right.y,2));
		}
	    }	
	    // Case 2: two points have the same parity
	    else {
		if(l_left_resid < 0 && l_right_resid < 0 ) {
		    return sqrt(pow((l_left.x - l_right.x),2) +
			   pow((ray.m*l_left.x - ray.b)-(ray.m*l_right.x) ,2));
		} else if(u_left_resid < 0 && u_right_resid < 0) {
		    return sqrt(pow((u_left.x - u_right.x),2) +
			   pow((ray.m*u_left.x - ray.b)-(ray.m*u_right.x) ,2));
		} else if(l_left_resid < 0 && u_right_resid < 0) {	
		    return sqrt(pow((l_left.x - u_right.x),2) +
			   pow((ray.m*l_left.x - ray.b)-(ray.m*u_right.x) ,2));
		} else if(l_left_resid < 0 && u_left_resid < 0) {
		    return sqrt(pow((l_left.x - u_left.x),2) +
			   pow((ray.m*l_left.x - ray.b)-(ray.m*u_left.x) ,2));
		} else if(l_right_resid < 0 && u_right_resid < 0) {
		    return sqrt(pow((l_right.x - u_right.x),2) +
			   pow((ray.m*l_right.x - ray.b)-(ray.m*u_right.x) ,2));
		} else { // lower right and upper left
		    return sqrt(pow((l_right.x - u_left.x),2) +
			   pow((ray.m*l_right.x - ray.b)-(ray.m*u_left.x) ,2));
		}
	    }
	}
	return -1; // Indicate that the ray does not intersect the cell
    }

private:

    Mesh& mesh;
    Tally& tally;

    uint32_t n_rays; 
    uint32_t n_tally_points;

    double mesh_start_x, mesh_start_y;
    double mesh_end_x, mesh_end_y;
    double mesh_dx, mesh_dy, mesh_dz;   

    double tally_x1, tally_y1;
    double tally_x2, tally_y2;

};
#endif
