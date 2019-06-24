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
#include "cell.h"

using namespace std;

    struct Point {
	double x;
	double y;
    };
    struct Ray{
	double m; // Slope of the ray
	double b; // Y-Intercept of the ray
    };
    struct Response_Cell {
	Point l_left; // Lower left point of box
	Point u_right; // Upper right point of box	
	int row; // row in 2D array view of grid
	int col; // col in 2D array view of grid
	double response; // Contribution to cell
    };
    struct Particle {
	Point pos; // Starting position of the particle
	Ray dir; // the line the particle is travelling along
	double weight; // Particle weight
    };

class Response_Funct {

// NOTE: Assumes the source is left of the tally!!!!

// Testing Results:
// 	The method to check if a line (y = mx + b) passes
// 	through a cell (box) is test verified to work as 
// 	expected. assumes the point (0,0) is the lower left corner

public:

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
	mesh_n_x_cells = mesh.get_n_x_cells();
	mesh_n_y_cells = mesh.get_n_y_cells();
	mesh_n_z_cells = mesh.get_n_z_cells();

	tally_x1 = tally.get_x1();
	tally_y1 = tally.get_y1();
	tally_x2 = tally.get_x2();
	tally_y2 = tally.get_y2();
    }   

    ~Response_Funct() {}

    // This uses each cell as a region to accumulate over
    // Add a return statement as a double**
    void get_response(int n_points, int n_angles) {
	// Set up a 'blank table': each response cell has no initial contr.	
	Response_Cell** resp_table = get_empty_resp_table(mesh_n_x_cells, mesh_n_y_cells, mesh_dx, mesh_dy);

	// Generate particles with a given start pt and dir
	Particle* particles = get_initial_particles(n_points, n_angles);

	// track particles updating the cell contr and particle weight
	for(int k = 0; k < n_points*n_angles; k++) {
	    Particle particle = particles[k];   
	
	    bool is_active = true;
	    while(is_active) {
		if(particle.pos.x <= 0 || particle.pos.y <= 0) {
		    is_active = false;
		} else {
		    move_particle(particle, resp_table, mesh_dx, mesh_dy);
		}
	    }
	}
    }


    // Find a way to get the sigma_a / sigma_s values!!!!
    bool move_particle(Particle& particle, Response_Cell** resp_table, double table_dx, double table_dy ) {
	int cur_cell_col = (int)((particle.pos.x - 1e-4) / table_dx); // get the col of the resp_cell the particle is in
	int cur_cell_row = (int)((particle.pos.y - 1e-4) / table_dy); // get the row of the resp_cell the particle is in

	if(cur_cell_col < 0) { cur_cell_col = 0; }
	if(cur_cell_row < 0) { cur_cell_row = 0; }

	Point int_pt = get_intersect_point(resp_table[cur_cell_row][cur_cell_col] , particle); // get the point where the particle intersects the next cell

	double dist = sqrt(pow(particle.pos.x-int_pt.x,2) + pow(particle.pos.y-int_pt.y,2)); // distance the particle traveled
	double sigma_a = 0; // GET THE REAL SIGMA_A VALUE FROM MESH INFO ... 
	double contr = get_contribution(particle, dist, sigma_a);
	resp_table[cur_cell_row][cur_cell_col].response += contr;
	
	particle.pos = int_pt;
	return true;
    }

    Point get_intersect_point(Response_Cell& rc, Particle& particle) {
	// Calculate particle residuals 	
	Ray ray = particle.dir;
	Point l_left = rc.l_left;
	Point u_right = rc.u_right;
	Point l_right = {u_right.x, l_left.y};
	Point u_left = {l_left.x, u_right.y};

	double l_left_resid = l_left.y - ray.m * l_left.x - ray.b;
	double l_right_resid = l_right.y - ray.m * l_right.x - ray.b;
	double u_left_resid = u_left.y - ray.m * u_left.x - ray.b;
	double u_right_resid = u_right.y - ray.m * u_right.x - ray.b;

	double min_resid = min(l_left_resid, min(l_right_resid, min(u_left_resid, u_right_resid)));
	double max_resid = max(l_left_resid, max(l_right_resid, max(u_left_resid, u_right_resid)));
	
	if(min_resid <= 0 && max_resid >= 0) { // ray intersects the cell
	    // Check if the particle is horiz
	    if(ray.m == 0) {
		Point int_pt = {l_left.x , particle.pos.y};
		return int_pt;
	    }

	    double resid_array[4] = {l_left_resid,l_right_resid,u_right_resid,u_left_resid};
	    sort(resid_array, resid_array + 4);
            
	    Point int_pt_1;
	    Point int_pt_2;
 
	    //Case 1: three of the residuals are either positive or negative
	    if(resid_array[0] < 0 && resid_array[1] > 0) {
		if(resid_array[0] == l_left_resid) {
		    int_pt_1 = {l_left.x , ray.m*l_left.x + ray.b};
		    int_pt_2 = {(l_left.y-ray.b)/ray.m , l_left.y};
		} else if(resid_array[0] == l_right_resid) {
		    int_pt_1 = {l_right.x , ray.m*l_right.x + ray.b};
		    int_pt_2 = {(l_right.y-ray.b)/ray.m , l_right.y};
		} else if(resid_array[0] == u_left_resid) {
		    int_pt_1 = {u_left.x , ray.m*u_left.x + ray.b};
		    int_pt_2 = {(u_left.y-ray.b)/ray.m , u_left.y};
		} else  {
		    int_pt_1 = {u_right.x , ray.m*u_right.x + ray.b};
		    int_pt_2 = {(u_right.y-ray.b)/ray.m , u_right.y};
		}
	    } else if (resid_array[2] < 0 && resid_array[3] > 0) {
		if(resid_array[3] == l_left_resid) {
		    int_pt_1 = {l_left.x , ray.m*l_left.x + ray.b};
		    int_pt_2 = {(l_left.y-ray.b)/ray.m , l_left.y};
		} else if(resid_array[3] == l_right_resid) {
		    int_pt_1 = {l_right.x , ray.m*l_right.x + ray.b};
		    int_pt_2 = {(l_right.y-ray.b)/ray.m , l_right.y};
		} else if(resid_array[3] == u_left_resid) {
		    int_pt_1 = {u_left.x , ray.m*u_left.x + ray.b};
		    int_pt_2 = {(u_left.y-ray.b)/ray.m , u_left.y};
		} else  {
		    int_pt_1 = {u_right.x , ray.m*u_right.x + ray.b};
		    int_pt_2 = {(u_right.y-ray.b)/ray.m , u_right.y};
		}
	    }	
	    // Case 2: two points have the same parity
	    else {
		if(l_left_resid <= 0 && l_right_resid <= 0 ) {
		    int_pt_1 = {l_left.x , ray.m*l_left.x - ray.b};
		    int_pt_2 = {l_right.x , ray.m*l_right.x};
		} else if(u_left_resid <= 0 && u_right_resid <= 0) {
		    int_pt_1 = {u_left.x , ray.m*u_left.x - ray.b};
		    int_pt_2 = {u_right.x , ray.m*u_right.x};
		} else if(l_left_resid <= 0 && u_right_resid <= 0) {	
		    int_pt_1 = {l_left.x , ray.m*l_left.x - ray.b};
		    int_pt_2 = {u_right.x , ray.m*u_right.x};
		} else if(l_left_resid <= 0 && u_left_resid <= 0) {
		    int_pt_1 = {l_left.x , ray.m*l_left.x - ray.b};
		    int_pt_2 = {u_left.x , ray.m*u_left.x};
		} else if(l_right_resid <= 0 && u_right_resid <= 0) {
		    int_pt_1 = {l_right.x , ray.m*l_right.x - ray.b};
		    int_pt_2 = {u_right.x , ray.m*u_right.x};
		} else { // lower right and upper left 
		    int_pt_1 = {l_right.x , ray.m*l_right.x - ray.b};
		    int_pt_2 = {u_left.x , ray.m*u_left.x};
		}	
	    }
	    return (int_pt_1.x < particle.pos.x) ? int_pt_1 : int_pt_2;
	}
	return {NULL , NULL}; // Indicate that the ray does not intersect the cell
    }

    double get_contribution(Particle& particle, double dist, double sigma_a) {
	double absorbed_E = particle.weight * (1 - exp(-sigma_a * dist));
	return absorbed_E;
    }

    Response_Cell** get_empty_resp_table(int n_cols, int n_rows, double dx, double dy) {
	Response_Cell** resp_table;
	resp_table = new Response_Cell*[n_rows];
	for(int row = 0; row < n_cols; row++) {
	    resp_table[row] = new Response_Cell[n_cols];
	    for(int col = 0; col < n_cols; col++) { 
		Point l_left = {col * dx , row * dy};
		Point u_right = {(col+1) * dx , (row + 1) * dy};
		resp_table[row][col] = {l_left, u_right, row, col, 0.0};
	    }
	}
	return resp_table;
    }
	
    Particle* get_initial_particles(int n_points, int n_angles) {
	Particle* particles = new Particle[n_angles * n_points];
	
	if(tally_x1 > tally_x2) {
	    if(tally_y1 > tally_y2) {
		// Case: tally_1 is poitively greater than tally_2
		double pos_dx = (tally_x1 - tally_x2) / n_points;
		double pos_dy = (tally_y1 - tally_y2) / n_points;
		double d_angle = (max_slope - min_slope) / n_angles; 
		
		int cur_index = 0;
		for(int pos = 0; pos < n_points; pos++) {
		Point pt = {tally_x2 + pos*pos_dx , tally_y2 + pos*pos_dy};
		    for(int angle = 0; angle < n_angles; angle++) {
			double slope = min_slope + (angle * d_angle);
			Ray ray = {slope, -slope * pt.x + pt.y};
			Particle particle = {pt , ray, particle_weight};
			particles[cur_index] = particle;
			cur_index++;
		    }
		}
	    } else {
		// Case: tally_1 is negatively greater than tally_2  
		double pos_dx = (tally_x1 - tally_x2) / n_points;
		double pos_dy = (tally_y2 - tally_y1) / n_points;
		double d_angle = (max_slope - min_slope) / n_angles; 
		
		int cur_index = 0;
		for(int pos = 0; pos < n_points; pos++) {
		Point pt = {tally_x2 + pos*pos_dx , tally_y1 + pos*pos_dy};
		    for(int angle = 0; angle < n_angles; angle++) {
			double slope = min_slope + (angle * d_angle);
			Ray ray = {slope, -slope * pt.x + pt.y};
			Particle particle = {pt , ray, particle_weight};
			particles[cur_index] = particle;
			particles[cur_index] = particle;
			cur_index++;
		    }
		}
	    }
	} else {
	    if(tally_y1 > tally_y2) {
		// Case: tally_2 is positively greater than tally_1
		double pos_dx = (tally_x2 - tally_x1) / n_points;
		double pos_dy = (tally_y2 - tally_y1) / n_points;
		double d_angle = (max_slope - min_slope) / n_angles; 
		
		int cur_index = 0;
		for(int pos = 0; pos < n_points; pos++) {
		Point pt = {tally_x1 + pos*pos_dx , tally_y1 + pos*pos_dy};
		    for(int angle = 0; angle < n_angles; angle++) {
			double slope = min_slope + (angle * d_angle);
			Ray ray = {slope, -slope * pt.x + pt.y};
			Particle particle = {pt , ray, particle_weight};
			particles[cur_index] = particle;
			cur_index++;
		    }
		}
	    } else {
		// Case: tally_2 is negatively greater than tally_2
		double pos_dx = (tally_x2 - tally_x1) / n_points;
		double pos_dy = (tally_y1 - tally_y2) / n_points;
		double d_angle = (max_slope - min_slope) / n_angles; 
		
		int cur_index = 0;
		for(int pos = 0; pos < n_points; pos++) {
		Point pt = {tally_x1 + pos*pos_dx , tally_y2 + pos*pos_dy};
		    for(int angle = 0; angle < n_angles; angle++) {
			double slope = min_slope + (angle * d_angle);
			Ray ray = {slope, -slope * pt.x + pt.y};
			Particle particle = {pt , ray, particle_weight};
			particles[cur_index] = particle;
			cur_index++;
		    }
		}
	    }
	}
	
	return particles;
    }

private:

    Mesh& mesh;
    Tally& tally;

    uint32_t n_rays; 
    uint32_t n_tally_points;

    double mesh_start_x, mesh_start_y;
    double mesh_end_x, mesh_end_y;
    double mesh_dx, mesh_dy, mesh_dz;   
    int mesh_n_x_cells, mesh_n_y_cells, mesh_n_z_cells;

    double tally_x1, tally_y1;
    double tally_x2, tally_y2;
	
    double min_slope = -10;
    double max_slope = 10;
    double particle_weight = 100; 

};
#endif
