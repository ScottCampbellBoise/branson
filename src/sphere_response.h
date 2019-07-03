/*
 * @file:    response_funct.h
 * @author:  Scott Campbell
 * @version: 18 June 2019
 * @desc:    This class stores/holds the response value for a each
 *           cell in a grid system.
 */

#ifndef sphere_response_h
#define sphere_response_h

#include <mpi.h>

#include "RNG.h"
#include "constants.h"
#include "info.h"
#include "mesh.h"
#include "photon.h"
#include "tally.h"

using namespace std;
using Constants::ELEMENT;
using Constants::bc_type;
using Constants::event_type;
using Constants::EXIT;
using Constants::KILL;

class Sphere_Response {

public:

    // Constructor: store the info about the mesh and tally
    Sphere_Response(Tally& tally, Mesh& mesh, IMC_State &imc_state) : tally(tally),
					     mesh(mesh), imc_state(imc_state) {}   

    ~Sphere_Response() {}

    // generate photon objects starting on the tally surface
    // track the photons outward,
    //     each cell the photon passes through, update the dist and sigma_dist
    //     if the particle reaches a barrier - kill it
    void generate_response(int n_particles) {
	if(!response_set) {
	    setup_response();
	    response_set = true;
	}

	phtn_total_sigma_dist = new double[n_particles];
	phtn_total_dist = new double[n_particles];
	
	for(int k = 0; k < n_particles; k++) {
	    Photon phtn_1, phtn_2;
            create_photon(phtn_1);
	    phtn_2 = phtn_1; // Make a copy of the phtn
	
	    phtn_total_dist[k] = 0.0;
	    phtn_total_sigma_dist[k] = 0.0;
   		
	    move_photon(phtn_1, k, true); 
	    move_photon(phtn_2, k, false);
	}

	uint32_t n_cells = mesh.get_n_local_cells();
	cell_response = new double[n_cells];
	
	for(int k = 0; k < n_cells; k++) {
	    cell_response[k] = (cell_total_sigma_dist[k] / cell_total_dist[k]);
	}
    }

    void reset_response() {
	response_set = false;	

	// Iterate through every cell in the mesh and reset its value	
 	uint32_t n_cell = mesh.get_n_local_cells();
	for(uint32_t k = 0; k < n_cell; ++k) {
	    cell_total_sigma_dist[k] = 0.0;
	    cell_total_dist[k] = 0.0;
	}
    }

    void create_photon(Photon& phtn) {
        // Create a photon on the tally surface	
        double phi = 2 * Constants::pi * rng->generate_random_number();
        double mu = 1 - 2 * rng->generate_random_number();
        double theta = acos(mu);
        double pos[3] = {tally.get_x1() + tally.get_radius()*(cos(phi)*sqrt((1-pow(mu,2)))),
       	     	         tally.get_y1() + tally.get_radius()*(sin(phi)*sqrt(1-pow(mu,2))),
        		 tally.get_z1() + tally.get_radius()*mu};

        // Cosine-distribution for angle
        double mu_r = sqrt(rng->generate_random_number());
        double phi_r = 2 * Constants::pi * rng->generate_random_number();
        double mu_theta = cos(phi_r) * sqrt(1 - pow(mu,2));
        double mu_phi = sin(phi_r) * sqrt(1 - pow(mu_r,2));
       
        double mu_x = sin(theta)*cos(phi)*mu_r + cos(theta)*cos(phi)*mu_theta
        	      - sin(phi)*mu_phi;
        double mu_y = sin(theta)*sin(phi)*mu_r + cos(theta)*sin(phi)*mu_theta
        	      + cos(phi)*mu_phi;
        double mu_z = cos(theta)*mu_r - sin(theta)*mu_theta;
        double angle[3] = {mu_x, mu_y, mu_z};
        
        // Set photon values
        phtn.set_position(pos);
        phtn.set_angle(angle);
        phtn.set_cell(get_photons_cell(pos[0],pos[1],pos[2]));
        phtn.set_group(floor(rng->generate_random_number() * double(BRANSON_N_GROUPS)));
    }

    void print_response() {
	ofstream outfile("RESP_TABLE.txt");

 	uint32_t n_cell = mesh.get_n_local_cells();
	for(uint32_t k = 0; k < n_cell; ++k) {
	    outfile <<  "\tcell resp: " << cell_response[k] << endl;
	}
	
	outfile.close();
    }

private:

    void setup_response() {
	tally_r = tally.get_radius();
	tally_x = tally.get_x1();
	tally_y = tally.get_y1();
	tally_z = tally.get_z1();

	// To set the cell opacity values
        mesh.calculate_photon_energy(imc_state);
	
	// Set the random number generator
	rng = new RNG();

        int cur_pos = 0;

 	uint32_t n_cell = mesh.get_n_local_cells();
	Cell temp_tally_cells[n_cell]; 

	cell_total_sigma_dist = new double[n_cell];
	cell_total_dist = new double[n_cell];
	
	for(uint32_t k = 0; k < n_cell; ++k) {

	    cell_total_sigma_dist[k] = 0.0;
	    cell_total_dist[k] = 0.0;

	    Cell cell = mesh.get_cell(k);
	    if(tally_intersects_cell(cell)) {
		temp_tally_cells[cur_pos++] = cell;
	    }
	}

	tally_cells = new Cell[cur_pos];
	for(int k = 0; k < cur_pos; k++) {
	    tally_cells[k] = temp_tally_cells[k];
	}	

	n_tally_cells = cur_pos + 1;
    }

    void move_photon(Photon &phtn, int phtn_pos, bool track_phtn) {
	uint32_t cell_id, next_cell;
	bc_type boundary_event;
	double dist_to_event;
	double sigma_a;
	double angle[3];

        uint32_t surface_cross = 0;
  	
        int group;
  	Cell cell;
  	cell_id = phtn.get_cell();
  	cell = mesh.get_cell(cell_id);
	
  	bool active = true;
  	while (active) {
            group = phtn.get_group();
            sigma_a = cell.get_op_a(group);
       
	    // get distance to event
            dist_to_event = cell.get_distance_to_boundary(
                phtn.get_position(), phtn.get_angle(), surface_cross);
       
	    if(track_phtn) {
		phtn_total_dist[phtn_pos] += dist_to_event;
	    	phtn_total_sigma_dist[phtn_pos] += (dist_to_event * sigma_a);
	    } else {
		cell_total_dist[cell_id] += dist_to_event;
	    	cell_total_sigma_dist[cell_id] += 
		    (phtn_total_sigma_dist[phtn_pos] / phtn_total_dist[phtn_pos]) * dist_to_event;
	    }
 
 	    // update position
            phtn.move(dist_to_event);

            boundary_event = cell.get_bc(surface_cross);
            if (boundary_event == ELEMENT) {
                next_cell = cell.get_next_cell(surface_cross);
                phtn.set_cell(next_cell);
                cell_id = next_cell;
                cell = mesh.get_cell(cell_id);
            } else {
	        active = false;
            }
        }   // end while alive
    }

    // Check if a point is inside one of the cells
    // that the tally surface is in
    uint32_t get_photons_cell(double x, double y, double z) {
	for(int k = 0; k < n_tally_cells; k++) {
	    Cell cell = tally_cells[k];
	    const double* cell_dim = cell.get_node_array();
	    if(x >= cell_dim[0] && x <= cell_dim[1] &&
	       y >= cell_dim[1] && y <= cell_dim[2] &&
	       z >= cell_dim[3] && z <= cell_dim[5]) {
		return tally_cells[k].get_ID();
	    } 
	}
	return tally_cells[0].get_ID();
    }

    // Check if the tally surface intersects the cell
    bool tally_intersects_cell(Cell& cell) {
	const double* cell_dim = cell.get_node_array();
        double dist_squared = pow(tally_r,2);
	
	if(tally_x < cell_dim[0])
	    dist_squared -= pow(tally_x - cell_dim[0], 2);
	else if(tally_x > cell_dim[1]) 
	    dist_squared -= pow(tally_x - cell_dim[1], 2);
	if(tally_y < cell_dim[2])
	    dist_squared -= pow(tally_y - cell_dim[2], 2);
	else if(tally_y > cell_dim[3]) 
	    dist_squared -= pow(tally_y - cell_dim[3], 2);
	if(tally_z < cell_dim[4])
	    dist_squared -= pow(tally_z - cell_dim[4], 2);
	else if(tally_z > cell_dim[5]) 
	    dist_squared -= pow(tally_z - cell_dim[5], 2);
	return dist_squared > 0;
    }

    //-------------------------------------------
    // Variables
    //-------------------------------------------

    bool response_set = false;

    Mesh& mesh;
    Tally& tally;
    IMC_State &imc_state;

    Cell* tally_cells;
    int n_tally_cells;

    RNG* rng;

    double tally_r, tally_x, tally_y, tally_z;

    double* cell_total_sigma_dist;
    double* cell_total_dist;
    double* phtn_total_sigma_dist;
    double* phtn_total_dist;

    double* cell_response;
};
#endif
