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
#include <iostream>
#include <fstream>
#include <math.h>
//#include <chrono>
	
#include "RNG.h"
#include "constants.h"
#include "info.h"
#include "mesh.h"
#include "photon.h"
#include "tally.h"

using namespace std::chrono;

using namespace std;
using Constants::ELEMENT;
using Constants::bc_type;
using Constants::event_type;
using Constants::EXIT;
using Constants::KILL;

class Sphere_Response {

public:

    // Constructor: store the info about the mesh and tally
    Sphere_Response(Tally*& tally, const Mesh& mesh, IMC_State &imc_state) : tally(tally),
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

	for(uint32_t k = 0; k < n_cell; k++) {
	    cell_total_sigma_dist[k] = 0;
            cell_total_dist[k] = 0;
	}

	Photon phtn;
	double total_phtn_dist;
	double total_phtn_sigma_dist;
 	double sigma_a;
	double cell_dist;
	Cell cell;
	    
	for(int k = 0; k < n_particles; k++) {
	    phtn = photon_deck[(int)(rng->generate_random_number() * photon_deck_size)];
	    total_phtn_dist = 0;
	    total_phtn_sigma_dist = 0;
		
	    for(uint32_t index = 0; index < phtn.get_n_resp_cells(); index++) {
 	        cell_dist = phtn.get_resp_cell_dist(index);
		cell = mesh.get_cell(phtn.get_resp_cell(index));
		sigma_a = cell.get_op_a(phtn.get_group());
		    
		total_phtn_dist += cell_dist;
		total_phtn_sigma_dist += sigma_a * cell_dist;
	        cell_total_dist += cell_dist;
		cell_total_sigma_dist += 
		    (total_phtn_sigma_dist / total_phtn_dist) * cell_dist;
	    }
	}
	response_generated = true;
    }

    double get_response(uint32_t cell_id) {
	double resp = cell_total_sigma_dist[cell_id] / cell_total_dist[cell_id];
	if(resp <= 0 || isnan(resp))
	    return std::numeric_limits<int>::max();
	return resp; 
    }

    void reset_response() { response_set = false; }

    void create_photon(Photon& phtn) {
        // Create a photon on the tally surface	
        double phi = 2 * Constants::pi * rng->generate_random_number();
        double mu = 1 - 2 * rng->generate_random_number();
        double theta = acos(mu);
        double pos[3] = {tally_x + tally_r*(cos(phi)*sqrt((1-pow(mu,2)))),
       	     	         tally_y + tally_r*(sin(phi)*sqrt(1-pow(mu,2))),
        		 tally_z + tally_r*mu};

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

	// Direct the photon into the tally       
	double temp[3] = {pos[0] + angle[0] * tally_r,
		      	  pos[1] + angle[1] * tally_r,
		       	  pos[2] + angle[2] * tally_r};
	double pos_to_center = sqrt(pow(tally_x-pos[0],2) + pow(tally_y-pos[1],2) + pow(tally_z-pos[2],2));
	double temp_to_center = sqrt(pow(tally_x-temp[0],2) + pow(tally_y-temp[1],2) + pow(tally_z-temp[2],2));	
	if(temp_to_center > pos_to_center) {
	    angle[0] *= -1;
	    angle[1] *= -1;
	    angle[2] *= -1;
	} 

        // Set photon values
        phtn.set_total_dist(0.0);
	phtn.set_total_sigma_dist(0.0);
        phtn.set_position(pos);
        phtn.set_angle(angle);
        phtn.set_cell(get_photons_cell(pos[0],pos[1],pos[2]));
        phtn.set_group(floor(rng->generate_random_number() * double(BRANSON_N_GROUPS)));
    }

    void print_response() {
	ofstream outfile("RESP_TABLE.txt");

 	uint32_t n_cell = mesh.get_n_local_cells();
	for(uint32_t k = 0; k < n_cell; ++k) {
	    outfile <<  "\tcell resp: " << cell_total_sigma_dist[k] / cell_total_dist[k] << endl;
	}
	
	outfile.close();
    }

private:

    void setup_response() {
	tally_r = tally->get_radius();
	tally_x = tally->get_x1();
	tally_y = tally->get_y1();
	tally_z = tally->get_z1();

	// Set the random number generator
	rng = new RNG();

        int cur_pos = 0;

	// Find all cells that the tally intersects
 	n_cell = mesh.get_n_local_cells();
	Cell temp_tally_cells[n_cell]; 

	for(uint32_t k = 0; k < n_cell; ++k) {
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

	//Create a 'deck' of photons to use for the sampling
	photon_deck = new Photon[photon_deck_size];
	for(int k = 0; k < photon_deck_size; k++) {
	    Photon phtn;
	    create_photon(phtn);
	    move_photon(phtn);
	    photon_deck[k] = phtn;
	}

	//Generate the variables to hold the resp info
	cell_total_sigma_dist = new double[n_cell];
        cell_total_dist = new double[n_cell];	
    }

    void move_photon(Photon &phtn) {
	uint32_t cell_id, next_cell;
	double dist_to_event;
	double sigma_a;

        uint32_t surface_cross = 0;
  	
  	Cell cell;
  	cell_id = phtn.get_cell();
  	cell = mesh.get_cell(cell_id);
	
	double* cell_ids = new double[n_cell];
	double* cell_dist = new double[n_cell];
	uint32_t pos = 0;
	    
  	bool active = true;
  	while (active) {
            sigma_a = cell.get_op_a(phtn.get_group());

/*	   if(!tally->is_inside_tally(phtn)) {
		active = false;
	    }      
*/
	    // get distance to event
            dist_to_event = cell.get_distance_to_boundary(
                phtn.get_position(), phtn.get_angle(), surface_cross);
   
	    cell_ids[pos] = cell_id;
	    cell_dist[pos] = dist_to_event;
	    pos++;	
	
	    // update position
            phtn.move(dist_to_event);

	    if(active && cell.get_bc(surface_cross) == ELEMENT) {
                next_cell = cell.get_next_cell(surface_cross);
                phtn.set_cell(next_cell);
                cell_id = next_cell;
                cell = mesh.get_cell(cell_id);
            } else {
	        active = false;
            }
         }   // end while alive
	    
	 //Truncate the arrays
	 double* trunc_ids = new double[pos];
	 double* truc_dists = new double[pos];
	 for(uint32_t k = 0; k < pos; ++k) {
	     trunc_ids[k] = cell_ids[k];
	     trunc_dists[k] = cell_dist[k];
	 }
	 phtn.set_resp_cell_ids(trunc_ids);
	 phtn.set_resp_cell_dist(trunc_dists);
	 phtn.set_n_resp_cells(pos);
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
    bool response_generated = false;

    const Mesh& mesh;
    Tally*& tally;
    IMC_State &imc_state;

    Cell* tally_cells;
    int n_tally_cells;
    uint32_t n_cell;

    RNG* rng;

    double tally_r, tally_x, tally_y, tally_z;

    double* cell_total_sigma_dist;
    double*  cell_total_dist;

    int photon_deck_size = 1000000;
    Photon* photon_deck;
};
#endif
