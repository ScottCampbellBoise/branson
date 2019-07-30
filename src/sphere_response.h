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
	
#include "RNG.h"
#include "constants.h"
#include "info.h"
#include "mesh.h"
#include "photon.h"
#include "tally.h"
#include "response_exception.h"

using namespace std;
using Constants::ELEMENT;
using Constants::bc_type;
using Constants::event_type;
using Constants::EXIT;
using Constants::KILL;

class Sphere_Response {
public:

    // Constructor: store the info about the mesh and tally
    Sphere_Response(Tally*& tally, const Mesh& mesh, IMC_State &imc_state, uint32_t n_particles) : tally(tally),
					     mesh(mesh), imc_state(imc_state), n_particles(n_particles) {}   

    ~Sphere_Response() {}

    // generate photon objects starting on the tally surface
    // track the photons outward,
    //     each cell the photon passes through, update the dist and sigma_dist
    //     if the particle reaches a barrier - kill it
    void generate_response() {

	if(!response_set) {
	    setup_response();
	    response_set = true;
	}
	
	cell_total_sigma_dist.assign(n_cell, 0.0);
        cell_total_dist.assign(n_cell, 0.0);
        cell_xp_total_sigma_dist.assign(n_cell, 0.0);
        cell_xp_total_dist.assign(n_cell, 0.0);
        cell_xm_total_sigma_dist.assign(n_cell, 0.0);
        cell_xm_total_dist.assign(n_cell, 0.0);
        cell_yp_total_sigma_dist.assign(n_cell, 0.0);
        cell_yp_total_dist.assign(n_cell, 0.0);
        cell_ym_total_sigma_dist.assign(n_cell, 0.0);
        cell_ym_total_dist.assign(n_cell, 0.0);
        cell_zp_total_sigma_dist.assign(n_cell, 0.0);
        cell_zp_total_dist.assign(n_cell, 0.0);
        cell_zm_total_sigma_dist.assign(n_cell, 0.0);
        cell_zm_total_dist.assign(n_cell, 0.0);
        num_respon_source.assign(n_cell, 0);

        uint32_t index;
	double pos[3];
        uint32_t cell_id;
	Photon phtn;

	for(int k = 0; k < n_particles; k++) {
	    index = (uint32_t)(rng->generate_random_number() * phtn_deck_size);
	    pos[0] = start_x[index];
	    pos[1] = start_y[index];
	    pos[2] = start_z[index];	
	    cell_id = start_cell_id[index];    
	    num_respon_source[cell_id]++;	    

	    create_photon(phtn, pos, cell_id);
	    move_photon(phtn); 
	}

	response_generated = true;
    }

    inline bool get_response_state() const { return response_generated; }

    void increase_response() {
        uint32_t index;
	double pos[3];
        uint32_t cell_id;
    Photon phtn;

    for(int k = 0; k < n_particles; k++) {
        index = (uint32_t)(rng->generate_random_number() * phtn_deck_size);
        pos[0] = start_x[index];
        pos[1] = start_y[index];
        pos[2] = start_z[index];    
        cell_id = start_cell_id[index];        

        create_photon(phtn, pos, cell_id);
        move_photon(phtn); 
    }

	response_generated = true;
    }
 
    double get_angle_response(uint32_t cell_id, const double angle[3]) const throw(Response_Exception){

    double xp = max(angle[0]*cell_xp_total_sigma_dist[cell_id]/cell_xp_total_dist[cell_id],0.0);
    double xm = max(-angle[0]*cell_xm_total_sigma_dist[cell_id]/cell_xm_total_dist[cell_id],0.0);
    double yp = max(angle[1]*cell_yp_total_sigma_dist[cell_id]/cell_yp_total_dist[cell_id],0.0);
    double ym = max(-angle[1]*cell_ym_total_sigma_dist[cell_id]/cell_ym_total_dist[cell_id],0.0);
    double zp = max(angle[2]*cell_zp_total_sigma_dist[cell_id]/cell_zp_total_dist[cell_id],0.0);
    double zm = max(-angle[2]*cell_zm_total_sigma_dist[cell_id]/cell_zm_total_dist[cell_id],0.0);

    double resp =  xp+xm+yp+ym+zp+zm;

    //double resp = cell_total_sigma_dist[cell_id] / cell_total_dist[cell_id];
    if(resp <= 0 || isnan(resp))
        throw Response_Exception();
    return resp; 
    }

    double get_response(uint32_t cell_id) const throw(Response_Exception){
	double resp = cell_total_sigma_dist[cell_id] / cell_total_dist[cell_id];
	if(resp <= 0 || isnan(resp))
	    throw Response_Exception();
	return resp; 
    }

    double get_dist(uint32_t cell_id) const { return cell_total_dist[cell_id]; }
    double get_n_sourced(uint32_t cell_id) const { return num_respon_source[cell_id]; }

    void reset_response() { response_set = false; }

    void create_photon(Photon& phtn, double pos[3], uint32_t cell_id) {
	double angle[3];
        get_uniform_angle(pos, angle);

        // Set photon values
        phtn.set_total_dist(0.0);
	phtn.set_total_sigma_dist(0.0);
        phtn.set_position(pos);
        phtn.set_angle(angle);
        phtn.set_cell(cell_id);
        phtn.set_group(floor(rng->generate_random_number() * double(BRANSON_N_GROUPS)));
    }

    void get_uniform_angle(double pos[3], double angle[3]) {
        // Cosine-distribution for angle
        double costheta = sqrt(rng->generate_random_number());
	double phi = 2.0 * Constants::pi * rng->generate_random_number();

	// calculate the inward normal.
	double r = sqrt(pow(tally_x-pos[0],2)+pow(tally_y-pos[1],2)+pow(tally_z-pos[2],2));
	double x_hat = (tally_x - pos[0]) / r;
	double y_hat = (tally_y - pos[1]) / r;
	double z_hat = (tally_z - pos[2]) / r;

	// helpful terms
	double sintheta = sqrt(1.0 - costheta*costheta);
	double factor = sqrt(fabs(1.0-z_hat*z_hat));
	double cosphi = cos(phi);
	double sinphi = sin(phi);
	double sintcosp = sintheta * cosphi;
	double sintsinp = sintheta * sinphi;
	double f_inv = 1.0/factor;
	double stcpDf = sintcosp*f_inv;
	double stspDf = sintsinp*f_inv;

	// scatter through the random angle
	if(factor<1.0e-6) {
	    angle[0] = sintcosp;
	    angle[1] = sintsinp;
	    angle[2] = ((z_hat<0)? - 1.0 : 1.0)*costheta;
	} else {
	    angle[0] = x_hat*costheta + z_hat*x_hat*stcpDf - y_hat*stspDf;
	    angle[1] = y_hat*costheta + z_hat*y_hat*stcpDf + x_hat*stspDf;
	    angle[2] = z_hat*costheta - factor*sintcosp;
	}	

	double norm = sqrt(angle[0]*angle[0]+angle[1]*angle[1]+angle[2]*angle[2]);
	angle[0] /= norm;
	angle[1] /= norm;
	angle[2] /= norm;

	// just set the direction in the outward normal procedures the analytic result
	/*
	angle[0] = -(tally_x - pos[0]) / r;
	angle[1] = -(tally_y - pos[1]) / r;
	angle[2] = -(tally_z - pos[2]) / r;
	*/
    }

    void get_start_pos(double pos[3]) {
        double phi = 2 * Constants::pi * rng->generate_random_number();
        double mu = 1 - 2 * rng->generate_random_number();
	double theta = acos(mu);       

	pos[0] = tally_x + tally_r*(cos(phi)*sqrt((1-pow(mu,2))));
	pos[1] = tally_y + tally_r*(sin(phi)*sqrt((1-pow(mu,2))));
        pos[2] = tally_z + tally_r*mu;
    }

    void print_response() {
	ofstream outfile("RESP_TABLE.txt");

 	uint32_t n_cell = mesh.get_n_local_cells();
	for(uint32_t k = 0; k < n_cell; ++k) {
	    outfile << "\tcell resp: " << cell_total_sigma_dist[k] / cell_total_dist[k] << endl;
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

	// Find all cells that the tally intersects
 	n_cell = mesh.get_n_local_cells();
	tally_cells.resize(n_cell);
        uint32_t cur_pos = 0;

	for(uint32_t k = 0; k < n_cell; ++k) {
	    Cell cell = mesh.get_cell(k);
	    if(tally_intersects_cell(cell)) {
		tally_cells[cur_pos++] = cell.get_ID();
	    }
	}

	tally_cells.shrink_to_fit();

	//Create a 'deck' of photons to use for the sampling
	start_x.resize(phtn_deck_size);
	start_y.resize(phtn_deck_size);
 	start_z.resize(phtn_deck_size);
	start_cell_id.resize(phtn_deck_size);

        double pos[3];

	for(uint32_t k = 0; k < phtn_deck_size; k++) {
	    get_start_pos(pos);

	    start_x[k] = pos[0];
	    start_y[k] = pos[1];
	    start_z[k] = pos[2];
	    start_cell_id[k] = get_photons_cell(pos[0], pos[1], pos[2]);
	}

	//Generate the variables to hold the resp info
	cell_total_sigma_dist.resize(n_cell);
        cell_total_dist.resize(n_cell);	
   	num_respon_source.resize(n_cell);
    }

    void move_photon(Photon &phtn) {
	uint32_t cell_id, next_cell;
	double dist_to_event;
	double sigma_a;
	const double* angle = phtn.get_angle();

        uint32_t surface_cross = 0;
  	
  	Cell cell;
  	cell_id = phtn.get_cell();
  	cell = mesh.get_cell(cell_id);
	
  	bool active = true;
  	while (active) {
            sigma_a = cell.get_op_a(phtn.get_group());

	    // get distance to event
            dist_to_event = cell.get_distance_to_boundary(
                phtn.get_position(), phtn.get_angle(), surface_cross);
       
        phtn.add_to_total_dist(dist_to_event, sigma_a);
         cell_total_dist[cell_id] += dist_to_event;
        
        double ave_sig = phtn.get_total_sigma_dist() / phtn.get_total_dist();
        cell_total_sigma_dist[cell_id] += 
            (ave_sig) * dist_to_event;
        
        double xp = max(-angle[0]*dist_to_event,0.0);
        double xm = max(angle[0]*dist_to_event,0.0);
        double yp = max(-angle[1]*dist_to_event,0.0);
        double ym = max(angle[1]*dist_to_event,0.0);
        double zp = max(-angle[2]*dist_to_event,0.0);
        double zm = max(angle[2]*dist_to_event,0.0);


         cell_xp_total_dist[cell_id] += xp;
         cell_xm_total_dist[cell_id] += xm;
         cell_yp_total_dist[cell_id] += yp;
         cell_ym_total_dist[cell_id] += ym;
         cell_zp_total_dist[cell_id] += zp;
         cell_zm_total_dist[cell_id] += zm;

         cell_xp_total_sigma_dist[cell_id] += ave_sig*xp;
         cell_xm_total_sigma_dist[cell_id] += ave_sig*xm;
         cell_yp_total_sigma_dist[cell_id] += ave_sig*yp;
         cell_ym_total_sigma_dist[cell_id] += ave_sig*ym;
         cell_zp_total_sigma_dist[cell_id] += ave_sig*zp;
         cell_zm_total_sigma_dist[cell_id] += ave_sig*zm;

    
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
    }

    // Check if a point is inside one of the cells
    // that the tally surface is in
    uint32_t get_photons_cell(double x, double y, double z) {
    for(uint32_t k = 0; k < mesh.get_n_local_cells(); k++) {
    //for(auto k = tally_cells.begin(); k != tally_cells.end(); k++) {
        Cell cell = mesh.get_cell(k);
        const double* cell_dim = cell.get_node_array();
        if(x >= cell_dim[0] && x <= cell_dim[1] &&
           y >= cell_dim[2] && y <= cell_dim[3] &&
           z >= cell_dim[4] && z <= cell_dim[5]) {
        return k;
        } 
    }
    // Don't allow us to not find the correct cell
    std::cout<<"ERROR: Possition non found -> "<<x<<" "<<y<<" "<<z<<std::endl;
    exit(0);
    }

    // Check if the tally surface intersects the cell
    bool tally_intersects_cell(Cell& cell) {
	const double* pos = cell.get_node_array();
   	// get the distance from all six points to the center
   	double dist1 = sqrt(pow(pos[0]-tally_x,2)+pow(pos[2]-tally_y,2)+pow(pos[4]-tally_z,2));	
   	double dist2 = sqrt(pow(pos[1]-tally_x,2)+pow(pos[2]-tally_y,2)+pow(pos[4]-tally_z,2));	
   	double dist3 = sqrt(pow(pos[1]-tally_x,2)+pow(pos[3]-tally_y,2)+pow(pos[4]-tally_z,2));	
   	double dist4 = sqrt(pow(pos[0]-tally_x,2)+pow(pos[3]-tally_y,2)+pow(pos[4]-tally_z,2));	
   	double dist5 = sqrt(pow(pos[0]-tally_x,2)+pow(pos[2]-tally_y,2)+pow(pos[5]-tally_z,2));	
   	double dist6 = sqrt(pow(pos[1]-tally_x,2)+pow(pos[2]-tally_y,2)+pow(pos[5]-tally_z,2));	
   	double dist7 = sqrt(pow(pos[0]-tally_x,2)+pow(pos[3]-tally_y,2)+pow(pos[5]-tally_z,2));	
   	double dist8 = sqrt(pow(pos[1]-tally_x,2)+pow(pos[3]-tally_y,2)+pow(pos[5]-tally_z,2));	

	return !((dist1 > tally_r && dist2 > tally_r && dist3 > tally_r && dist4 > tally_r &&
	     dist5 > tally_r && dist6 > tally_r && dist7 > tally_r && dist8 > tally_r) ||
            (dist1 < tally_r && dist2 < tally_r && dist3 < tally_r && dist4 < tally_r &&
	     dist5 < tally_r && dist6 < tally_r && dist7 < tally_r && dist8 < tally_r));
    }

    //-------------------------------------------
    // Variables
    //-------------------------------------------

    const Mesh& mesh;
    Tally*& tally;
    IMC_State &imc_state;
    RNG* rng;
    
    uint32_t n_particles;
    uint32_t n_cell;
    double tally_r, tally_x, tally_y, tally_z;
    
    bool response_set = false;
    bool response_generated = false;

    vector<uint32_t> tally_cells;

    vector<uint32_t> num_respon_source;
    vector<double> cell_total_sigma_dist;
    vector<double> cell_xp_total_sigma_dist;
    vector<double> cell_xm_total_sigma_dist;
    vector<double> cell_yp_total_sigma_dist;
    vector<double> cell_ym_total_sigma_dist;
    vector<double> cell_zp_total_sigma_dist;
    vector<double> cell_zm_total_sigma_dist;

    vector<double> cell_total_dist;
    vector<double> cell_xp_total_dist;
    vector<double> cell_xm_total_dist;
    vector<double> cell_yp_total_dist;
    vector<double> cell_ym_total_dist;
    vector<double> cell_zp_total_dist;
    vector<double> cell_zm_total_dist;

    uint32_t phtn_deck_size = 100000;
    vector<double> start_x;
    vector<double> start_y;
    vector<double> start_z;
    vector<uint32_t> start_cell_id;
};
#endif
