#ifndef response_transport_h_
#define resposne_transport_h_

#include <algorithm>
#include <functional>
#include <iostream>
#include <mpi.h>
#include <numeric>
#include <vector>

#include "RNG.h"
#include "constants.h"
#include "info.h"
#include "mesh.h"
#include "photon.h"
#include "sampling_functions.h"
#include "source.h"
#include "tally.h"
#include "sphere_response.h"

Constants::event_type resp_transport_photon(Photon &phtn, const Mesh &mesh, RNG *rng,
                                       double &next_dt, double &exit_E,
                                       double &census_E,
                                       std::vector<double> &rank_abs_E,
                                       std::vector<double> &rank_track_E,
				       Tally* tally, Sphere_Response* resp) {
  using Constants::ELEMENT;
  using Constants::REFLECT;
  using Constants::VACUUM;
  // events
  using Constants::bc_type;
  using Constants::c;
  using Constants::CENSUS;
  using Constants::event_type;
  using Constants::EXIT;
  using Constants::KILL;
  using std::min;

  uint32_t cell_id, next_cell;
  bc_type boundary_event;
  event_type event;
  double dist_to_scatter, dist_to_boundary, dist_to_census, dist_to_event;
  double sigma_a, sigma_s, f, absorbed_E, ew_factor;
  double angle[3];
  int group;
  Cell cell;

  uint32_t surface_cross = 0;
  double cutoff_fraction = 0.01; // note: get this from IMC_state

  cell_id = phtn.get_cell();
  cell = mesh.get_on_rank_cell(cell_id);
  bool active = true;



  // generate the response table
  resp->generate_response(1000);
 


  // transport this photon
  while (active) {
    group = phtn.get_group();
    sigma_a = cell.get_op_a(group);
    sigma_s = cell.get_op_s(group);
    f = cell.get_f();


    //------------------------------------------------------------------------------------------
    // Added methods for the sphere response method
    //------------------------------------------------------------------------------------------

    // Get the distance of the photon from the tally surface
    double dist_to_tally = tally->get_dist_to_tally(phtn);

    // calculate the contribution to the tally 
    double tally_contr = phtn.get_E() * 
	exp(-(resp->get_response(cell_id) + 1 / phtn.get_distance_remaining()) * dist_to_tally);   
  
    // Add the contribution to the tally
    tally->add_weight(tally_contr);	

    //------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------


    // get distance to event
    dist_to_scatter =
        -log(rng->generate_random_number()) / ((1.0 - f) * sigma_a + sigma_s);

    dist_to_boundary = cell.get_distance_to_boundary(
        phtn.get_position(), phtn.get_angle(), surface_cross);
    dist_to_census = phtn.get_distance_remaining();

    // select minimum distance event
    dist_to_event = min(dist_to_scatter, min(dist_to_boundary, dist_to_census));

    // calculate energy absorbed by material, update photon and material energy
    ew_factor = exp(-sigma_a * f * dist_to_event);
    absorbed_E = phtn.get_E() * (1.0 - ew_factor);

    rank_track_E[cell_id] += absorbed_E / (sigma_a * f);
    rank_abs_E[cell_id] += absorbed_E;

    phtn.set_E(phtn.get_E() - absorbed_E);

    // update position
    phtn.move(dist_to_event);

    // apply variance/runtime reduction
    if (phtn.below_cutoff(cutoff_fraction)) {
      rank_abs_E[cell_id] += phtn.get_E();
      active = false;
      event = KILL;
    }
    // or apply event
    else {
      // EVENT TYPE: SCATTER
      if (dist_to_event == dist_to_scatter) {
        get_uniform_angle(angle, rng);
        phtn.set_angle(angle);
        if (rng->generate_random_number() >
            (sigma_s / ((1.0 - f) * sigma_a + sigma_s)))
          phtn.set_group(sample_emission_group(rng, cell));
      }
      // EVENT TYPE: BOUNDARY CROSS
      else if (dist_to_event == dist_to_boundary) {
        boundary_event = cell.get_bc(surface_cross);
        if (boundary_event == ELEMENT) {
          next_cell = cell.get_next_cell(surface_cross);
          phtn.set_cell(next_cell);
          cell_id = next_cell;
          cell = mesh.get_on_rank_cell(cell_id);
        } else if (boundary_event == VACUUM) {
          exit_E += phtn.get_E();
          active = false;
          event = EXIT;
        } else
          phtn.reflect(surface_cross);
      }
      // EVENT TYPE: REACH CENSUS
      else if (dist_to_event == dist_to_census) {
        phtn.set_distance_to_census(c * next_dt);
        active = false;
        event = CENSUS;
        census_E += phtn.get_E();
      }
    } // end event loop
  }   // end while alive
  return event;
}

std::vector<Photon> response_transport(Source &source, const Mesh &mesh,
                                         IMC_State &imc_state,
                                         std::vector<double> &rank_abs_E,
                                         std::vector<double> &rank_track_E,
					 Tally* tally, Sphere_Response* resp) {
  using Constants::CENSUS;
  using Constants::event_type;
  using Constants::EXIT;
  using Constants::KILL;
  using Constants::WAIT;
  using std::cout;
  using std::endl;
  using std::vector;

  double census_E = 0.0;
  double exit_E = 0.0;
  double next_dt = imc_state.get_next_dt(); //! Set for census photons
  double dt = imc_state.get_next_dt();      //! For making current photons

  RNG *rng = imc_state.get_rng();


  // timing
  Timer t_transport;
  t_transport.start_timer("timestep transport");

  // replicated transport does not require the global photon count
  uint64_t n_local = source.get_n_photon();

  //------------------------------------------------------------------------//
  // main transport loop
  //------------------------------------------------------------------------//

  vector<Photon> census_list;   //! End of timestep census list
  uint64_t n_local_sourced = 0; //! Photons pulled from source object
  Photon phtn;
  event_type event;

  uint32_t count = 0;

  //------------------------------------------------------------------------//
  // Transport photons from source
  //------------------------------------------------------------------------//
  while (n_local_sourced < n_local) {

//   std::cout << "\tTRACKING NEW PARTICLE" << std::endl;

    count++;

    if((count % 1000) == 0) 
    	std::cout << "\t\tCOUNT: " << count << std::endl;

    phtn = source.get_photon(rng, dt);
    n_local_sourced++;

    event = resp_transport_photon(phtn, mesh, rng, next_dt, exit_E,
                             census_E, rank_abs_E, rank_track_E, tally, resp);

    switch (event) {
    // this case should never be reached
    case WAIT:
      break;
    case KILL:
      break;
    case EXIT:
      break;
    case CENSUS:
      census_list.push_back(phtn);
      break;
    }
	
  } // end while

  
  cout << "\t# OF TRACKED PARTICLES: " << count << endl;


  // record time of transport work for this rank
  t_transport.stop_timer("timestep transport");

  // wait for all ranks to finish
  MPI_Barrier(MPI_COMM_WORLD);

  std::sort(census_list.begin(), census_list.end());

  // set diagnostic quantities
  imc_state.set_exit_E(exit_E);
  imc_state.set_post_census_E(census_E);
  imc_state.set_census_size(census_list.size());
  imc_state.set_rank_transport_runtime(
      t_transport.get_time("timestep transport"));

  return census_list;
}

#endif // def transport_replicated_h_
//---------------------------------------------------------------------------//
// end of transport_replicated.h
//---------------------------------------------------------------------------//
