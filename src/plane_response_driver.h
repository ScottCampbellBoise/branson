#ifndef plane_response_driver_h_
#define plane_response_driver_h_

#include <functional>
#include <iostream>
#include <mpi.h>
#include <vector>
#include <fstream>

#include "census_creation.h"
#include "imc_parameters.h"
#include "imc_state.h"
#include "info.h"
#include "mesh.h"
#include "mpi_types.h"
#include "plane_response_transport.h"
#include "timer.h"
#include "write_silo.h"

#include "plane_tally.h"

void imc_plane_response_driver(Mesh &mesh, IMC_State &imc_state,
                         const IMC_Parameters &imc_parameters,
                         const MPI_Types &mpi_types, const Info &mpi_info,
                         Plane_Tally* tally, uint32_t n_resp_particles) {
    
    
    bool write_flux = false;
    bool write_cubanova = true;
    
    
    using std::vector;
    vector<double> abs_E(mesh.get_global_num_cells(), 0.0);
    vector<double> track_E(mesh.get_global_num_cells(), 0.0);
    vector<Photon> census_photons;
    int rank = mpi_info.get_rank();
    constexpr double fake_mpi_runtime = 0.0;
    
    uint64_t max_census_size = static_cast<uint64_t>(1.1*imc_parameters.get_n_user_photon());
    
    
    // Create an object to hold the sphere response class
    Plane_Response* resp = new Plane_Response(tally, mesh, imc_state, n_resp_particles);
    
    
 /*
    if (imc_parameters.get_write_silo_flag()) {
        // write SILO file
        write_silo(mesh, *resp, imc_state.get_time(), imc_state.get_step(),
                   imc_state.get_rank_transport_runtime(), fake_mpi_runtime, rank,
                   mpi_info.get_n_rank());
    }
*/
    
    double sourced_E = imc_state.get_pre_census_E();
    
    ofstream flux_file("flux_file.csv");
    flux_file << "Time,Reg. Fluence,Resp. Fluence,Reg. Flux,Resp. Flux" << endl;
    
    ofstream cuba_file("plane_cubanova_results.csv", ios::app);
    cuba_file << "time, reg flux, resp flux" << endl;
    
    double prev_resp = 0.0;
    double prev_reg = 0.0;
    
    while (!imc_state.finished()) {
        if (rank == 0)
            imc_state.print_timestep_header();
        
        // set opacity, Fleck factor, all energy to source
        mesh.calculate_photon_energy(imc_state);
        
        // all reduce to get total source energy to make correct number of
        // particles on each rank
        double global_source_energy = mesh.get_total_photon_E();
        MPI_Allreduce(MPI_IN_PLACE, &global_source_energy, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        
        imc_state.set_pre_census_E(get_photon_list_E(census_photons));
        
        // setup source
        Source source(mesh, imc_state, imc_parameters.get_n_user_photon(),
                      global_source_energy, census_photons);
        // no load balancing in particle passing method, just call the method
        // to get accurate count and map census to work correctly
        source.post_lb_prepare_source();
        
        imc_state.set_transported_particles(source.get_n_photon());
        
        // transport photons, return the particles that reached census
        census_photons = plane_response_transport(source, mesh, imc_state, abs_E,
                                            track_E, max_census_size, tally, resp, sourced_E);
        
        // reduce the abs_E and the track weighted energy (for T_r)
        MPI_Allreduce(MPI_IN_PLACE, &abs_E[0], mesh.get_global_num_cells(),
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &track_E[0], mesh.get_global_num_cells(),
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        mesh.update_temperature(abs_E, track_E, imc_state);
        
        // for replicated let root do conservation checks, zero out your tallies
        if (rank) {
            imc_state.set_absorbed_E(0.0);
            imc_state.set_pre_mat_E(0.0);
            imc_state.set_post_mat_E(0.0);
        }
        
        imc_state.print_conservation();
        
        //Print out tally information
        cout << "Regular Tally Info: " << endl;
        cout << "\tRegular Total Flux: \t\t" << tally->get_regular_E() << endl;
        cout << "\t# of crossings: \t" << tally->get_regular_hits() << endl;
        cout << "Response Tally Info: " << endl;
        cout << "\tResponse Total Flux: \t\t" << tally->get_response_E() << endl;
        cout << "\t# of crossings: \t" << tally->get_response_hits() << endl;
        
        if(write_cubanova) {
            cuba_file << imc_state.get_time() << "," << tally->get_regular_E() << "," << tally->get_response_E() << endl;
            tally->reset_regular_E();
            tally->reset_response_E();
        }
        
        if(write_flux) {
            double reg_flux = (tally->get_regular_E() - prev_reg) / imc_state.get_next_dt();
            double resp_flux = (tally->get_response_E() - prev_resp) / imc_state.get_next_dt();
            prev_reg = tally->get_regular_E();
            prev_resp = tally->get_response_E();
            
            
            flux_file << imc_state.get_time() << "," << tally->get_regular_E()
            << "," << tally->get_response_E() << "," << reg_flux << ","
            << resp_flux << endl;
        }
        
        tally->reset_regular_hits();
        tally->reset_response_hits();
/*
        if (imc_parameters.get_write_silo_flag() &&
            !(imc_state.get_step() % imc_parameters.get_output_frequency())) {
            // write SILO file
            write_silo(mesh, *resp, imc_state.get_time(), imc_state.get_step(),
                       imc_state.get_rank_transport_runtime(), fake_mpi_runtime, rank,
                       mpi_info.get_n_rank());
        }
*/
        // update time for next step
        imc_state.next_time_step();
    }
    
    if(write_flux)
        flux_file.close();
    if(write_cubanova) {
        cuba_file << "\n\n";
        cuba_file.close();
    }
}

#endif // response_driver_h_

//---------------------------------------------------------------------------//
// end of response_driver.h
//---------------------------------------------------------------------------//
