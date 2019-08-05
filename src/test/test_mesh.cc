//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_mesh.h
 * \author Alex Long
 * \date   April 7 2016
 * \brief  Test region assignment after mesh construction
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include "../response_driver.h"
#include "../sphere_response.h"
#include "../tally.h"

#include "../plane_response_driver.h"
#include "../plane_response.h"
#include "../plane_tally.h"

#include "../replicated_driver.h"
#include "../constants.h"
#include "../input.h"
#include "../mpi_types.h"
#include "../imc_parameters.h"
#include "../proto_mesh.h"
#include "../mesh.h"
#include "testing_functions.h"

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  int nfail = 0;

  // scope for MPI type commit and frees (needs to be done before MPI_Finalize)
  {
    const Info mpi_info;

    MPI_Types mpi_types;	

    using std::cout;
    using std::endl;
    using std::string;

    // Test that a simple mesh (one division in each dimension) is constructed
    // correctly from the input file (simple_input.xml) and that each cell
    // is assigned the correct region
    {
      string filename("/users/campbell_s/branson/src/test/simple_input.xml");
      Input input(filename, mpi_types);
      Proto_Mesh mesh(input, mpi_types, mpi_info);

      bool simple_mesh_pass = true;

      uint32_t n_cell = mesh.get_n_local_cells();
      if (n_cell != 10 * 20 * 30)
        simple_mesh_pass = false;

      Proto_Cell cell;
      for (uint32_t i = 0; i < n_cell; i++) {
        cell = mesh.get_pre_window_allocation_cell(i);
        if (cell.get_region_ID() != 6)
          simple_mesh_pass = false;
      }

      if (simple_mesh_pass)
        cout << "TEST PASSED: simple mesh construction" << endl;
      else {
        cout << "TEST FAILED: simple mesh construction" << endl;
        nfail++;
      }
    }

    // Test that a multi-region mesh is constructed correctly from the input file
    // (three_region_input_mesh.xml) and that each cell is assigned the correct
    // region
    {
      bool three_region_mesh_pass = true;
      // first test large particle input file

      string filename("/users/campbell_s/branson/src/test/three_region_mesh_input.xml");
      Input three_reg_input(filename, mpi_types);

      Proto_Mesh mesh(three_reg_input, mpi_types, mpi_info);

      uint32_t n_cell = mesh.get_n_local_cells();
      if (n_cell != 21 * 10)
        three_region_mesh_pass = false;

      Proto_Cell cell;
      const double *coor;
      double x_low;
      // check the lower x position of the cell to see if the region matches
      // the divisions set in the input file
      for (uint32_t i = 0; i < n_cell; i++) {
        cell = mesh.get_pre_window_allocation_cell(i);
        coor = cell.get_node_array();
        x_low = coor[0];
        //cells in the first region
        if (x_low < 4.0) {
          if (cell.get_region_ID() != 230)
            three_region_mesh_pass = false;
        } else if (x_low >= 4.0 && x_low < 8.0) {
          if (cell.get_region_ID() != 177)
            three_region_mesh_pass = false;
        } else if (x_low >= 8.0) {
          if (cell.get_region_ID() != 11)
            three_region_mesh_pass = false;
        } else {
          // this should not occur, test fails
          three_region_mesh_pass = false;
        }
      }

      if (three_region_mesh_pass)
        cout << "TEST PASSED: three region mesh construction" << endl;
      else {
        cout << "TEST FAILED: three region mesh construction" << endl;
        nfail++;
      }
    }
   
    // test the tally surface and counter
    {
        using std::string;
        bool tally_surface_pass = true;

        double prev_pos_1[3] = {0,0,0};
        double prev_pos_2[3] = {2.12,1.6433,1};
        double prev_pos_3[3] = {1,3,1};

        double pos_1[3] = {5,5,6};       
        double pos_2[3] = {4.7,8.4,6.3};       
        double pos_3[3] = {7,7,7};
       
        double prev_pos_4[3] = {5.9323,5.0234,5};
        double prev_pos_5[3] = {5.4852,7.9212,7};
        double prev_pos_6[3] = {1.0112,0.5634,1};

        double pos_4[3] = {6.0293,6.7652,8}; 
        double pos_5[3] = {4.1724,5.4621,9};       
        double pos_6[3] = {1.2189,0.5126,1};       

        string filename("/users/campbell_s/branson/src/test/simple_input.xml");
	const Info mpi_info;
	MPI_Types mpi_types;
	Input input(filename, mpi_types);
	IMC_Parameters imc_p(input);
	
	Mesh mesh(input, mpi_types, mpi_info, imc_p); // Create a mesh
	Tally* tally = new Tally(5, 0, 0, 0, mesh); // Create a spherical tally

        Photon phtn; 
	
	// Test the basic tally operation
       	phtn.set_position(pos_1);
	phtn.set_prev_position(prev_pos_1);
       	if(!tally->hit_tally(phtn)) { tally_surface_pass = false; }
       	
	phtn.set_position(pos_2);
	phtn.set_prev_position(prev_pos_2);
       	if(!tally->hit_tally(phtn)) { tally_surface_pass = false; }
       	
	phtn.set_position(pos_3);
	phtn.set_prev_position(prev_pos_3);
	if(!tally->hit_tally(phtn)) { tally_surface_pass = false; }
	
	// Test the cases where the tally shouldn't be triggered
       	phtn.set_position(pos_4);
       	phtn.set_prev_position(prev_pos_4);
	if(tally->hit_tally(phtn)) { tally_surface_pass = false; }
       
	phtn.set_position(pos_5); 
	phtn.set_prev_position(prev_pos_5);
	if(tally->hit_tally(phtn)) { tally_surface_pass = false; } 
       	
	phtn.set_position(pos_6);	
	phtn.set_prev_position(prev_pos_6);
	if(tally->hit_tally(phtn)) { tally_surface_pass = false; }

 	// Check the dist_to_tally() method
 	double pos_7[3] = {0,0,0};
	double ang_7[3] = {1,0,0};
 	Photon phtn2;
	phtn2.set_position(pos_7);
	phtn2.set_angle(ang_7);
 	double dist = tally->get_dist_to_tally(phtn2);
	cout << "distance to tally (5): " << dist << endl;
	if(dist != 5) { tally_surface_pass = false; }

	double pos_8[3] = {2,0,0};
	double ang_8[3] = {1,0,0};
	phtn2.set_position(pos_8);
	phtn2.set_angle(ang_8);
 	dist = tally->get_dist_to_tally(phtn2);
	cout << "distance to tally (3): " << dist << endl;
	if(dist != 3) { tally_surface_pass = false; }

	double pos_9[3] = {6,6,6};
	double prev_pos_9[3] = {6,16,6};
	phtn2.set_position(pos_9);
	phtn2.set_prev_position(prev_pos_9);
 	dist = tally->get_dist_to_tally(phtn2);
	cout << "distance to tally (inf): " << dist << endl;
	if(!dist >= 1e6) { tally_surface_pass = false; }

        double pos_10[3] = {-15,0,0};
	double ang_10[3] = {1,0,0};
	phtn2.set_position(pos_10);
	phtn2.set_angle(ang_10);
 	dist = tally->get_dist_to_tally(phtn2);
	cout << "distance to tally (20): " << dist << endl;
	if(dist != 20) { tally_surface_pass = false; }


	if (tally_surface_pass) {
	  cout << "TEST PASSED: tally surface count" << endl;
        } else {
	  cout << "TEST FAILED: tally surface count" << endl;
          nfail++;
        }
    } 
/*
    {
	bool passed = true;

	string filename("/users/campbell_s/branson/src/test/point_source.xml");
	const Info mpi_info;
	MPI_Types mpi_types;
	Input input(filename, mpi_types);
	IMC_Parameters imc_p(input);	
        IMC_State imc_state(input, mpi_info.get_rank());

	Mesh mesh(input, mpi_types, mpi_info, imc_p); // Create a mesh
        mesh.initialize_physical_properties(input); // Initialize the physical props (T)

	Tally* tally = new Tally(.99, 1e-6, 1e-6, 1e-6, mesh); // Tally for point_source.xml
	
    	imc_response_driver(mesh, imc_state, imc_p, mpi_types, mpi_info, tally, 10000);

	// PRINT OUT THE TALLY INFORMATION
    	cout << "\n\tTally energy for Regular: \t" << tally->get_regular_E() << endl;
	cout << "\tTally energy for Response: \t" << tally->get_response_E() << endl << endl;
    }
*/
 
/*
    {
	bool passed = true;

	string filename("/users/campbell_s/branson/run/cubanova.xml");
	const Info mpi_info;
	MPI_Types mpi_types;
	Input input(filename, mpi_types);
	IMC_Parameters imc_p(input);	
        IMC_State imc_state(input, mpi_info.get_rank());

	Mesh mesh(input, mpi_types, mpi_info, imc_p); // Create a mesh
        mesh.initialize_physical_properties(input); // Initialize the physical props (T)

	Tally* tally = new Tally(1.99, 1e-6, 1e-6, 1e-6, mesh); // Tally for point_source.xml
	
    	imc_response_driver(mesh, imc_state, imc_p, mpi_types, mpi_info, tally, 10000);

	// PRINT OUT THE TALLY INFORMATION
    	cout << "\n\tTally energy for Regular: \t" << tally->get_regular_E() << endl;
	cout << "\tTally energy for Response: \t" << tally->get_response_E() << endl << endl;
    }
*/ 

/*
    // Testing the cubanova with SPHERICAL TALLY	  
    {
	int num_files = 10;	
	
	string files[] = {
	    "/users/campbell_s/branson/run/Nova_Files/cubanova1.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova2.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova3.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova4.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova5.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova6.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova7.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova8.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova9.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova10.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova11.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova12.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova13.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova14.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova15.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova16.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova17.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova18.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova19.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova20.xml", 
	};


	cout << "\n\nRUNNING CUBANOVA VARIANCE TEST WITH *** SPHERICAL TALLY *** ... \n\n";

	for(int k = 0; k < num_files; k++) {

 	   cout << "\nRunning File #" << (k+1) << "\n" << endl;

	   string filename = files[k]; 
	   const Info mpi_info;
	   MPI_Types mpi_types;
	   Input input(filename, mpi_types);
	   IMC_Parameters imc_p(input);	
           IMC_State imc_state(input, mpi_info.get_rank());

	   Mesh mesh(input, mpi_types, mpi_info, imc_p); // Create a mesh
           mesh.initialize_physical_properties(input); // Initialize the physical props (T)

	   Tally* tally = new Tally(1.99, 1e-6, 1e-6, 1e-6, mesh); // Tally for point_source.xml

    	   imc_response_driver(mesh, imc_state, imc_p, mpi_types, mpi_info, tally, 100000);
	}
	
	cout << "FINISHED RUNNING ALL TEST FILES" << endl;
    }
*/


    // Testing the cubanova with PLANAR TALLY	  
    {
	int num_files = 10;	
	
	string files[] = {
	    "/users/campbell_s/branson/run/Nova_Files/cubanova1.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova2.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova3.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova4.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova5.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova6.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova7.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova8.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova9.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova10.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova11.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova12.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova13.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova14.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova15.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova16.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova17.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova18.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova19.xml", 
	    "/users/campbell_s/branson/run/Nova_Files/cubanova20.xml", 
	};


	cout << "\n\nRUNNING CUBANOVA VARIANCE TEST WITH *** PLANAR TALLY *** ... \n\n";

	for(int k = 0; k < num_files; k++) {

 	   cout << "\nRunning File #" << (k+1) << "\n" << endl;

	   string filename = files[k]; 
	   const Info mpi_info;
	   MPI_Types mpi_types;
	   Input input(filename, mpi_types);
	   IMC_Parameters imc_p(input);	
           IMC_State imc_state(input, mpi_info.get_rank());

	   Mesh mesh(input, mpi_types, mpi_info, imc_p); // Create a mesh
           mesh.initialize_physical_properties(input); // Initialize the physical props (T)

	   double v1[3] = {1.9, 1.9, -1.9};
	   double v2[3] = {1.9, 1.9, 1.9};
	   double v3[3] = {-1.9, -1.9, -1.9};
	   double v4[3] = {-1.9, -1.9, 1.9};
	   Plane_Tally* tally = new Plane_Tally(v1, v2, v3, v4, mesh); // Tally for point_source.xml

    	   imc_plane_response_driver(mesh, imc_state, imc_p, mpi_types, mpi_info, tally, 100000);
	}
	
	cout << "FINISHED RUNNING ALL TEST FILES" << endl;
    }
	  
  }

  MPI_Finalize();

  return nfail;
}
//----------------------------------------------------------------------------//
// end of test_mesh.cc
//----------------------------------------------------------------------------//
