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

#include "../response_funct.h"
#include "../photon.h"
#include "../tally.h"
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
      string filename("simple_input.xml");
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
      string three_reg_filename("three_region_mesh_input.xml");
      Input three_reg_input(three_reg_filename, mpi_types);

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

        double prev_pos_1[3] = {1.2394,1.3421,1};
        double prev_pos_2[3] = {2.1245,1.6433,1};
        double prev_pos_3[3] = {1.2234,3.9863,1};

        double pos_1[3] = {3.3253,5.2523,1};       
        double pos_2[3] = {2.7356,5.3532,1};       
        double pos_3[3] = {3.8623,3.8325,1};
       
        double prev_pos_4[3] = {3.9323,5.0234,1};
        double prev_pos_5[3] = {2.4852,5.9212,1};
        double prev_pos_6[3] = {1.0112,0.5634,1};

        double pos_4[3] = {6.0293,6.7652,1}; 
        double pos_5[3] = {4.1724,5.4621,1};       
        double pos_6[3] = {1.2189,0.5126,1};       

	string filename("simple_input.xml");
	const Info mpi_info;
	MPI_Types mpi_types;
	Input input(filename, mpi_types);
	IMC_Parameters imc_p(input);
	
	Mesh mesh(input, mpi_types, mpi_info, imc_p); // Create a mesh
	Tally* tally = new Tally(1, 5, 3, 1, mesh); // Create a line tally

	ofstream outfile("RESULTS.txt");
	outfile << "Mesh start x: " << mesh.get_start_x() << endl;
	outfile << "Mesh start y: " << mesh.get_start_y() << endl;
	outfile << "Mesh start z: " << mesh.get_start_z() << endl;
	outfile << "Mesh end x: " << mesh.get_end_x() << endl;
	outfile << "Mesh end y: " << mesh.get_end_y() << endl;
	outfile << "Mesh end z: " << mesh.get_end_z() << endl;
	outfile << "Mesh dx: " << mesh.get_dx() << endl;
	outfile << "Mesh dy: " << mesh.get_dy() << endl;        
	outfile << "Mesh dz: " << mesh.get_dz() << endl;
	outfile.close();
   
        Photon phtn; 
	// Test the mesh get mesh dimension methods
	if(mesh.get_start_x() != 0) { tally_surface_pass = false; }
	if(mesh.get_end_x() != 10) { tally_surface_pass = false; }

	if(mesh.get_start_y() != 0) { tally_surface_pass = false; }
	if(mesh.get_end_y() != 40) { tally_surface_pass = false; }

	if(mesh.get_start_z() != 0) { tally_surface_pass = false; }
	if(mesh.get_end_z() != 90) { tally_surface_pass = false; }

	// Test the mesh get_dx ... functions
	if(mesh.get_dx() != 1) { tally_surface_pass = false; }
	if(mesh.get_dy() != 2) { tally_surface_pass = false; }
	if(mesh.get_dz() != 3) { tally_surface_pass = false; }

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

	// Test the motion filter of the tally
	tally->set_motion_filter(tally->POSITIVE_ONLY_FILTER);
	
	phtn.set_position(pos_1);
	phtn.set_prev_position(prev_pos_1);

	if(!tally->hit_tally(phtn)) { tally_surface_pass = false; }
	
	phtn.set_position(prev_pos_1);
	phtn.set_prev_position(pos_1);
	if(tally->hit_tally(phtn)) { tally_surface_pass = false; }

	tally->set_motion_filter(tally->NEGATIVE_ONLY_FILTER);
	
	phtn.set_position(pos_1);
	phtn.set_prev_position(prev_pos_1);
	if(tally->hit_tally(phtn)) { tally_surface_pass = false; }

	phtn.set_position(prev_pos_1);
	phtn.set_prev_position(pos_1);
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

       if (tally_surface_pass) {
	 cout << "TEST PASSED: tally surface count" << endl;
       } else {
	 cout << "TEST FAILED: tally surface count" << endl;
         nfail++;
       }
    }  
   
    // Test the response funtion class
    {
	bool response_pass = true;
	// Need to create a mesh
	string filename("simple_input.xml");
	const Info mpi_info;
	MPI_Types mpi_types;
	Input input(filename, mpi_types);
	IMC_Parameters imc_p(input);
	
	Mesh mesh(input, mpi_types, mpi_info, imc_p); // Create a mesh	
	Tally tally(1, 5, 3, 1, mesh); // Create a line tally
	Response_Funct resp(tally, mesh);
	
	//Define a box: lower_left = (0,0), and upper_right = (10,10)
	//v1 - v3 are false cases, v4 - v8 are true cases
	double v1 = resp.intersect_cell(resp.make_point(0,0), 
		         resp.make_point(10,10), resp.make_ray(1,11));
	double v2 = resp.intersect_cell(resp.make_point(0,0), 
		         resp.make_point(10,10), resp.make_ray(0,-1));
	double v3 = resp.intersect_cell(resp.make_point(0,0), 
		         resp.make_point(10,10), resp.make_ray(-1,21));

	double v4 = resp.intersect_cell(resp.make_point(0,0), 
		         resp.make_point(10,10), resp.make_ray(1,0));
	double v5 = resp.intersect_cell(resp.make_point(0,0), 
		         resp.make_point(10,10), resp.make_ray(-1,8));
	double v6 = resp.intersect_cell(resp.make_point(0,0), 
		         resp.make_point(10,10), resp.make_ray(0,5));
	double v7 = resp.intersect_cell(resp.make_point(0,0), 
		         resp.make_point(10,10), resp.make_ray(0,0));
	double v8 = resp.intersect_cell(resp.make_point(0,0), 
		         resp.make_point(10,10), resp.make_ray(1e8,0));


	ofstream outfile("RESULTS.txt");
	outfile << "v1 (-1): " << v1 << endl;	
	outfile << "v2 (-1): " << v2 << endl;	
	outfile << "v3 (-1): " << v3 << endl;	
	outfile << "v4 (14.142135): " << v4 << endl;	
	outfile << "v5 (11.3137): " << v5 << endl;
	outfile << "v6 (10): " << v6 << endl;
	outfile << "v7 (10): " << v7 << endl;
	outfile << "v8 (10): " << v8 << endl;
	outfile.close();
   

	if(v1 != -1) { response_pass = false; }
	if(v2 != -1) { response_pass = false; }
	if(v3 != -1) { response_pass = false; }
	
	if(v4 < 14.1420 || v4 > 14.1422) { response_pass = false; }
	if(v5 < 11.312 || v5 > 11.314) { response_pass = false; }
	if(v6 != 10) { response_pass = false; }
	if(v7 != 10) { response_pass = false; }
	if(v8 != 10) { response_pass = false; }

	if(response_pass) {
      	    cout << "TEST PASSED: Response Function" << endl;		
	} else {
	    cout << "TEST FAILED: Response Function" << endl;
	    nfail++;
	}
    }


  } // need to call destructors for mpi_types before MPI_Finalize

  MPI_Finalize();

  return nfail;
}
//----------------------------------------------------------------------------//
// end of test_mesh.cc
//----------------------------------------------------------------------------//
