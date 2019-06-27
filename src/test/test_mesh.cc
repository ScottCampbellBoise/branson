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

#include "../sphere_response.h"
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

	string filename("simple_input.xml");
	const Info mpi_info;
	MPI_Types mpi_types;
	Input input(filename, mpi_types);
	IMC_Parameters imc_p(input);
	
	Mesh mesh(input, mpi_types, mpi_info, imc_p); // Create a mesh
	Tally* tally = new Tally(5, 0, 0, 0, mesh); // Create a spherical tally
/*
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
*/ 
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

       if (tally_surface_pass) {
	 cout << "TEST PASSED: tally surface count" << endl;
       } else {
	 cout << "TEST FAILED: tally surface count" << endl;
         nfail++;
       }
    }  
  
 
    // Test the response funtion class
/*    {
	bool response_pass = true;
	// Need to create a mesh
	string filename("simple_input.xml");
	const Info mpi_info;
	MPI_Types mpi_types;
	Input input(filename, mpi_types);
	IMC_Parameters imc_p(input);
	
	Mesh mesh(input, mpi_types, mpi_info, imc_p); // Create a mesh	
	Tally tally(9, 20, 9, 10, mesh); // Create a line tally
	Response_Funct resp(tally, mesh);

	//Define a response_cell: lower_left = (0,0), and upper_right = (10,10)
	Point l_left = {0,0};
	Point u_right = {10,10};
	Response_Cell rc = {l_left, u_right, 0, 0, 0};
	
	// Define particles in that cell
	Particle p1 = { {5,5} , {1,0} };
	Particle p2 = { {3,0} , {0,0} };
	Particle p3 = { {5,10} , {1,5} };
	Particle p4 = { {10,13} , {0,12} };	

	// Test the get intersect point method
	Point int_pt = resp.get_intersect_point(rc, p1);
	if(int_pt.x != 0 && int_pt.y != 0) { response_pass = false; }
	int_pt = resp.get_intersect_point(rc, p2);
	if(int_pt.x != 0 && int_pt.y != 0) { response_pass = false; }	
	int_pt = resp.get_intersect_point(rc, p3);
	if(int_pt.x != 0 && int_pt.y != 5) { response_pass = false; }
	int_pt = resp.get_intersect_point(rc, p4);
	if(int_pt.x != NULL && int_pt.y != NULL) { response_pass = false; }

	// Test the move Particle method
	double table_dx = 10, table_dy = 10;
	Response_Cell** table = resp.get_empty_resp_table(5,5,table_dx,table_dy);

	Particle p5 = {{15,15} , {1,0}};
	if(!resp.move_particle(p5, table, table_dx, table_dy)) {response_pass = false;} 
	if(p5.pos.x != 10 && p5.pos.y != 10) { response_pass = false; }

	Particle p6 = {{5,5} , {1,0}};
	bool res6 = resp.move_particle(p6, table, table_dx, table_dy);
	if(!res6) { response_pass = false; } 	
	if(p6.pos.x != 0 && p6.pos.y != 0) { response_pass = false; }
	
	// Make a table of response_values
	resp.get_response(3,3);

	
	if(response_pass) {
      	    cout << "TEST PASSED: Response Function" << endl;		
	} else {
	    cout << "TEST FAILED: Response Function" << endl;
	    nfail++;
	}
   }
*/

  } // need to call destructors for mpi_types before MPI_Finalize

  MPI_Finalize();

  return nfail;
}
//----------------------------------------------------------------------------//
// end of test_mesh.cc
//----------------------------------------------------------------------------//
