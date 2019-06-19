//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_cell.cc
 * \author Alex Long
 * \date   January 11 2016
 * \brief  Test cell check_in_cell and distance_to_boundary functions
 * \note   Copyright (C) 2018 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

#include "../cell.h"
#include "../constants.h"
#include "testing_functions.h"
#include <iostream>
#include <vector>

#include "../photon.h"
#include "../tally.h"
#include "../response_funct.h"

int main(void) {

  using Constants::bc_type;
  using Constants::dir_type;
  using std::cout;
  using std::endl;
  using std::vector;

  int nfail = 0;

  //setup cell
  Cell cell;
  // simple cube of size 1.0
  double x_low = 0.0;
  double x_high = 1.0;
  double y_low = 0.0;
  double y_high = 1.0;
  double z_low = 0.0;
  double z_high = 1.0;

  cell.set_coor(x_low, x_high, y_low, y_high, z_low, z_high);

  // test the check_in_cell function
  {
    double true_pos_1[3] = {0.5, 0.5, 0.5};
    double true_pos_2[3] = {0.9, 0.9, 0.9};
    double true_pos_3[3] = {0.1, 0.1, 0.1};
    double true_pos_4[3] = {0.1, 0.9, 0.1};
    double false_pos_1[3] = {-1.0, -1.0, -1.0};
    double false_pos_2[3] = {1.0, -1.0, -1.0};
    double false_pos_3[3] = {1.0, 1.0, -1.0};
    double false_pos_4[3] = {-1.0, 1.0, -1.0};

    bool check_in_cell_pass = true;
    //positions in cell
    if (!cell.check_in_cell(true_pos_1))
      check_in_cell_pass = false;
    if (!cell.check_in_cell(true_pos_2))
      check_in_cell_pass = false;
    if (!cell.check_in_cell(true_pos_3))
      check_in_cell_pass = false;
    if (!cell.check_in_cell(true_pos_4))
      check_in_cell_pass = false;
    //positions out of cell
    if (cell.check_in_cell(false_pos_1))
      check_in_cell_pass = false;
    if (cell.check_in_cell(false_pos_2))
      check_in_cell_pass = false;
    if (cell.check_in_cell(false_pos_3))
      check_in_cell_pass = false;
    if (cell.check_in_cell(false_pos_4))
      check_in_cell_pass = false;

    if (check_in_cell_pass)
      cout << "TEST PASSED: check_in_cell function" << endl;
    else {
      cout << "TEST FAILED: check_in_cell function" << endl;
      nfail++;
    }
  }

  // test distance to boundary function
  {
    double tolerance = 1.0e-8;

    double pos[3] = {0.5, 0.5, 0.5};

    double angle_1[3] = {0.999, 0.031614, 0.031614};
    double angle_2[3] = {0.031614, 0.999, 0.031614};
    double angle_3[3] = {0.031614, 0.31614, 0.999};

    double angle_4[3] = {-0.999, 0.031614, 0.031614};
    double angle_5[3] = {0.031614, -0.999, 0.031614};
    double angle_6[3] = {0.031614, 0.31614, -0.999};

    unsigned int surface_cross;
    bool distance_to_bound_pass = true;

    // tests for true
    double distance_1 =
        cell.get_distance_to_boundary(pos, angle_1, surface_cross);
    if (!soft_equiv(distance_1, 0.5 / 0.999, tolerance))
      distance_to_bound_pass = false;

    double distance_2 =
        cell.get_distance_to_boundary(pos, angle_2, surface_cross);
    if (!soft_equiv(distance_2, 0.5 / 0.999, tolerance))
      distance_to_bound_pass = false;

    double distance_3 =
        cell.get_distance_to_boundary(pos, angle_3, surface_cross);
    if (!soft_equiv(distance_3, 0.5 / 0.999, tolerance))
      distance_to_bound_pass = false;

    double distance_4 =
        cell.get_distance_to_boundary(pos, angle_4, surface_cross);
    if (!soft_equiv(distance_4, 0.5 / 0.999, tolerance))
      distance_to_bound_pass = false;

    double distance_5 =
        cell.get_distance_to_boundary(pos, angle_5, surface_cross);
    if (!soft_equiv(distance_5, 0.5 / 0.999, tolerance))
      distance_to_bound_pass = false;

    double distance_6 =
        cell.get_distance_to_boundary(pos, angle_6, surface_cross);
    if (!soft_equiv(distance_6, 0.5 / 0.999, tolerance))
      distance_to_bound_pass = false;

    if (distance_to_bound_pass)
      cout << "TEST PASSED: distance_to_boundary function" << endl;
    else {
      cout << "TEST FAILED: distance_to_boundary function" << endl;
      nfail++;
    }
  }

  // test get_volume function
  {
    bool get_volume_pass = true;

    double tolerance = 1.0e-8;

    if (!soft_equiv(cell.get_volume(), 1.0, tolerance))
      get_volume_pass = false;

    //make an oblong cell
    Cell oblong_cell;
    oblong_cell.set_coor(0.01, 0.02, 0.0, 10.0, -0.1, 0.1);

    if (!soft_equiv(oblong_cell.get_volume(), 0.02, tolerance))
      get_volume_pass = false;

    if (get_volume_pass)
      cout << "TEST PASSED: get_volume function" << endl;
    else {
      cout << "TEST FAILED: get_volume function" << endl;
      nfail++;
    }
  }

  // test construction from Proto_Cell
  {
    bool proto_construction_pass = true;

    const vector<dir_type> dirs{Constants::X_NEG, Constants::X_POS,
                                Constants::Y_NEG, Constants::Y_POS,
                                Constants::Z_NEG, Constants::Z_POS};

    //setup cell
    Proto_Cell proto_cell;

    uint32_t cell_ID = 3271733928; // 64-bit cell ID
    uint32_t region_ID = 12;
    uint32_t silo_index = 1231;

    vector<Constants::bc_type> bcs{Constants::REFLECT, Constants::VACUUM,
                                   Constants::ELEMENT, Constants::REFLECT,
                                   Constants::ELEMENT, Constants::VACUUM};
    vector<uint32_t> neighbors{3500000000, 3500000001, 3500000002,
                               3500000003, 3500000004, 3500000005};

    // simple cube of size 1.0
    double x_low = 0.0;
    double x_high = 1.0;
    double y_low = 0.0;
    double y_high = 1.0;
    double z_low = 0.0;
    double z_high = 1.0;
    vector<double> coords = {x_low, x_high, y_low, y_high, z_low, z_high};

    // set values
    proto_cell.set_coor(x_low, x_high, y_low, y_high, z_low, z_high);
    proto_cell.set_ID(cell_ID);
    proto_cell.set_grip_ID(cell_ID);
    proto_cell.set_region_ID(region_ID);
    proto_cell.set_silo_index(silo_index);
    for (auto i : dirs) {
      proto_cell.set_bc(i, bcs[i]);
      proto_cell.set_neighbor(i, neighbors[i]);
    }

    Cell from_proto = Cell(proto_cell);

    // test the get methods
    const double *cell_coords = from_proto.get_node_array();
    if (from_proto.get_ID() != cell_ID)
      proto_construction_pass = false;
    if (from_proto.get_grip_ID() != cell_ID)
      proto_construction_pass = false;
    if (from_proto.get_region_ID() != region_ID)
      proto_construction_pass = false;
    if (from_proto.get_silo_index() != silo_index)
      proto_construction_pass = false;
    for (int i = 0; i < 6; ++i) {
      if (from_proto.get_next_cell(i) != neighbors[i])
        proto_construction_pass = false;
      if (from_proto.get_bc(i) != bcs[i])
        proto_construction_pass = false;
      if (cell_coords[i] != coords[i])
        proto_construction_pass = false;
    }

    if (proto_construction_pass)
      cout << "TEST PASSED: construction from Proto_Cell" << endl;
    else {
      cout << "TEST FAILED: construction from proto cell" << endl;
      nfail++;
    }

    // test the tally surface and counter
    {
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

       Tally* tally = new Tally(1, 5, 3, 1); // Create a line tally
       Photon phtn;      

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


    // test response_funct.h class for proper getting and setting ...
    {
       bool resp_func_pass = true;

       Response_Funct* resp = new Response_Funct(); // generate a blank table
       
       // Check that the table is 'blank' - all 1's
       for(int r=0; r<3; r++) {
	 for(int c=0; c<3; c++) {
	   resp_func_pass = (resp->get_response(r,c) == 1);
	 }
       }

       if(resp_func_pass) {
	 cout << "TEST PASSED: response function get/set" << endl;
       } else {
	 cout << "TEST FAILED: response function get/set" << endl;
         nfail++;
       }
      
    }  

  }

  return nfail; // CHANGE TO nfail;
}
//---------------------------------------------------------------------------//
// end of test_cell.cc
//---------------------------------------------------------------------------//
