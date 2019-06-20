/*
 * @file:    response_funct.h
 * @author:  Scott Campbell
 * @version: 18 June 2019
 * @desc:    This class stores/holds the response value for a each
 *           cell in a grid system.
 */

#ifndef response_funct_h
#define response_funct_h

using namespace std;

class Response_Funct {

// TO DO:
//     * Add a reference to a tally surface in const
//     * Add a way to do ray tracing towards cells, 
//       based on position on the tally line and angle

public:

    Response_Funct(Tally& tally, uint32_t n_rays, uint32_t n_tally_points) : 
		   tally(tally), n_rays(n_rays), n_tally_points(n_tally_points) {
	mesh_start_x = tally.get_mesh_start_x();
	mesh_end_x = tally.get_mesh_end_x();
	mesh_start_y = tally.get_mesh_start_y();
	mesh_end_y = tally.get_mesh_end_y();
	
	tally_x1 = tally.get_x1();
	tally_y1 = tally.get_y1();
	tally_x2 = tally.get_x2();
	tally_y2 = tally.get_y2();
    }   

    ~Response_Funct() {}

private:

    Tally& tally;

    uint32_t n_rays; 
    uint32_t n_tally_points;

    double mesh_start_x, mesh_start_y;
    double mesh_end_x, mesh_end_y;
    double tally_x1, tally_y1;
    double tally_x2, tally_y2;

};
#endif
