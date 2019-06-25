#ifndef response_h
#define response_h

#include <cmath>

using namespace std;


// Need a way to take a particles (x,y,z) and get the
// 	response cell it is in


struct Response_Cell {
    double* node_array; // holds the x_min, x_max ... 
    double total_d; // total distance of all particles
    double total_s_d; // total sigma_a * dist of all particles
    double sigma_a;
};

class Response {

public:
    Response(bool state) : n_x(10), n_y(10), n_z(10), dx(1), dy(1), dz(1) {
        
	// Dynamically allocate a 3D array to store the response cells
	resp_cells = new Response_Cell[n_x * n_y * n_z];
	for(int x_pos = 0; x_pos < n_x; x_pos++) {
	    for(int y_pos = 0; y_pos < n_y; y_pos++) {
		for(int z_pos = 0; z_pos < n_z; z_pos++) {
		    double pos[6] = {x_pos*dx,(x_pos+1)*dx,y_pos*dy,(y_pos+1)*dy,z_pos*dz,(z_pos+1)*dz};
		    Response_Cell rc = {pos, 0.0, 0.0, 0.0};
		    *(resp_cells + x_pos*n_y*n_z + y_pos*n_z + z_pos) = rc;
		}
	    }
	}
    }

    ~Response() {}

    Response_Cell& get_response_cell(double x, double y, double z) {
	int x_pos = (int)(x / dx);
	int y_pos = (int)(y / dy);
	int z_pos = (int)(z / dz);

	return *(resp_cells + x_pos*n_y*n_z + y_pos*n_z + z_pos);
    }
	
private:
 
Response_Cell* resp_cells;

int n_x, n_y, n_z;

double dx, dy, dz;


};

#endif



