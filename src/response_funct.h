/*
 * @file:    response_funct.h
 * @author:  Scott Campbell
 * @version: 18 June 2019
 * @desc:    This class stores/holds the response value for a each
 *           cell in a grid system.
 */

#ifndef response_funct_h
#deine response_funct_h

using namespace std;

public:

Response_Funct() {
    generate_blank_table();
}

Response_Funct(int n_rows, int n_cols) : num_rows(n_rows), num_cols(n_cols) {
    generate_blank_table();
}

inline double get_response(int row, int col) {
    if(row >= 0 && row < num_rows && col >= 0 && col < num_cols) {
        return response_table[row][col];
    }
    return -1;
}

inline bool set_response(int row, int col, double new_response) {
    if(row >= 0 && row < num_rows && col >= 0 && col < num_cols) {
        response_table[row][col] = new_response;
        return true;
    }
    return false;
}

void set_response_table(int n_rows, int n_cols, double** new_table) {
    num_rows = n_rows;
    num_cols = n_cols;
    response_table = new_table;
}

void generate_blank_table() {
    response_table = new double[num_rows];
    for(int row = 0; row < num_rows; row++) {
        response_table[row] = new double[num_cols];
        for(int col = 0; col < num_cols; col++) {
            response_table[row][col] = 1;
        }
    }
}

private:

int num_rows = 3;
int num_cols = 3;
double** response_table; // Holds the avg. response of each cell

#endif
