#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>

/**
 * @brief Matrix datatype, simplifying handling of matrices in C++ code
 * 
 */
class matrix_t {

/* TODO:(DONE) Implement
	* Constructor for any matrix
	* Constructor for square matrix
	* Destructor
	* Private data members for width, height
  Do not forget to allocate memory of the data array!
 */

public:
    matrix_t(int m,int n);  // constructor
    matrix_t(int n);
    ~matrix_t();                     // destructor
	friend std::ostream & operator<< (std::ostream & stream, 
            const matrix_t & matrix); // prints a matrix  

	double * data;

private: 
    int width;      // columns
    int height;     // rows
};
