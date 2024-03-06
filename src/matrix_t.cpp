#include <iostream>
#include <iomanip>
#include <fstream>
#include "matrix_t.hpp"

std::ostream & operator<< (std::ostream & stream, const matrix_t & matrix) {
	for (int row = 0; row < matrix.height; row++) {
		for (int col = 0; col < matrix.width; col++) {
			stream << std::fixed << std::setw(matrix.width) << std::setfill(' ') 
                << matrix.data[row * matrix.width + col] << " ";
		}
		stream << "\n";
	}
	stream << "\n";
	return stream;
}

matrix_t::matrix_t(int m, int n) {
    width = n;
    height = m;
    data = (double *) malloc(m * n * sizeof(double)); // C style declaration
}

matrix_t::matrix_t(int n) {
    width = n;
    height = n;
    data = (double *) malloc(n * n * sizeof(double)); 
}

matrix_t::~matrix_t() {
    free(data); // C style deallocation
}

