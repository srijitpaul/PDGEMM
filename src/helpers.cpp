#include <iostream>
#include <stdexcept>
#include <string>
#include "helpers.hpp"

void getParams(params_t & params, int argc, char ** argv) {
	params.matrixSize = 4;  // default matrix size (if no command line param)
	if (argc == 2) {
        std::string input(argv[1]); 
        params.matrixSize = std::stoi(input, nullptr); 
	} else {
        throw std::invalid_argument("Invalid arguments!"); 
    }
}
