#pragma once // serves as inclusion guards

// ###############
// ### GENERAL
// ###############

/**
 * @brief Holds all matrix-related info
 */
struct params_t {
	int matrixSize = -1;
};

/**
 * @brief Fills matrix info struct with content
 */
void getParams(params_t & params, int argc, char ** argv); 

/**
 * @brief Fills content of matrix according to some @p func
 */
template <typename T>
void fillMatrix(const int rows, const int cols, double * matrix, T func) {
	for (int row = 0; row < rows; row++) {
		for (int col = 0; col < cols; col++) {
            // matrix[row * cols + col] = func(row, col); // TODO  //@perf:o
            matrix[row * cols + col] = func(); // TODO  //@perf:o
		}
	}
}
