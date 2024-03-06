#pragma once
#include <mpi.h>
#include <string>
#include <iostream>
#include <cmath>

/**
 * @brief Holds info on MPI stuff (pretty much world size and world rank).
 */
struct mpi_t {
	int world_size = -1; 
    int world_rank = -1;

};

/**
 * @brief Simple struct to convert between a process' rank and the associated 
 *        (x,y) coordinates in the grid of tiles.
 *
 * @remark Notice that x (the horizontal coordinate) is given first in the 
 *         constructors, whereas in usual matrix notation the row index 
 *         usually comes first.
*/
struct rankCoordinates {

	rankCoordinates(int _rank, int _ntiles);
	rankCoordinates(int _x, int _y, int _ntiles);

	int rank;
    int ntiles;
	int x;
    int y;

	rankCoordinates left(int steps = 1);
	rankCoordinates right(int steps = 1);
	rankCoordinates top(int steps = 1);
	rankCoordinates bottom(int steps = 1);
};

/**
 * @brief Holds all matrix-related info.
 */
struct params_t {

	int matrixSize;
	int nTilesPerRow;
    int nTilesPerCol;
	int tileWidth;
};

/**
 * @brief Fills MPI info struct.
 */
void setupMpi(mpi_t &mpi);


/**
 * @brief Fills matrix info struct with content.
 *
 * @details Needs @mpi to figure out number of tiles.
 */
void getParams(params_t & params, mpi_t & mpi, int argc, char ** argv);

/**
 * @brief Fills content of matrix according to some @p func.
 */
template <typename T>
void fillMatrix(const int rows, const int cols, double * matrix, T func) {
	for (int row = 0; row < rows; row++) {
		for (int col = 0; col < cols; col++) {
            matrix[row * cols + col] = func(row, col); 
		}
	}
}
