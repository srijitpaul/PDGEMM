#include <mpi.h>
#include "mpi-helpers.hpp"
#include "gemm.hpp"

#ifndef VERBOSE
#define VERBOSE 0  // Handy to create different levels of verbose output
#endif

// Can be used to print to cerr output
#define VERBOUT(LVL, MESSAGE) if (VERBOSE > LVL) { std::cerr << MESSAGE }  

/**
 * @brief Send/receive combination
 * 
 * @param toSend Matrix to send to @p targetRank
 * @param toReceive Matrix to receive from @p sourceRank
 * @param targetRank Rank number of to send @p toSend to
 * @param sourceRank Rank number to receive @p toReceive from
 * @param tag Some tag, can be ignored or used for the direction 
 *        or something similar
 */

int main(int argc, char **argv) {
    MPI_Init(& argc, & argv);

    mpi_t mpiInfo;
    setupMpi(mpiInfo);

    params_t params;
    getParams(params, mpiInfo, argc, argv);

    double alpha = 1.0;
    double beta = 1.0;

    double *subA = new double[params.tileWidth*params.tileWidth];
    double *subB = new double[params.tileWidth*params.tileWidth];
    double *subC = new double[params.tileWidth*params.tileWidth];
   
    for (int i = 0; i < params.tileWidth * params.tileWidth; i++) {
        subA[i] = 1;
        subB[i] = 2;
        subC[i] = 3;
    }
    
    // Matrix to write into.
    double *subAtemp = new double[params.tileWidth * params.tileWidth];
    double *subBtemp = new double[params.tileWidth * params.tileWidth];

    // Initialising communicator.
    int ndims = 2; // number of dimensions
    int *dims = new int[2] {params.nTilesPerRow, params.nTilesPerCol};
    int *periods = new int[2] {1 , 1};  // periodic in both cols and rows.
    int reorder = 1;
    MPI_Comm proc_grid;
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &proc_grid);
   
    // Initialising coordinates.
    int *coords = new int[ndims];
    MPI_Cart_coords(proc_grid, mpiInfo.world_rank, ndims, coords);
    int y = coords[1];
    int x = coords[0];

    // Destination and source ranks of process in communication.
    int source_rank;
    int target_rank;

    // Step 1 -> skewing.
    if (y != 0) { // shifting rows of A
        int shiftA = y;
        MPI_Cart_shift(proc_grid, 0, -shiftA, &source_rank, &target_rank); 
        MPI_Sendrecv(subA, params.tileWidth * params.tileWidth, MPI_DOUBLE, 
            target_rank, 0, subAtemp, params.tileWidth * params.tileWidth,
            MPI_DOUBLE, source_rank, 0, proc_grid, MPI_STATUS_IGNORE);
             
        std::swap(subA, subAtemp);        
    }

    if (x != 0) { // shifting cols of B
        int shiftB = x;
        MPI_Cart_shift(proc_grid, 1, -shiftB, &source_rank, &target_rank); 
        MPI_Sendrecv(subB, params.tileWidth * params.tileWidth, MPI_DOUBLE, 
            target_rank, 0, subBtemp, params.tileWidth * params.tileWidth,
            MPI_DOUBLE, source_rank, 0, proc_grid, MPI_STATUS_IGNORE);

        std::swap(subB, subBtemp);
    }

    DGEMM(params.tileWidth, params.tileWidth, params.tileWidth, alpha,
            subA, params.tileWidth, subB, params.tileWidth, beta, 
            subC, params.tileWidth);

    // Step 2 -> iterative multiplication of tiles. 
    for (int i = 1; i < params.nTilesPerCol; i++) { 
        MPI_Cart_shift(proc_grid, 0, -1, &source_rank, &target_rank);
        MPI_Sendrecv(subA, params.tileWidth * params.tileWidth, MPI_DOUBLE, 
            target_rank, 0, subAtemp, params.tileWidth * params.tileWidth,
            MPI_DOUBLE, source_rank, 0, proc_grid, MPI_STATUS_IGNORE);

        MPI_Cart_shift(proc_grid, 1, -1, &source_rank, &target_rank);
        MPI_Sendrecv(subB, params.tileWidth * params.tileWidth, MPI_DOUBLE, 
            target_rank, 0, subBtemp, params.tileWidth * params.tileWidth,
            MPI_DOUBLE, source_rank, 0, proc_grid, MPI_STATUS_IGNORE);

        std::swap(subA, subAtemp);        
        std::swap(subB, subBtemp);

        DGEMM(params.tileWidth, params.tileWidth, params.tileWidth, alpha,
                subA, params.tileWidth, subB, params.tileWidth, beta,
                subC, params.tileWidth);
    }

    // Step 3 -> collect 'subC' in 'unsortedC'.
    double *unsortedC = new double[params.matrixSize * params.matrixSize];
    MPI_Gather(subC, params.tileWidth * params.tileWidth, MPI_DOUBLE,
            unsortedC, params.tileWidth * params.tileWidth, MPI_DOUBLE, 0,
            MPI_COMM_WORLD);

    // Step 4 -> sort 'C'. 
    if (mpiInfo.world_rank == 0) {
        double *C = new double[params.matrixSize*params.matrixSize];

        int tileCindex;
        int tileUnsortedCindex;

        // Iterates over all tiles
        for (int i = 0; i < params.nTilesPerRow; i++) {
            for(int j = 0; j < params.nTilesPerCol; j++) {

                tileCindex = i * params.tileWidth * params.matrixSize + 
                    j * params.tileWidth; 
                tileUnsortedCindex = (i * params.nTilesPerRow + j) 
                    * params.tileWidth * params.tileWidth;

                // Iterates within tile     
                for (int k = 0; k < params.tileWidth * params.tileWidth; k++) {
                    C[tileCindex + k / params.tileWidth *
                        params.matrixSize + k % params.tileWidth] =
                        unsortedC[tileUnsortedCindex + k];     
                }
            }
        }

        // Prints 'C' if not too large.
        if (params.matrixSize < 20) {
            for (int i = 0; i < params.matrixSize; i++) {
                for (int j = 0; j < params.matrixSize; j++) {
                    std::cout << C[i * params.matrixSize + j] << "\t";
                }
                std::cout << std::endl;
            }
        }     
    } 

    MPI_Finalize();
}
