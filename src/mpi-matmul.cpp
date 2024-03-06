#include <iostream>
#include <mpi.h>
#include <utility>
#include <cstdlib>
 #include <algorithm>
#include <iterator>
#include "mpi-helpers.hpp"
#include "matrix_t.hpp"
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
void communicate(params_t &params, matrix_t &toSend, matrix_t &toReceive, 
        int targetRank, int sourceRank, int tag) {
    MPI_Request request;
    MPI_Status status;

    MPI_Issend(toSend.data, params.tileWidth * params.tileWidth, MPI_DOUBLE, 
            targetRank, tag, MPI_COMM_WORLD, &request);

    MPI_Recv(toReceive.data, params.tileWidth * params.tileWidth, MPI_DOUBLE, 
            sourceRank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Wait(&request, &status);
}

int main(int argc, char **argv) {
    MPI_Init(& argc, & argv);

    mpi_t mpiInfo;
    setupMpi(mpiInfo);

    params_t params;
    getParams(params, mpiInfo, argc, argv);

    matrix_t subA{params.tileWidth};
    matrix_t subB{params.tileWidth};
    matrix_t subC{params.tileWidth};

    std::for_each(subA.data, subA.data+params.tileWidth*params.tileWidth,[](double& data) {
        data = 1;
        });
    std::for_each(subB.data, subB.data+params.tileWidth*params.tileWidth,[](double& data) {
        data = 2;
        });
    std::for_each(subC.data, subC.data+params.tileWidth*params.tileWidth,[](double& data) {
	data = 0;
        });

    double alpha = 1.0;
    double beta = 1.0;

    // Coordinates of each proc.
    rankCoordinates coords{mpiInfo.world_rank, params.nTilesPerCol}; 

    // Matrix to write into.
    matrix_t subAtemp{params.tileWidth};
    matrix_t subBtemp{params.tileWidth};

    // Step 1 -> skewing.
    if (coords.y != 0) { // shifting rows of A
        int shiftA = coords.y;
        communicate(params, subA, subAtemp, coords.left(shiftA).rank,
                coords.right(shiftA).rank, 0); 
        std::swap(subA.data, subAtemp.data);        
    }

    if (coords.x != 0) { // shifting cols of B
        int shiftB = coords.x;
        communicate(params, subB, subBtemp, coords.top(shiftB).rank,
                coords.bottom(shiftB).rank, 0); 
        std::swap(subB.data, subBtemp.data);
    }

    DGEMM(params.tileWidth, params.tileWidth, params.tileWidth, alpha,
            subA.data, params.tileWidth, subB.data, params.tileWidth, beta, 
            subC.data, params.tileWidth);

    // Step 2 -> iterative multiplication of tiles. 
    for (int i = 1; i < params.nTilesPerCol; i++) { 

        communicate(params, subA, subAtemp, coords.left().rank,
                coords.right().rank, 0); 
        communicate(params, subB, subBtemp, coords.top().rank,
                coords.bottom().rank, 0); 

        std::swap(subA.data, subAtemp.data);        
        std::swap(subB.data, subBtemp.data);

        DGEMM(params.tileWidth, params.tileWidth, params.tileWidth, alpha,
                subA.data, params.tileWidth, subB.data, params.tileWidth, beta,
                subC.data, params.tileWidth);
    }

    // Step 3 -> collect 'subC' in 'unsortedC'.
    //delete &subA;
    //delete &subB;
    matrix_t unsortedC{params.matrixSize};
    MPI_Gather(subC.data, params.tileWidth * params.tileWidth, MPI_DOUBLE,
            unsortedC.data, params.tileWidth * params.tileWidth, MPI_DOUBLE, 0,
            MPI_COMM_WORLD);
    // TODO: free up subC?
    //delete &subC;

    // Step 4 -> sort 'C'. 
    if (mpiInfo.world_rank == 0) {
        matrix_t C{params.matrixSize};

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
                    C.data[tileCindex + k / params.tileWidth *
                        params.matrixSize + k % params.tileWidth] =
                        unsortedC.data[tileUnsortedCindex + k];     
                }
            }
        }

        // Prints 'C' if not too large.
        if (params.matrixSize < 20) {
            std::cout << "C \n" << C << std::endl;  
        }    
    } 

    MPI_Finalize();
}
