#include "mpi-helpers.hpp"

void setupMpi(mpi_t &mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, & mpi.world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, & mpi.world_size);
}

void getParams(params_t & params, mpi_t & mpi, int argc, char ** argv) {
    params.matrixSize = 4; // default value

    if (argc > 1) {  
        std::string input(argv[1]); 
        params.matrixSize = std::stoi(input, nullptr); 
    }

    params.nTilesPerCol = std::lround(sqrt(mpi.world_size));
    params.nTilesPerRow = params.nTilesPerCol;
    params.tileWidth = params.matrixSize / params.nTilesPerCol;

    if (mpi.world_rank == 0 ) {
        if (mpi.world_size % params.nTilesPerCol != 0) {
            std::cout << "Wrong input! The square root of the number of"
                "processes must be a natural number. This insures that the"
                "matrix of processes is square."; 

            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (params.matrixSize % params.tileWidth != 0) {
            std::cout << "The square root of the number of processes must "
                "divide the size of the matrix without residual. This "
                "insures that each process stores equally big tiles of the "
                "original matrix.";

            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
}

rankCoordinates::rankCoordinates(int _rank, int _ntiles) {
    rank = _rank;
    ntiles = _ntiles;
    x = rank % _ntiles;  // row major
    y = rank / _ntiles;
};

rankCoordinates::rankCoordinates(int _x, int _y, int _ntiles) {
    x = _x;
    y = _y;
    ntiles = _ntiles;
    rank = _ntiles * y + x; // row major
};

rankCoordinates rankCoordinates::left(int steps) {
    int x_new = x - steps < 0 ? x - steps + ntiles : x - steps;
    return rankCoordinates(x_new, y, ntiles);
};

rankCoordinates rankCoordinates::right(int steps) {
    int x_new = (x + steps) % ntiles;
    return rankCoordinates(x_new, y, ntiles);
};

rankCoordinates rankCoordinates::top(int steps) {
    int y_new = y - steps < 0 ? y - steps + ntiles : y - steps;
    return rankCoordinates(x , y_new, ntiles);
};

rankCoordinates rankCoordinates::bottom(int steps) {
    int y_new = (y + steps) % ntiles;
    return rankCoordinates(x, y_new, ntiles);
};
