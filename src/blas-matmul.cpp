#include <iostream>
#include <stdexcept>
#include <cstdlib>
 #include <algorithm>
#include <iterator>
#include "helpers.hpp"
#include "matrix_t.hpp"
#include "cblas-dgemm.hpp"
#include <chrono> //@perf:n


#ifndef VERBOSE
	#define VERBOSE 0  // Handy to create different levels of verbose output, 
                       // decided at compile time.
#endif

// Can be used to print to cerr output
#define VERBOUT(LVL, MESSAGE) if (VERBOSE > LVL) { std::cerr << MESSAGE }  

int main(int argc, char * argv[]) {

    //@perf:n{ @STARTING-TIMER  ... initialization
        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    //@perf:n}

	params_t params; //@perf:o  this was just after main(){ ...

    try {
        getParams(params, argc, argv); 
    } catch (std::exception e) {
        printf("Invalid argument!\n");
        std::cerr << "Usage: " << argv[0] << " <matrixSize>" << std::endl;
        return 1;
    }

    matrix_t A {params.matrixSize};
    matrix_t B {params.matrixSize};
    matrix_t C {params.matrixSize};

    //@perf:n{ 
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        std::chrono::duration<double> time_span_t1t2 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << time_span_t1t2.count() << "s";
    //@perf:n} @ENDING-TIMER


    //@perf:n{ @STARTING-TIMER ... matrix filling
        std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
    //@perf:n}
	srand48(100);
	// TODO: Fill matrices with lambda function
        std::for_each(A.data, A.data+params.matrixSize*params.matrixSize,[](double& data) {
        data = drand48();
        });
        std::for_each(B.data, B.data+params.matrixSize*params.matrixSize,[](double& data) {
        data = drand48();
        });
        std::for_each(C.data, C.data+params.matrixSize*params.matrixSize,[](double& data) {
	data = 0;
        });
    

    //@perf:n{
        std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
        std::chrono::duration<double> time_span_t3t4 = std::chrono::duration_cast<std::chrono::duration<double>>(t4 - t3);
        std::cout << "\t" << time_span_t3t4.count() << "s";
    //@perf:n} @ENDING-TIMER

    // Printing matrices @perf:n
    // std::cout << "A\n" << A << std::endl; //@perf:n
    // std::cout << "B\n" << B << std::endl; //@perf:n
    // std::cout << "C\n" << C << std::endl; //@perf:n

    // TODO: What value to put here?
    double alpha = 1.0;
    double beta = 1.0;

    // matrix multiplication
    try {

        //@perf:n{ @STARTINNG-TIMER ... GEMM
            std::chrono::steady_clock::time_point t5 = std::chrono::steady_clock::now();
        //@perf:n} 
        
        CBLAS_DGEMM(params.matrixSize, params.matrixSize, params.matrixSize, alpha, 
              A.data, params.matrixSize, B.data, params.matrixSize, 
              beta, C.data, params.matrixSize);
        
        //@perf:n{
            std::chrono::steady_clock::time_point t6 = std::chrono::steady_clock::now();
            std::chrono::duration<double> time_span_t5t6 = std::chrono::duration_cast<std::chrono::duration<double>>(t6 - t5);
            std::cout << "\t" << time_span_t5t6.count() << "s";
        //@perf:n} @ENDING-TIMER
            std::cout << std::endl; // @END-ALL-TIMINGS

    } catch (std::exception e) {
        std::cout << e.what() << std::endl;
        return 1;
    }

    // std::cout << "C -> result \n" << C << std::endl;//@perf:n
}
