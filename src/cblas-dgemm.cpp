#include "cblas-dgemm.hpp"
#include <stdexcept>
#include <iostream>
#include <fstream>
#include "mkl.h"


void CBLAS_DGEMM(const int m, const int n, const int k, 
			const double alpha,
			const double * __restrict__ A,
			const int lda,
			const double * __restrict__ B, 
			const int ldb,
			const double beta, 
			double * __restrict__ C,
			const int ldc) {
 
                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                m, n, k, alpha, A, k, B, n, beta, C, n);

}
