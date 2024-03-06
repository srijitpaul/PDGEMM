#include "gemm.hpp"
#include <stdexcept>
#include <omp.h>
#include <iostream>
#include <immintrin.h>
using std::min;

void DGEMM(const int m, const int n, const int k, 
        const double alpha,
        const double * __restrict__ A,
        const int lda ,
        const double * __restrict__ B, 
        const int ldb,
        const double beta, 
        double * __restrict__ C,
        const int ldc) {
    // NOTE: Ignore lda, ldb, ldc inside the function
    // TODO: Include OpenMP directive

    for (int i = 0; i < m; i++) {

        for(int j = 0; j < n; j++) {

            C[i * n + j] *= beta;			

            for (int l = 0; l < k; l++) {

                C[i * n + j] += A[i * k + l] * B[l * n + j] * alpha; 
                // TODO: (Andy) This flattened array indexing is "optimial" for
                // A: since the memory of A is read from one element to the next
                // adjacent element, in strides of 1. [l=0][l=1][l=2][l=3][...],
                // but not for B: since the memory of B is read from one 
                // element to the next, hopping to non-adjacent elements in 
                // strides of n. [l=0][--][--][l=1][--][--][l=2][...].
                // SOLUTION: perhaps we can flatten B, so columns are stacked
                // vertically, rather than in A, where rows are stacked 
                // horizontally. E.g. by transposing B, then indexing B[j * k + l]
            }

        }
    }
}

void DGEMM_omp(const int m, const int n, const int k, 
        const double alpha,
        const double * __restrict__ A,
        const int lda,
        const double * __restrict__ B, 
        const int ldb,
        const double beta, 
        double * __restrict__ C,
        const int ldc) {
    // NOTE: Ignore lda, ldb, ldc inside the function
    // TODO: Include OpenMP directive

    if (k != lda || n != ldb || n != ldc)
        throw std::invalid_argument("Invalid argument!"); 

#pragma omp parallel for 
    for (int i = 0; i < m; i++) {

        for(int j = 0; j < n; j++) {

            C[i * n + j] *= beta;			

            for (int l = 0; l < k; l++) {

                C[i * n + j] += A[i * k + l] * B[l * n + j] * alpha; 
            }

        }
    }
}


void DGEMM_loop(const int m, const int n, const int k, 
        const double alpha,
        const double * __restrict__ A,
        const int lda,
        const double * __restrict__ B, 
        const int ldb,
        const double beta, 
        double * __restrict__ C,
        const int ldc) {
    //int nthreads = omp_get_num_threads();
    ///printf("Number of threads = %d\n", nthreads);

    double dum; 
    if (k != lda || n != ldb || n != ldc)
        throw std::invalid_argument("Invalid argument!"); 

#pragma omp parallel for
    for (int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            C[i * n + j] *= beta;
        }

        //#pragma omp parallel for
        for(int l = 0; l < k; l++) {
            dum = A[i * k + l] * alpha;
            // #pragma omp  simd
            for (int j = 0; j < n; j++) {
                //  if(l == 0){C[i * n + j] *= beta;}
                C[i * n + j] += dum * B[l * n + j]; 
            }
            //C[i * n + l] *= beta;			
            //#pragma reduction(+:C[i * n + j])
        }
    }

}

void DGEMM_simd(const int m, const int n, const int k, 
        const double alpha,
        const double * __restrict__ A,
        const int lda,
        const double * __restrict__ B, 
        const int ldb,
        const double beta, 
        double * __restrict__ C,
        const int ldc) {
    // TODO: Include OpenMP directive

    //int nthreads = omp_get_num_threads();
    ///printf("Number of threads = %d\n", nthreads);

    double dum; 
    if (k != lda || n != ldb || n != ldc)
        throw std::invalid_argument("Invalid argument!"); 

#pragma omp parallel for
    for (int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            C[i * n + j] *= beta;
        }

        //#pragma omp parallel for
        for(int l = 0; l < k; l++) {
            dum = A[i * k + l] * alpha;
#pragma omp  simd
            for (int j = 0; j < n; j++) {
                //  if(l == 0){C[i * n + j] *= beta;}
                C[i * n + j] += dum * B[l * n + j]; 
            }
            //C[i * n + l] *= beta;			
            //#pragma reduction(+:C[i * n + j])
        }
    }
}

/*
   void DGEMM_vec(const int m, const int n, const int k, 
   const double alpha,
   const double * __restrict__ A,
   const int lda,
   const double * __restrict__ B, 
   const int ldb,
   const double beta, 
   double * __restrict__ C,
   const int ldc) {

//C[i * n + l] *= beta;                 

//#pragma omp parallel for
for(int l = 0; l < k; l++) {
//#pragma omp parallel for 
for (int j = 0; j < n;) {
int cow=n/8;
if (cow != 0){
//j += 8;
if(j==cow) {
C[i * n + j] += A[i * k + l] * B[l * n + j] * alpha; 
j++;
} else {
//C[i * n + j] += A[i * k + l] * B[l * n + j] * alpha; 

__m256 c = _mm256_load_ps(& C[i * n + j]);
__m256 a = _mm256_load_ps(& A[i * k + l]);
__m256 b = _mm256_load_ps(& B[l * n + j]);
__m256 ALPHA = _mm256_load_ps(& alpha);
register __m256 s0, s1;

s0 = _mm256_mul_ps(ALPHA, b);
s1 = _mm256_mul_ps(a, s0);
s1 = _mm256_add_ps(c, s1);

_mm256_store_ps(&C[i * n + j], s1);

j+=8;
}
} else {

C[i * n + j] += A[i * k + l] * B[l * n + j] * alpha; 
j++;
}
}
} //#pragma reduction(+:C[i * n + j])
}
*/


//! DEPRECATED
//void DGEMM_asimd(const int m, const int n, const int k, 
//        const double alpha,
//        const double * __restrict__ A,
//        const int lda,
//        const double * __restrict__ B, 
//        const int ldb,
//        const double beta, 
//        double * __restrict__ C,
//        const int ldc) {
//    // Aligns memory for simd
//    //int nthreads = omp_get_num_threads();
//    ///printf("Number of threads = %d\n", nthreads);
//    double dum; 
//    if (k != lda || n != ldb || n != ldc)
//        throw std::invalid_argument("Invalid argument!"); 
//
//#pragma omp parallel for 
//    for (int i = 0; i < m; i++) {
//        //#pragma omp simd aligned(C:256) 
//        //        for(int j = 0; j < n; j++) {
//        //            C[i * n + j] *= beta;
//        //        }
//        //#pragma omp parallel for
//        for(int l = 0; l < k; l++) {
//            //dum = A[i * k + l] * alpha;
//#pragma omp simd
//            for (int j = 0; j < n; j++) {
//                //  if(l == 0){C[i * n + j] *= beta;}
//                //C[i * n + j] += dum * B[l * n + j]; 
//                C[i * n + j] += A[i * k + l] * alpha * B[l * n + j]; 
//            }
//            //C[i * n + l] *= beta;			
//            //#pragma reduction(+:C[i * n + j])
//        }
//    }
//}

void DGEMM_block(const int m, const int n, const int k, 
        const double alpha,
        const double * __restrict__ A,
        const int lda,
        const double * __restrict__ B, 
        const int ldb,
        const double beta, 
        double * __restrict__ C,
        const int ldc) {

    const int b = 10; // TODO: set value for caching

    // If one of the matrix dimensions are less than 'b', 
    // calls an alternative method.
    if (min(b,n) != b || min(b,m) != b || min(b,k) != b ||
            n % b != 0 || m % b != 0 || k % b != 0) {
        DGEMM_simd(m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
        return;
    } 

    double dum; 
    #pragma omp parallel for
    for (int i = 0; i < m; i+=b) {
        for(int j = 0; j < n; j+=b) {
            for(int l = 0; l < k; l+=b) {

                for (int ii = 0; ii < b; ii++) {
                    for(int jj = 0; jj < b; jj++) {
                        C[(ii + i) * n + jj + j] *= beta;
                    }

                    for (int ll = 0; ll < b; ll++) {
                        dum = A[(ii + i) * k + ll + l] * alpha;
                        #pragma omp simd
                        for (int jj = 0; jj < b; jj++) {
                            C[(ii + i)* n + jj + j] += 
                                dum * B[(ll + l) * n + jj + j];       
                        } 
                    } 
                }
            }
        }
    }
}
