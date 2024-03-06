#include <gtest/gtest.h>
#include "../src/gemm.hpp"

int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

//! Tests as usual,
TEST(dgemm, DGEMM_block) {  
    int m {50};
    int n {50};
    int k {50};
    double alpha {1.0};
    double beta {1.0};
    int lda {k};
    int ldb {n};
    int ldc {n};
    double *A = new double[n*n];
    double *B = new double[n*n];
    double *C = new double[n*n]; 
    double *correctC = new double[n*n];

    srand48(100);
    for (int i = 0; i < n*n; i++) {
        A[i] = drand48();
        B[i] = drand48();
        C[i] = drand48();
        correctC[i] = C[i];
    } 

    DGEMM_block(m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    DGEMM_simd(m, n, k, alpha, A,lda, B, ldb, beta, correctC, ldc);

    // Checks if values are the same.
    for (int i = 0; i < n * n; i++) {
        ASSERT_DOUBLE_EQ(correctC[i], C[i]);
    }
}

TEST(dgemm, DGEMM_omp) {
    int m {50};
    int n {50};
    int k {50};
    double alpha {1.0};
    double beta {1.0};
    int lda {k};
    int ldb {n};
    int ldc {n};
    double *A = new double[n*n];
    double *B = new double[n*n];
    double *C = new double[n*n]; 
    double *correctC = new double[n*n];

    srand48(100);
    for (int i = 0; i < n*n; i++) {
        A[i] = drand48();
        B[i] = drand48();
        C[i] = drand48();
        correctC[i] = C[i];
    } 

    DGEMM(m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    DGEMM_omp(m, n, k, alpha, A, lda, B, ldb, beta, correctC, ldc);

    // Checks if values are the same.
    for (int i = 0; i < n * n; i++) {
        ASSERT_DOUBLE_EQ(C[i], correctC[i]);
    }
}

TEST(dgemm, DGEMM_simd) {
    int m {50};
    int n {50};
    int k {50};
    double alpha {1.0};
    double beta {1.0};
    int lda {k};
    int ldb {n};
    int ldc {n};
    double *A = new double[n*n];
    double *B = new double[n*n];
    double *C = new double[n*n]; 
    double *correctC = new double[n*n];

    srand48(100);
    for (int i = 0; i < n*n; i++) {
        A[i] = drand48();
        B[i] = drand48();
        C[i] = drand48();
        correctC[i] = C[i];
    } 

    DGEMM(m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    DGEMM_simd(m, n, k, alpha, A, lda, B, ldb, beta, correctC, ldc);

    // Checks if values are the same.
    for (int i = 0; i < n * n; i++) {
        ASSERT_DOUBLE_EQ(C[i], correctC[i]);
    }
}
