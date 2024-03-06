#pragma once
 /**
  * @brief Matrix multiplication of doubles
  * @details Notation is chosen to reflect a classic BLAS DGEMM interface, although transposing-descriptive keywords are neglected.
  * 
  * A has m rows and k cols (m x k); B has k rows and n cols (k x n); C then has m rows and n cols (m x n).
  * 
  * The equation computed is
  * @f[
  * \mathbf{C} = \alpha \mathbf{A} \mathbf{B} + \beta \mathbf{C}
  * @f]
  * 
  * @param[in] m Number of rows of input matrix @p A (and, since, you know, matrix multiplication, also the number of rows of output matrix @p C)
  * @param[in] n Number of cols of input matrix @p B (and, matrix multiplication, number of cols of output matrix @p C)
  * @param[in] k Number of cols of input matrix @p A and, also, number of rows of input matrix @p B
  * @param[in] alpha Scaling factor applied to the multiplication of @p A and @p B (@f$\alpha \mathbf{A} \mathbf{B}@f$)
  * @param[in] A First input matrix of dimension @p m x @p lda
  * @param[in] lda Leading dimension of A, should be @p k
  * @param[in] B Second input matrix of dimension @p k x @p ldb
  * @param[in] ldb Leading dimension of B, should be @p n
  * @param[in] beta Scaling factor applied to old values of @p C
  * @param[in,out] C Output matrix of dimension @p m x @p ldc; the result of the multiplication @p alpha * @p A * @p B is added to the contents C (with an additional factor @p beta multiplied to the old @p C values)
  * @param[in] ldc Leading dimension of C, should be @p n
 **/
void DGEMM(const int m, const int n, const int k, 
			const double alpha,
			const double * __restrict__ A,
			const int lda,
			const double * __restrict__ B, 
			const int ldb,
			const double beta, 
			double * __restrict__ C,
			const int ldc); 

/**
 * OMP implementation of DGEMM.
 *   
 *  
 **/
void DGEMM_omp(const int m, const int n, const int k, 
      const double alpha,
      const double * __restrict__ A,
      const int lda,
      const double * __restrict__ B, 
      const int ldb,
      const double beta, 
      double * __restrict__ C,
      const int ldc); 

void DGEMM_loop(const int m, const int n, const int k, 
      const double alpha,
      const double * __restrict__ A,
      const int lda,
      const double * __restrict__ B, 
      const int ldb,
      const double beta, 
      double * __restrict__ C,
      const int ldc);

void DGEMM_simd(const int m, const int n, const int k, 
      const double alpha,
      const double * __restrict__ A,
      const int lda,
      const double * __restrict__ B, 
      const int ldb,
      const double beta, 
      double * __restrict__ C,
      const int ldc);
 
void DGEMM_asimd(const int m, const int n, const int k, 
      const double alpha,
      const double * __restrict__ A,
      const int lda,
      const double * __restrict__ B, 
      const int ldb,
      const double beta, 
      double * __restrict__ C,
      const int ldc); 

void DGEMM_vec(const int m, const int n, const int k, 
      const double alpha,
      const double * __restrict__ A,
      const int lda,
      const double * __restrict__ B, 
      const int ldb,
      const double beta, 
      double * __restrict__ C,
      const int ldc); 

void DGEMM_block(const int m, const int n, const int k, 
        const double alpha,
        const double * __restrict__ A,
        const int lda,
        const double * __restrict__ B, 
        const int ldb,
        const double beta, 
        double * __restrict__ C,
        const int ldc);
