#ifndef _SPMV_CSR_H
#define _SPMV_CSR_H

#include "csr_matrix.h"

// bind the memory of blocked vector
int csr_vec_numa_bind(struct csr_cont_t *csr_cont, FLOAT **vec, int vec_len);

// sparse matrix-vector multiplication in csr format, y = mat * x
int spmv_csr(struct csr_mat_t *mat, FLOAT *x, FLOAT *y);

// parallel compute sparse matrix-vector multiplication in csr matrix
int spmv_csrs(struct csr_cont_t *csr_cont, FLOAT *x, FLOAT *y, FLOAT **local_y);

// sparse matrix-vector multiplication in csr format, y = mat * x, x and y are 4-interleaved vectors
int spmv_csr_sse_4(struct csr_mat_t *mat, FLOAT *x, FLOAT *y);

// parallel compute sparse matrix-vector multiplication in csr matrix
int spmv_csrs_sse_4(struct csr_cont_t *csr_cont, FLOAT *x, FLOAT **y);

// sparse matrix-vector multiplication in csr format, y = mat * x, x and y are 8-interleaved vectors
int spmv_csr_sse_8(struct csr_mat_t *mat, FLOAT *x, FLOAT *y);

// parallel compute sparse matrix-vector multiplication in csr matrix
int spmv_csrs_sse_8(struct csr_cont_t *csr_cont, FLOAT *x, FLOAT **y);

#endif

