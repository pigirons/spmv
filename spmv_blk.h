#ifndef _SPMV_BLK_H
#define _SPMV_BLK_H

#include "blk_matrix.h"

// bind the memory of blocked vector
int blk_vec_numa_bind(struct blk_cont_t *blk_cont, FLOAT **vec, int vec_len);

// sparse matrix-vector multiplication in blocked-mix format, y = mat * x
int spmv_blk(struct blk_mat_t *mat, FLOAT *x, FLOAT *y);

// parallel sparse matrix-vector multiplication
int spmv_blks(struct blk_cont_t *blk_cont, FLOAT *x, FLOAT *y, FLOAT **local_y);

//sse multi-vector sparse matrix-vector multiplication in blocked-mix format, y = mat * x, x and y must be aligned 16 interleaved vectors
int spmv_blk_sse_4(struct blk_mat_t *mat, FLOAT *x, FLOAT *y);

// parallel sparse matrix-vector multiplication using sse
int spmv_blks_sse_4(struct blk_cont_t *blk_cont, FLOAT *x, FLOAT **y);

// sse multi-vector sparse matrix-vector multiplication in blocked-mix format, y = mat * x, x and y must be aligned 16 interleaved vectors
int spmv_blk_sse_8(struct blk_mat_t *mat, FLOAT *x, FLOAT *y);

// parallel sparse matrix-vector multiplication using sse
int spmv_blks_sse_8(struct blk_cont_t *blk_cont, FLOAT *x, FLOAT **y);

#endif

