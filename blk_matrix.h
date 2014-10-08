#ifndef _BLK_MATRIX_H
#define _BLK_MATRIX_H

#include "utils.h"
#include "csr_matrix.h"

// define the row and column's size, must smaller than 65536
#define BLOCK_SIZE (32 * 1024)

// in a block, if non_zeros number bigger than BLOCK_SIZE * CSR_THRESHOLD, then use compressed CSR format, else use compressed COO format
#define CSR_THRESHOLD 4 

// define the type of basic struct of matrix storage
typedef enum 
{
    BLK_CSR,  // block stored in compressed csr format
    BLK_COO,  // block stored in compressed coo format
} blk_type_t;

// matrix storage in blocked-mix format
struct blk_mat_t
{
    int rows;
    int cols;
    INT64 non_zeros;

    blk_type_t *types;

    DWORD *row_id;
    WORD *row_info;
    WORD *col_idx;
    FLOAT *vals;
};

// blocked-mix matrix container
struct blk_cont_t
{
    split_dir_t dir;
    int count;
    int *split_idx;

    struct blk_mat_t *blks;
};

// release blocked-mix matrix
void release_blk_mat(struct blk_mat_t *mat);

// release blocked-mix matrix container
void release_blk_cont(struct blk_cont_t *blk_cont);

// transform csr matrix format to blocked-mix matrix format
int csr_to_blk(struct csr_mat_t *csr, struct blk_mat_t *blk);

// transform csr matrix container to blocked-mix matrix container
int csr_cont_to_blk_cont(struct csr_cont_t *csr_cont, struct blk_cont_t *blk_cont);

// bind the memory of blocked-mix matrix container to its local node
int blk_cont_numa_bind(struct blk_cont_t *blk_cont);

#endif

