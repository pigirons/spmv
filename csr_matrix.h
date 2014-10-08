#ifndef _CSR_MATRIX_H
#define _CSR_MATRIX_H

#include "utils.h"

// matrix storage in csr format
struct csr_mat_t
{
    int rows;
    int cols;
    INT64 non_zeros;

    DWORD *row_ptr;
    int *col_idx;
    FLOAT *vals;
};

// csr container split direction
typedef enum
{
    SPLIT_HORIZON,
    SPLIT_VERTICAL,
} split_dir_t;

// csr matrix container
struct csr_cont_t
{
    split_dir_t dir;
    int count;
    int *split_idx;

    struct csr_mat_t *csrs;
};

// release csr matrix
void release_csr_mat(struct csr_mat_t *mat);

// release the container of csr matrix
void release_csr_cont(struct csr_cont_t *csr_cont);

// reading data from local disk to csr format
int read_csr_mat(const char *file_name, struct csr_mat_t *mat);

// transpose the csr matrix
int csr_transpose(struct csr_mat_t *csr, struct csr_mat_t *csr_t);

// reorder csr matrix by rows
int csr_reorder(struct csr_mat_t *csr, struct csr_mat_t *csr_re, int *reorder_map);

// csr matrix split into count smaller csr matrices, satisfied load-balance by non_zeros
int split_csr_lb_nz(struct csr_mat_t *csr, struct csr_cont_t *csr_cont, int count, split_dir_t dir);

// bind the memory of csr matrix container to its local node
int csr_cont_numa_bind(struct csr_cont_t *csr_cont);

#endif

