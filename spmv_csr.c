#include "spmv_csr.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int spmv_csr(struct csr_mat_t *mat, FLOAT *x, FLOAT *y)
{
    DWORD *row_ptr = mat->row_ptr;
    int *col_idx = mat->col_idx;
    FLOAT *vals  = mat->vals;

    int rows = mat->rows;

    int i, j = 0;
    FLOAT tmp;
    for (i = 0; i < rows; i++)
    {
        tmp = 0.0;
        for (j = row_ptr[i]; j < row_ptr[i + 1]; j++)
        {
            tmp += vals[j] * x[col_idx[j]];
        }
        y[i] = tmp;
    }

    return 0;
}

int spmv_csrs(struct csr_cont_t *csr_cont, FLOAT *x, FLOAT *y, FLOAT **local_y)
{
    if (csr_cont->dir == SPLIT_HORIZON)
    {
        int i;
#pragma omp parallel for
        for (i = 0; i < csr_cont->count; i++)
        {
            spmv_csr(csr_cont->csrs + i, x, y + csr_cont->split_idx[i]);
        }
    }
    else if (csr_cont->dir == SPLIT_VERTICAL)
    {
        if (local_y == NULL)
        {
            fprintf(stderr, "local_y isn't provided in SPLIT_VERTICAL mode.\n");
            exit(0);
        }
        int i, j;
        int rows = csr_cont->csrs[0].rows;
        int count = csr_cont->count;
#pragma omp parallel for
        for (i = 0; i < count; i++)
        {
            spmv_csr(csr_cont->csrs + i, x + csr_cont->split_idx[i], local_y[i]);
        }
        FLOAT tmp;
#pragma omp parallel for private(i, j, tmp)
        for (i = 0; i < rows; i++)
        {
            tmp = 0.0;
            for (j = 0; j < count; j++)
            {
                tmp += local_y[j][i];
            }
            y[i] = tmp;
        }
    }
    return 0;
}

