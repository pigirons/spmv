#include "spmv_blk.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
//#include <sched.h>

int spmv_blk(struct blk_mat_t *mat, FLOAT *x, FLOAT *y)
{
    int blk_row = (mat->rows + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int blk_col = (mat->cols + BLOCK_SIZE - 1) / BLOCK_SIZE;

    blk_type_t *types = mat->types;
    DWORD *row_id = mat->row_id;
    WORD *row_info = mat->row_info;

    memset(y, 0, mat->rows * sizeof(FLOAT));

    int i, j, k;
    FLOAT *x_base, *y_base = y, tmp;
    int idx, end, blk_id = 0;
    WORD *cur_row_info;
    WORD *cur_row_info_end = row_info;
    WORD *cur_col_idx = mat->col_idx;
    FLOAT *cur_vals = mat->vals;
    for (i = 0; i < blk_row; i++)
    {
        x_base = x;
        for (j = 0; j < blk_col; j++)
        {
            cur_row_info = cur_row_info_end;
            cur_row_info_end = row_info + row_id[blk_id + 1];

            if (types[blk_id] == BLK_CSR)
            {
                idx = 0;
                for (k = 0; cur_row_info + k < cur_row_info_end; k++)
                {
                    end = idx + cur_row_info[k];
                    tmp = 0.0;
                    for (; idx < end; idx++)
                    {
                        tmp += x_base[cur_col_idx[idx]] * cur_vals[idx];
                    }
                    y_base[k] += tmp;
                }
                cur_col_idx += idx;
                cur_vals += idx;
            }
            else
            {
                for (k = 0; cur_row_info + k < cur_row_info_end; k++)
                {
                    y_base[cur_row_info[k]] += x_base[cur_col_idx[k]] * cur_vals[k];
                }
                cur_col_idx += k;
                cur_vals += k;
            }
            x_base += BLOCK_SIZE;
            blk_id++;
        }
        y_base += BLOCK_SIZE;
    }

    return 0;
}

int spmv_blks(struct blk_cont_t *blk_cont, FLOAT *x, FLOAT *y, FLOAT **local_y)
{
    int i, j;
    FLOAT tmp;

    if (blk_cont->dir == SPLIT_HORIZON)
    {
        int *split_idx = blk_cont->split_idx;

#pragma omp parallel for
        for (i = 0; i < blk_cont->count; i++)
        {
            spmv_blk(blk_cont->blks + i, x, y + split_idx[i]);
        }
    }
    else if (blk_cont->dir == SPLIT_VERTICAL)
    {
        if (local_y == NULL)
        {
            fprintf(stderr, "Error: local_y is null while in SPLIT_VERTICAL mode.\n");
            exit(0);
        }
        int *split_idx = blk_cont->split_idx;

#pragma omp parallel for
        for (i = 0; i < blk_cont->count; i++)
        {
            spmv_blk(blk_cont->blks + i, x + split_idx[i], local_y[i]);
        }

        int rows = blk_cont->blks[0].rows;
#pragma omp parallel for private(i, j, tmp)
        for (i = 0; i < rows; i++)
        {
            tmp = 0.0;
            for (j = 0; j < blk_cont->count; j++)
            {
                tmp += local_y[j][i];
            }
            y[i] = tmp;
        }
    }

    return 0;
}

