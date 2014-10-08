#include "blk_matrix.h"

#include <stdlib.h>
#include <stdio.h>
#include <sched.h>
#include <numa.h>
#include <numaif.h>

void release_blk_mat(struct blk_mat_t *mat)
{
    int blk_row = (mat->rows + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int blk_col = (mat->cols + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int blk_num = blk_row * blk_col;

    numa_free(mat->types, blk_num * sizeof(blk_type_t));
    numa_free(mat->row_info, mat->row_id[blk_num] * sizeof(WORD));
    numa_free(mat->col_idx, mat->non_zeros * sizeof(WORD));
    numa_free(mat->vals, mat->non_zeros * sizeof(FLOAT));
    numa_free(mat->row_id, (blk_num + 1) * sizeof(DWORD));
}

void release_blk_cont(struct blk_cont_t *blk_cont)
{
    int i;
    for (i = 0; i < blk_cont->count; i++)
    {
        release_blk_mat(blk_cont->blks + i);
    }
    free(blk_cont->split_idx);
}

int csr_to_blk(struct csr_mat_t *csr, struct blk_mat_t *blk)
{
    blk->rows = csr->rows;
    blk->cols = csr->cols;
    blk->non_zeros = csr->non_zeros;

    int blk_row = (csr->rows + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int blk_col = (csr->cols + BLOCK_SIZE - 1) / BLOCK_SIZE;

    printf("Notify: blk_row = %d, blk_col = %d.\n", blk_row, blk_col);

    int blk_num = blk_row * blk_col;

    blk->types = (blk_type_t*)numa_alloc_local(blk_num * sizeof(blk_type_t));
    blk->row_id = (DWORD*)numa_alloc_local((blk_num + 1) * sizeof(DWORD));
    blk->col_idx = (WORD*)numa_alloc_local(blk->non_zeros * sizeof(DWORD));
    blk->vals = (FLOAT*)numa_alloc_local(blk->non_zeros * sizeof(FLOAT));
    blk->row_id[0] = 0;

    int *blk_cnt = (int*)calloc(blk_num, sizeof(int));
    int idx;
    int i, j;
    int x, y;
    for (i = 0; i < csr->rows; i++)
    {
        x = i / BLOCK_SIZE;
        for (j = csr->row_ptr[i]; j < csr->row_ptr[i + 1]; j++)
        {
            y = csr->col_idx[j] / BLOCK_SIZE;
            idx = x * blk_col + y;
            blk_cnt[idx]++;
        }
    }

    int block_size;
    for (i = 0; i < blk_row; i++)
    {
        for (j = 0; j < blk_col; j++)
        {
            // get the real block size
            block_size = blk->rows - i * BLOCK_SIZE;
            if (block_size > BLOCK_SIZE)
            {
                block_size = BLOCK_SIZE;
            }

            idx = i * blk_col + j;

            if (blk_cnt[idx] >= block_size * CSR_THRESHOLD)
            {
                blk->row_id[idx + 1] = blk->row_id[idx] + block_size;
                blk->types[idx]= BLK_CSR;
            }
            else
            {
                blk->row_id[idx + 1] = blk->row_id[idx] + blk_cnt[idx];
                blk->types[idx]= BLK_COO;
            }
        }
    }
    blk->row_info = (WORD*)numa_alloc_local(blk->row_id[blk_num] * sizeof(WORD));
    memset(blk->row_info, 0, blk->row_id[blk_num] * sizeof(WORD));

    int *blk_idx = (int*)calloc(blk_num, sizeof(int));
    INT64 *blk_pos = (INT64*)malloc(blk_num * sizeof(INT64));
    blk_pos[0] = 0;
    for (i = 1; i < blk_num; i++)
    {
        blk_pos[i] = blk_pos[i - 1] + blk_cnt[i - 1];
    }

    WORD *cur_row_info;
    WORD *cur_col_idx;
    FLOAT *cur_vals;
    for (i = 0; i < csr->rows; i++)
    {
        x = i / BLOCK_SIZE;
        for (j = csr->row_ptr[i]; j < csr->row_ptr[i + 1]; j++)
        {
            y = csr->col_idx[j] / BLOCK_SIZE;
            idx = x * blk_col + y;

            cur_row_info = blk->row_info + blk->row_id[idx];
            cur_col_idx = blk->col_idx + blk_pos[idx];
            cur_vals = blk->vals + blk_pos[idx];

            if (blk->types[idx] == BLK_COO)
            {
                cur_row_info[blk_idx[idx]] = i - x * BLOCK_SIZE;
            }
            else if (blk->types[idx] == BLK_CSR)
            {
                cur_row_info[i - x * BLOCK_SIZE]++;
            }
            cur_col_idx[blk_idx[idx]] = csr->col_idx[j] - y * BLOCK_SIZE;
            cur_vals[blk_idx[idx]] = csr->vals[j];
            blk_idx[idx]++;
        }
    }

    free(blk_cnt);
    free(blk_idx);
    free(blk_pos);

    return 0;
}

int csr_cont_to_blk_cont(struct csr_cont_t *csr_cont, struct blk_cont_t *blk_cont)
{
    int i;
    blk_cont->dir = csr_cont->dir;
    blk_cont->count = csr_cont->count;
    blk_cont->split_idx = (int*)malloc((csr_cont->count + 1) * sizeof(int));
    memcpy(blk_cont->split_idx, csr_cont->split_idx, (csr_cont->count + 1) * sizeof(int));

    blk_cont->blks = (struct blk_mat_t*)malloc(blk_cont->count * sizeof(struct blk_mat_t));
#pragma omp parallel for private(i)
    for (i = 0; i < csr_cont->count; i++)
    {
        csr_to_blk(csr_cont->csrs + i, blk_cont->blks + i);
    }
    return 0;
}

