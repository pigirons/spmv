#include "csr_matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <numa.h>
#include <numaif.h>

void release_csr_mat(struct csr_mat_t *mat)
{
    numa_free(mat->row_ptr, (mat->rows + 1) * sizeof(DWORD));
    numa_free(mat->col_idx, mat->non_zeros * sizeof(int));
    numa_free(mat->vals, mat->non_zeros * sizeof(FLOAT));
}

void release_csr_cont(struct csr_cont_t *csr_cont)
{
    int i;
    for (i = 0; i < csr_cont->count; i++)
    {
        release_csr_mat(csr_cont->csrs + i);
    }
    numa_free(csr_cont->split_idx, (csr_cont->count + 1) * sizeof(int));
    numa_free(csr_cont->csrs, csr_cont->count * sizeof(struct csr_mat_t));
}

int read_csr_mat(const char *file_name, struct csr_mat_t *mat)
{
    FILE *fp = fopen(file_name, "rb");
    if (fp == NULL)
    {
        return -1;
    }

    fread(&mat->rows, sizeof(int), 1, fp);
    fread(&mat->cols, sizeof(int), 1, fp);
    fread(&mat->non_zeros, sizeof(INT64), 1, fp);

    mat->row_ptr = (DWORD*)numa_alloc((mat->rows + 1) * sizeof(DWORD));
    mat->col_idx = (int*)numa_alloc(mat->non_zeros * sizeof(int));
    mat->vals    = (FLOAT*)numa_alloc(mat->non_zeros * sizeof(FLOAT));

    fread(mat->row_ptr, sizeof(DWORD), mat->rows + 1, fp);
    fread(mat->col_idx, sizeof(int), mat->non_zeros, fp);
    fread(mat->vals, sizeof(FLOAT), mat->non_zeros, fp);

    printf("Row x Column: %d x %d\n", mat->rows, mat->cols);
    printf("Non-zero elements number: %ld\n", mat->non_zeros);

    return 0;
}

int csr_transpose(struct csr_mat_t *csr, struct csr_mat_t *csr_t)
{
    csr_t->cols = csr->rows;
    csr_t->rows = csr->cols;
    csr_t->non_zeros = csr->non_zeros;

    csr_t->row_ptr = (DWORD*)numa_alloc((csr_t->rows + 1) * sizeof(DWORD));
    csr_t->col_idx = (int*)numa_alloc(csr_t->non_zeros * sizeof(int));
    csr_t->vals = (FLOAT*)numa_alloc(csr_t->non_zeros * sizeof(FLOAT));
    memset(csr_t->row_ptr, 0, (csr_t->rows + 1) * sizeof(DWORD));

    int i, j;
    for (i = 0; i < csr->rows; i++)
    {
        for (j = csr->row_ptr[i]; j < csr->row_ptr[i + 1]; j++)
        {
            csr_t->row_ptr[csr->col_idx[j] + 1]++;
        }
    }
    for (i = 1; i <= csr_t->rows; i++)
    {
        csr_t->row_ptr[i] += csr_t->row_ptr[i - 1];
    }

    int *row_start = (int*)malloc(csr_t->rows * sizeof(int));
    memcpy(row_start, csr_t->row_ptr, csr_t->rows * sizeof(int));

    for (i = 0; i < csr->rows; i++)
    {
        for (j = csr->row_ptr[i]; j < csr->row_ptr[i + 1]; j++)
        {
            int row = row_start[csr->col_idx[j]];
            csr_t->col_idx[row] = i;
            csr_t->vals[row] = csr->vals[j];
            row_start[csr->col_idx[j]]++;
        }
    }

    free(row_start);
    return 0;
}

// define the basic position operation used in heap sort
#define PARENT(A) (((A) - 1) >> 1)
#define LEFT(A) ((((A) + 1) << 1) - 1)
#define RIGHT(A) (((A) + 1) << 1)

static inline void swap_pos(int *array, int pos1, int pos2)
{
    int tmp     = array[pos1];
    array[pos1] = array[pos2];
    array[pos2] = tmp;
}

static void min_heapify(int *row_len, int *reorder_map, int rows, int pos)
{
    int min, tmp;
    int left, right;
    while (pos < rows)
    {
        // find minimum value from current position and its children
        min = pos;
        if ((left = LEFT(pos)) < rows && row_len[left] < row_len[pos])
        {
            min = left;
        }
        if ((right = RIGHT(pos)) < rows && row_len[right] < row_len[min])
        {
            min = right;
        }

        // if current position is the minimum
        if (pos == min)
        {
            break;
        }

        // swap row length
        swap_pos(row_len, pos, min);

        // swap reorder map
        swap_pos(reorder_map, pos, min);

        pos = min;
    }
}

static void row_sort(int *row_len, int *reorder_map, int rows)
{
    int i;

    // build the initial heap first
    for (i = (rows - 2) >> 1; i >= 0; i--)
    {
        min_heapify(row_len, reorder_map, rows, i);
    }

    // move heap top(minimum) to the end and reheapify
    for (i = rows - 1; i > 0; i--)
    {
        swap_pos(row_len, 0, i);
        swap_pos(reorder_map, 0, i);
        min_heapify(row_len, reorder_map, i, 0);
    }
}

int csr_reorder(struct csr_mat_t *csr, struct csr_mat_t *csr_re, int *reorder_map)
{
    int *row_len = (int*)malloc(csr->rows * sizeof(int));
    int i, j;
    for (i = 0; i < csr->rows; i++)
    {
        reorder_map[i] = i;
        row_len[i] = csr->row_ptr[i + 1] - csr->row_ptr[i];
    }
    row_sort(row_len, reorder_map, csr->rows);

    csr_re->rows = csr->rows;
    csr_re->cols = csr->cols;
    csr_re->non_zeros = csr->non_zeros;

    csr_re->row_ptr = (DWORD*)numa_alloc((csr_re->rows + 1) * sizeof(DWORD));
    csr_re->col_idx = (int*)numa_alloc(csr_re->non_zeros * sizeof(int));
    csr_re->vals = (FLOAT*)numa_alloc(csr_re->non_zeros * sizeof(FLOAT));

    int idx = 0;
    csr_re->row_ptr[0] = 0;
    for (i = 0; i < csr_re->rows; i++)
    {
        memcpy(csr_re->col_idx + idx, csr->col_idx + csr->row_ptr[reorder_map[i]], row_len[i] * sizeof(int));
        memcpy(csr_re->vals + idx, csr->vals + csr->row_ptr[reorder_map[i]], row_len[i] * sizeof(FLOAT));
        idx += row_len[i];
        csr_re->row_ptr[i + 1] = idx;
    }

    free(row_len);

    return 0;
}

int split_csr_lb_nz(struct csr_mat_t *csr, struct csr_cont_t *csr_cont, int count, split_dir_t dir)
{
    int i, j;

    csr_cont->dir = dir;
    csr_cont->count = count;
    csr_cont->split_idx = (int*)numa_alloc((count + 1) * sizeof(int));
    csr_cont->csrs = (struct csr_mat_t*)numa_alloc(count * sizeof(struct csr_mat_t));

    if (dir == SPLIT_HORIZON)
    {
        struct csr_mat_t *csrs = csr_cont->csrs;
        int *split_idx = csr_cont->split_idx;
        split_idx[0] = 0;

        int avg_ele = csr->non_zeros / count, split_val;
        for (i = 1, j = 1; i < count; i++)
        {
            split_val = i * avg_ele;
            while (csr->row_ptr[j] < split_val)
            {
                j++;
            }
            if (csr->row_ptr[j] - split_val > split_val - csr->row_ptr[j - 1])
            {
                j--;
            }
            split_idx[i] = j;
        }
        split_idx[i] = csr->rows;

        int item_idx = 0;
        for (i = 0; i < count; i++)
        {
            csrs[i].rows = split_idx[i + 1] - split_idx[i];
            printf("csrs[%d].rows = %d\n", i, csrs[i].rows);
            csrs[i].cols = csr->cols;
            csrs[i].non_zeros = csr->row_ptr[split_idx[i + 1]] - csr->row_ptr[split_idx[i]];

            csrs[i].row_ptr = (DWORD*)numa_alloc((csrs[i].rows + 1) * sizeof(DWORD));
            csrs[i].col_idx = (int*)numa_alloc(csrs[i].non_zeros * sizeof(int));
            csrs[i].vals = (FLOAT*)numa_alloc(csrs[i].non_zeros * sizeof(FLOAT));

            memcpy(csrs[i].row_ptr, csr->row_ptr + split_idx[i], (csrs[i].rows + 1) * sizeof(DWORD));
            memcpy(csrs[i].col_idx, csr->col_idx + item_idx, csrs[i].non_zeros * sizeof(int));
            memcpy(csrs[i].vals, csr->vals + item_idx, csrs[i].non_zeros * sizeof(FLOAT));

            for (j = 0; j <= csrs[i].rows; j++)
            {
                csrs[i].row_ptr[j] -= csr->row_ptr[split_idx[i]];
            }
            item_idx += csrs[i].non_zeros;
        }
    }
    else if (dir == SPLIT_VERTICAL)
    {
        struct csr_mat_t *csrs = csr_cont->csrs;
        int *split_idx = csr_cont->split_idx;
        split_idx[0] = 0;

        INT64 *col_cnt = (INT64*)calloc((csr->cols + 1), sizeof(INT64));
        for (i = 0; i < csr->rows; i++)
        {
            for (j = csr->row_ptr[i]; j < csr->row_ptr[i + 1]; j++)
            {
                col_cnt[csr->col_idx[j] + 1]++;
            }
        }

        int avg_ele = csr->non_zeros / count, split_val;
        int cur_col = 0;
        for (i = 1, j = 1; i < count; i++)
        {
            split_val = i * avg_ele;
            do 
            {
                cur_col += col_cnt[j++];
            }
            while (cur_col < split_val);
            if (cur_col - split_val > split_val - (cur_col - col_cnt[j]))
            {
                cur_col -= col_cnt[j--];
            }
            split_idx[i] = j;
        }
        split_idx[i] = csr->cols;

        for (i = 0; i < csr->cols; i++)
        {
            col_cnt[i + 1] += col_cnt[i];
        }

        for (i = 0; i < count; i++)
        {
            csrs[i].rows = csr->rows;
            csrs[i].cols = split_idx[i + 1] - split_idx[i];
            csrs[i].non_zeros = col_cnt[split_idx[i + 1]] - col_cnt[split_idx[i]];
            csrs[i].row_ptr = (DWORD*)numa_alloc((csrs[i].rows + 1) * sizeof(DWORD));
            csrs[i].col_idx = (int*)numa_alloc(csrs[i].non_zeros * sizeof(int));
            csrs[i].vals = (FLOAT*)numa_alloc(csrs[i].non_zeros * sizeof(FLOAT));
            memset(csrs[i].row_ptr, 0, (csrs[i].rows + 1) * sizeof(DWORD));
        }

        int col, k;
        for (i = 0; i < csr->rows; i++)
        {
            for (j = 0; j < count; j++)
            {
                csrs[j].row_ptr[i + 1] = csrs[j].row_ptr[i];
            }
            for (j = csr->row_ptr[i]; j < csr->row_ptr[i + 1]; j++)
            {
                col = csr->col_idx[j];
                for (k = 0; k < count; k++)
                {
                    if (col < split_idx[k + 1])
                    {
                        break;
                    }
                }
                csrs[k].col_idx[csrs[k].row_ptr[i + 1]] = csr->col_idx[j] - split_idx[k];
                csrs[k].vals[csrs[k].row_ptr[i + 1]] = csr->vals[j];
                csrs[k].row_ptr[i + 1]++;
            }
        }
        
        free(col_cnt);
    }

    return 0;
}

