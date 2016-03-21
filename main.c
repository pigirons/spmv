#include "utils.h"
#include "spmv_csr.h"
#include "spmv_blk.h"

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <pthread.h>
#include <time.h>
#include <sched.h>
#include <numa.h>
#include <numaif.h>

#define LOOP_TIME 1024

static double get_sec(struct timespec *start, struct timespec *end)
{
    return (double)(end->tv_sec - start->tv_sec + (end->tv_nsec - start->tv_nsec) * 1e-9);
}

static void thread_bind(int cpu)
{
    cpu_set_t cpu_set;
    CPU_ZERO(&cpu_set);
    CPU_SET(cpu, &cpu_set);
    if (pthread_setaffinity_np(pthread_self(),
            sizeof(cpu_set_t), &cpu_set) != 0)
    {
        fprintf(stderr, "Error: cpu[%d] bind failed.\n", cpu);
        exit(0);
    }
}

int result_file(FLOAT *y, int len)
{
    static int rec_id = 1;
    static char file_name[32];

    memcpy(file_name, "result_", 7);
    sprintf(file_name + 7, "%d.bin", rec_id);

    FILE *fp = fopen(file_name, "wb");
    fwrite(y, len, sizeof(FLOAT), fp);
    fclose(fp);

    rec_id++;

    return 0;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        fprintf(stderr, "usage: %s csr_matrix_file\n", argv[0]);
        exit(0);
    }

    int i, j, k;

    struct timespec start, end;

    int num_threads = 1;
#pragma omp parallel
    {
#pragma omp master
        {
            num_threads = omp_get_num_threads();
        }
    }
    printf("Thread number: %d.\n", num_threads);

#pragma omp parallel for
    for (i = 0; i < num_threads; i++)
    {
        int cpu = omp_get_thread_num();
        thread_bind(cpu);
    }

    FILE *fp;

    struct csr_mat_t csr, csr_re, csr_t, csr_t_re, csr_t_t;
    struct blk_mat_t blk;

    struct csr_cont_t csr_h, csr_v;
    struct blk_cont_t blk_h, blk_t_h;

    read_csr_mat(argv[1], &csr);
    int rows = csr.rows;
    int cols = csr.cols;
    INT64 non_zeros = csr.non_zeros;

    csr_transpose(&csr, &csr_t);
    release_csr_mat(&csr);
    int *reorder_map = (int*)malloc(cols * sizeof(int));
    csr_reorder(&csr_t, &csr_re, reorder_map);
    release_csr_mat(&csr_t);
    csr_transpose(&csr_re, &csr_t_t);
    release_csr_mat(&csr_re);
    split_csr_lb_nz(&csr_t_t, &csr_h, num_threads, SPLIT_HORIZON);
    release_csr_mat(&csr_t_t);
    csr_cont_to_blk_cont(&csr_h, &blk_h);
    release_csr_cont(&csr_h);

    printf("Notify: finished the preprocessing.\n");

    FLOAT *x = (FLOAT*)numa_alloc(cols * sizeof(FLOAT));
    FLOAT *y = (FLOAT*)numa_alloc(rows * sizeof(FLOAT));

    for (i = 0; i < cols; i++)
    {
        x[i] = 1.0;
    }

    // warm up
    spmv_blks(&blk_h, x, y, NULL);

    printf("Notify: begin csr spmv.\n");
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (i = 0; i < LOOP_TIME; i++)
    {
        spmv_blks(&blk_h, x, y, NULL);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    double time = get_sec(&start, &end) / LOOP_TIME;
    double gflops = 2.0 * non_zeros / time * 1e-9;
    printf("Notify: blk spmv time = %lfs, perf = %lf GFLOPS.\n", time, gflops);
    // result_file(y, rows);

    return 0;
}

