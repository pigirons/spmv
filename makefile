CC = gcc
OBJECTS = main.o spmv_csr.o csr_matrix.o spmv_blk.o blk_matrix.o
CFLAGS = -O3 -c -lpthread -fopenmp
OFLAGS = -O3
#-Wl,--start-group ./libnuma.a -Wl,--end-group
LIBS = -lpthread -fopenmp -lnuma

spmv: $(OBJECTS)
	$(CC) $(LIBS) -o spmv $(OBJECTS) $(OFLAGS)

main.o: main.c
	$(CC) $(LIBS) $(CFLAGS) main.c

spmv_blk.o: spmv_blk.c spmv_blk.h
	$(CC) $(LIBS) $(CFLAGS) spmv_blk.c

blk_matrix.o: blk_matrix.c blk_matrix.h
	$(CC) $(LIBS) $(CFLAGS) blk_matrix.c

spmv_csr.o: spmv_csr.c spmv_csr.h
	$(CC) $(LIBS) $(CFLAGS) spmv_csr.c

csr_matrix.o: csr_matrix.c csr_matrix.h
	$(CC) $(LIBS) $(CFLAGS) csr_matrix.c

clean:
	rm $(OBJECTS) spmv 
