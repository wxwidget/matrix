CC=g++
CCFLAG=-g 
INCLUDE=-I../

all:lu_main matrix_main qr_main svd_main vector_main
	#cholesky_main matrix_main qr_main svd_main vector_main
clean:
	rm lu_main matrix_main qr_main svd_main vector_main
lu_main:lu_main.cpp
	$(CC) $(INCLUDE) $(CCFLAG) -o $@ $^

matrix_main:matrix_main.cpp
	$(CC) $(INCLUDE) $(CCFLAG) -o $@  $^
qr_main:qr_main.cpp
	$(CC) $(INCLUDE) $(CCFLAG) -o $@  $^
svd_main:svd_main.cpp
	$(CC) $(INCLUDE) $(CCFLAG) -o $@  $^
vector_main:vector_main.cpp
	$(CC) $(INCLUDE) $(CCFLAG) -o $@  $^

