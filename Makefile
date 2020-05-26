SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -std=c++0x -O3 -DNDEBUG -pthread -march=native
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -std=c++0x -g2 -pthread
AUX_PAR_FLAGS = -fopenmp
#AUX_DISK_FLAGS = -DMONITOR_DISK_USAGE

all: compute_bwt_sequential compute_bwt_parallel

compute_bwt_sequential:
	$(CC) $(CFLAGS) -o compute_bwt_sequential src/pem_bwt_src/utils.cpp src/main.cpp $(AUX_DISK_FLAGS)

compute_bwt_parallel:
	$(CC) $(CFLAGS) $(AUX_PAR_FLAGS) -o compute_bwt_parallel src/pem_bwt_src/utils.cpp src/main.cpp $(AUX_DISK_FLAGS)

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f compute_bwt_sequential compute_bwt_parallel *.o
