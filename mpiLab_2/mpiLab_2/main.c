#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main(int argc, char* argv[]) {

    int rank, size;
    int mess = 0, order = 0;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    for (int i = 0; i < atoi(argv[1]) * size; i++) {
        MPI_Bcast(&mess, 1, MPI_INT, order, MPI_COMM_WORLD);
        if (order != rank)
            printf("%d -> %d: %d\n", order, rank, mess);
        mess += rank;
        MPI_Barrier(MPI_COMM_WORLD);
        order++;
        order %= 3;
    }

    MPI_Finalize();

    return 0;
}
