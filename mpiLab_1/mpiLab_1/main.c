#include <stdio.h>
#include "mpi.h"

int main(int argc, char* argv[]) {

    int rank, size;
    int mess, recv = 0;
    MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        mess = recv;
        MPI_Send(&mess, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    }
    
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < size; j++) {
            if (rank == j && ((rank + 1) < size)) {
                if (rank == 0) MPI_Recv(&recv, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, &status);
                else MPI_Recv(&recv, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);
                
                printf("%d -> %d: %d\n", status.MPI_SOURCE, rank, recv);

                mess = recv + rank;
                MPI_Send(&mess, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
            } else if (rank == j && ((rank + 1) == size)) {
                MPI_Recv(&recv, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);
                printf("%d -> %d: %d\n", status.MPI_SOURCE, rank, recv);

                mess = recv + rank;
                MPI_Send(&mess, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
        }

    MPI_Finalize();

    return 0;
}
