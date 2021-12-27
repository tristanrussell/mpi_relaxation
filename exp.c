#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char **argv) {
    int rc, myrank, nproc;
    char name[MPI_MAX_PROCESSOR_NAME];
    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS) {
        printf ("Error starting MPI program\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    int size = 100000;
    double *arr = (double *) calloc(size, sizeof(double));

    if (myrank == 0) {
        for (int i = 0; i < size; i++) {
            arr[i] = 1.0;
        }
        MPI_Send(arr, size, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    } else {
        MPI_Status stat;
        MPI_Recv(arr, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &stat);
    }

    printf("Success: %d\n", myrank);

    MPI_Finalize();
    return 0;
}