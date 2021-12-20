#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv)
{
    int rc, myrank, nproc;
    char name[MPI_MAX_PROCESSOR_NAME];
    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS) {
        printf ("Error starting MPI program\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    double arr[1];

    if (myrank == 1) {
        arr[0] = 1.0;
        MPI_Send(arr, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        arr[0] = 2.0;
        MPI_Send(arr, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    if (myrank == 3) {
        arr[0] = 1.0;
        MPI_Send(arr, 1, MPI_DOUBLE, 4, 0, MPI_COMM_WORLD);
    }

    int cont[1];
    if (myrank == 2) cont[0] = 1;
    MPI_Bcast(cont, 1, MPI_INT, 2, MPI_COMM_WORLD);

    if (myrank == 0) {
        MPI_Status stat;
        MPI_Recv(arr, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &stat);
        printf("Process 0, Received 1: %1.2f\n", arr[0]);
        arr[0] = 0.0;
        MPI_Status stat2;
        MPI_Recv(arr, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &stat2);
        printf("Process 0, Received 2: %1.2f\n", arr[0]);
    }

    if (myrank == 4) {
        MPI_Status stat;
        MPI_Recv(arr, 1, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD, &stat);
        printf("Process: 4, Received 1: %1.2f\n", arr[0]);
        arr[0] = 0.0;
        MPI_Status stat2;
        MPI_Recv(arr, 1, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD, &stat2);
        printf("Process 4, Received 2: %1.2f\n", arr[0]);
    }

    MPI_Finalize();
    return 0;
}
