#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/**
 * Builds a 2d array of zeros.
 *
 * @return : The 2d array of zeros.
 */
double **createArray(int height, int width)
{
    double **arr = (double**)malloc((unsigned long)height * sizeof(double*));
    for (int i = 0; i < height; i++) {
        arr[i] = (double*)calloc(width, sizeof(double));

        for (int j = 0; j < width; j++) {
            arr[i][j] = 0.00;
        }
    }

    return arr;
}

void printArray(double **arr, int height, int width)
{
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            printf("%1.2f,", arr[i][j]);
        }
        printf("\n");
    }
}

int calculate(double **in, double **out, int height, int width, double accuracy)
{
    int changed = 0;
    out[0] = in[0];
    out[height - 1] = in[height - 1];
    for (int i = 1; i < height - 1; i++) {
        out[i][0] = in[i][0];
        out[i][width - 1] = in[i][width - 1];
        for (int j = 1; j < width - 1; j++) {
            double val = in[i - 1][j];
            val += in[i + 1][j];
            val += in[i][j - 1];
            val += in[i][j + 1];
            val /= 4;
            out[i][j] = val;

            if (!changed && (fabs((val - in[i][j])) > accuracy)) changed = 1;
        }
    }
    return changed;
}

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

    int width = 50;
    double accuracy = 0.01;

    int lineCount = (width - 2) / nproc;
    int rowsCovered = lineCount * nproc;
    int rowsLeft = (width - 2) - rowsCovered;

    int height = lineCount + 2;
    if (myrank < rowsLeft) height++;

    double **in = createArray(height, width);
    double **out = createArray(height, width);

    if (myrank == 0) {
        for (int i = 0; i < width; i++) {
            in[0][i] = 1.0;
        }
        for (int i = 1; i < height; i++) {
            in[i][0] = 1.0;
            in[i][width - 1] = 1.0;
            for (int j = 1; j < width - 1; j++) {
                in[i][j] = 0.0;
            }
        }
    } else if (myrank == nproc - 1) {
        for (int i = 0; i < width; i++) {
            in[height - 1][i] = 1.0;
        }
        for (int i = 0; i < height - 1; i++) {
            in[i][0] = 1.0;
            in[i][width - 1] = 1.0;
            for (int j = 1; j < width - 1; j++) {
                in[i][j] = 0.0;
            }
        }
    } else {
        for (int i = 0; i < height; i++) {
            in[i][0] = 1.0;
            in[i][width - 1] = 1.0;
            for (int j = 1; j < width - 1; j++) {
                in[i][j] = 0.0;
            }
        }
    }

    calculate(in, out, height, width, accuracy);

    int cont[1];
    cont[0] = 1;

    int loop = 0;
    while (cont[0]) {
        loop++;
        cont[0] = calculate(in, out, height, width, accuracy);
        if (myrank > 0) MPI_Send(out[1], width, MPI_DOUBLE, myrank - 1, loop, MPI_COMM_WORLD);
        if (myrank < nproc - 1) MPI_Send(out[height - 2], width, MPI_DOUBLE, myrank + 1, loop, MPI_COMM_WORLD);
        MPI_Status upStat;
        if (myrank > 0) MPI_Recv(out[0], width, MPI_DOUBLE, myrank - 1, loop, MPI_COMM_WORLD, &upStat);
        MPI_Status downStat;
        if (myrank < nproc - 1) MPI_Recv(out[height - 1], width, MPI_DOUBLE, myrank + 1, loop, MPI_COMM_WORLD, &downStat);
        if (myrank > 0) MPI_Send(cont, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        double **tmp = out;
        out = in;
        in = tmp;
        if (myrank == 0) {
            int tmpVal = cont[0];
            for (int i = 1; i < nproc; i++) {
                MPI_Status stat;
                MPI_Recv(cont, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &stat);
                if (!tmpVal && cont[0]) tmpVal = 1;
            }
            cont[0] = tmpVal;
        }
        MPI_Bcast(cont, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    if (myrank == 0) printf("Loops: %d\n", loop);

    if (myrank < 6) {
        printf("Process: %d\n", myrank);
        printArray(in, height, width);
    }

    MPI_Finalize();
    return 0;
}
