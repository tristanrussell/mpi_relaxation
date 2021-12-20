#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

/**
 * The structure used for carrying the program arguments back to the main
 * function. This avoids the use of global variables.
 */
typedef struct arguments {
    int size;
    double accuracy;
} ARGUMENTS;

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

/**
 * Prints out an array of values to the console.
 *
 * @param arr : The array to print.
 * @param height : The height of the array.
 * @param width : The width of the array.
 */
void printArray(double **arr, int height, int width)
{
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            printf("%1.2f,", arr[i][j]);
        }
        printf("\n");
    }
}

/**
 * Takes the input arguments and set the global variables appropriately.
 *
 * The available command line arguments are:
 *      -size (or -s)     : The size of the array to compute e.g. -s=1000 would
 *                          be a 1000 x 1000 array.
 *      -accuracy (or -a) : The accuracy the array must reach.
 *
 * E.g. ./main -s=1000 -t=10 -a=0.001
 *
 * @param argc : The argc value for the program.
 * @param argv : The argv value for the program.
 * @return : The arguments in an ARGUMENT structure.
 */
ARGUMENTS *processArgs(int argc, char *argv[]) {
    ARGUMENTS *ar = (ARGUMENTS*)malloc(sizeof(ARGUMENTS));

    // Set defaults
    ar->size = 500;
    ar->accuracy = 0.01;

    char *sizeStr = "-size=";
    char *sStr = "-s=";
    char *accStr = "-accuracy=";
    char *aStr = "-a=";
    for (int i = 1; i < argc; i++) {
        int len = (int)strlen(argv[i]);
        if (len > 3 && strncmp(argv[i], sStr, 3) == 0) {
            char val[len - 2];
            for (int j = 3; j < len; j++) {
                val[j - 3] = argv[i][j];
            }
            val[len - 3] = '\0';
            ar->size = (int) strtoll(val, (char **)NULL, 10);
        }
        if (len > 6 && strncmp(argv[i], sizeStr, 6) == 0) {
            char val[len - 5];
            for (int j = 6; j < len; j++) {
                val[j - 6] = argv[i][j];
            }
            val[len - 6] = '\0';
            ar->size = (int) strtoll(val, (char **)NULL, 10);
        }
        if (len > 3 && strncmp(argv[i], aStr, 3) == 0) {
            char val[len - 2];
            for (int j = 3; j < len; j++) {
                val[j - 3] = argv[i][j];
            }
            val[len - 3] = '\0';
            ar->accuracy = strtod(val, (char**)NULL);
        }
        if (len > 10 && strncmp(argv[i], accStr, 10) == 0) {
            char val[len - 9];
            for (int j = 10; j < len; j++) {
                val[j - 10] = argv[i][j];
            }
            val[len - 10] = '\0';
            ar->accuracy = strtod(val, (char**)NULL);
        }
    }

    return ar;
}

/**
 * Calculates the values for the current iteration from the input array.
 *
 * @param in : The input array.
 * @param out : The output array.
 * @param height : The height of the input and output arrays.
 * @param width : The width of the input and output arrays.
 * @param accuracy : The accuracy to work to.
 * @return : Whether any values changed by more than the accuracy.
 */
int calculate(double **in, double **out, int height, int width, double accuracy)
{
    int changed = 0;
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

/**
 * Sends the results to the processes that needs the results and retrieves the
 * results from other processes and replaces the current values the output
 * array. The messages are divided into 30000 chunks (if greater) because
 * messages larger than around 31000 would cause the send to fail, so 30000
 * was chosen to leave a little headroom and as the closest round number.
 *
 * @param arr : The output array from the current iteration.
 * @param height : The height of the input array.
 * @param width : The width of the input array.
 * @param myrank : The rank of the current process.
 * @param nproc : The number of processes.
 * @param loop : The current iteration loop count.
 * @return : 0 if successful and 1 otherwise.
 */
int sendAndReceiveResults(double **arr, int height, int width, int myrank, int nproc, int loop)
{
//    Original:
//    if (myrank > 0) MPI_Send(arr[1], width, MPI_DOUBLE, myrank - 1, loop, MPI_COMM_WORLD);
//    if (myrank < nproc - 1) MPI_Send(arr[height - 2], width, MPI_DOUBLE, myrank + 1, loop, MPI_COMM_WORLD);
//    MPI_Status upStat;
//    if (myrank > 0) MPI_Recv(arr[0], width, MPI_DOUBLE, myrank - 1, loop, MPI_COMM_WORLD, &upStat);
//    MPI_Status downStat;
//    if (myrank < nproc - 1) MPI_Recv(arr[height - 1], width, MPI_DOUBLE, myrank + 1, loop, MPI_COMM_WORLD, &downStat);
//    return 0;

    int w = width;
    int acc = 0;

    while (w > 0) {
        if (w > 30000) {
            if (myrank > 0) MPI_Send(arr[1] + acc, 30000, MPI_DOUBLE, myrank - 1, loop, MPI_COMM_WORLD);
            if (myrank < nproc - 1) MPI_Send(arr[height - 2] + acc, 30000, MPI_DOUBLE, myrank + 1, loop, MPI_COMM_WORLD);
            w -= 30000;
            acc += 30000;
        } else {
            if (myrank > 0) MPI_Send(arr[1] + acc, w, MPI_DOUBLE, myrank - 1, loop, MPI_COMM_WORLD);
            if (myrank < nproc - 1) MPI_Send(arr[height - 2] + acc, w, MPI_DOUBLE, myrank + 1, loop, MPI_COMM_WORLD);
            w -= w;
        }
    }

    w = width;
    acc = 0;

    while (w > 0) {
        if (w > 30000) {
            MPI_Status upStat;
            if (myrank > 0) MPI_Recv(arr[0] + acc, 30000, MPI_DOUBLE, myrank - 1, loop, MPI_COMM_WORLD, &upStat);
            MPI_Status downStat;
            if (myrank < nproc - 1) MPI_Recv(arr[height - 1] + acc, 30000, MPI_DOUBLE, myrank + 1, loop, MPI_COMM_WORLD, &downStat);
            w -= 30000;
            acc += 30000;
        } else {
            MPI_Status upStat;
            if (myrank > 0) MPI_Recv(arr[0] + acc, w, MPI_DOUBLE, myrank - 1, loop, MPI_COMM_WORLD, &upStat);
            MPI_Status downStat;
            if (myrank < nproc - 1) MPI_Recv(arr[height - 1] + acc, w, MPI_DOUBLE, myrank + 1, loop, MPI_COMM_WORLD, &downStat);
            w -= w;
        }
    }

    return 0;
}

/**
 * This program takes a square array of values and iterates across the array
 * replacing the cells with the sum of its 4 neighbours. The edge values are
 * fixed.
 *
 * @param argc
 * @param argv
 * @return
 */
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

    ARGUMENTS *ar = processArgs(argc, argv);

    int width = ar->size;
    double accuracy = ar->accuracy;

    free(ar);

    if (myrank == 0) printf("Size:%d,Accuracy:%f\n", width, accuracy);

    int lineCount = (width - 2) / nproc;
    int rowsCovered = lineCount * nproc;
    int rowsLeft = (width - 2) - rowsCovered;

    int height = lineCount + 2;
    if (myrank < rowsLeft) height++;

    double **in = createArray(height, width);
    double **out = createArray(height, width);
    free(out[0]);
    free(out[height - 1]);
    out[0] = in[0];
    out[height - 1] = in[height - 1];

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

    int cont[1];
    cont[0] = 1;

    int loop = 0;
    while (cont[0]) {
        loop++;
        cont[0] = calculate(in, out, height, width, accuracy);
        int sendRes = sendAndReceiveResults(out, height, width, myrank, nproc, loop);
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
    printf("Rank: %d, height: %d\n", myrank, height);

//    if (myrank < 6) {
//        printf("Process: %d\n", myrank);
//        printArray(in, height, width);
//    }

    MPI_Finalize();
    return 0;
}
