#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

// The size of the array to compute (excluding edges)
int SIZE = 50;
// The accuracy to compute the array to
double ACCURACY = 0.01;

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
 * Free's the memory of the given 2 dimensional array. This free's both the
 * sub-arrays and the outer array.
 *
 * @param a : The 2 dimensional array to free.
 */
void freeArray(double **a, int height) {
    for (int i = 0; i < height; i++) {
        free(a[i]);
    }
    free(a);
}

/**
 * Builds a 2d array of zeros.
 *
 * @return : The 2d array of zeros.
 */
double **createArray(int height, int width) {
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
 * Takes the input arguments and set the global variables appropriately.
 *
 * The available command line arguments are:
 *      -size (or -s)     : The size of the array to compute e.g. -s=1000 would
 *                          be a 1000 x 1000 array.
 *      -accuracy (or -a) : The accuracy the array must reach.]
 *
 * E.g. ./main -s=1000 -t=10 -a=0.001
 *
 * @param argc : The argc value for the program.
 * @param argv : The argv value for the program.
 */
void processArgs(int argc, char *argv[], int myrank) {
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
            SIZE = (int) strtoll(val, (char **)NULL, 10);
        }
        if (len > 6 && strncmp(argv[i], sizeStr, 6) == 0) {
            char val[len - 5];
            for (int j = 6; j < len; j++) {
                val[j - 6] = argv[i][j];
            }
            val[len - 6] = '\0';
            SIZE = (int) strtoll(val, (char **)NULL, 10);
        }
        if (len > 3 && strncmp(argv[i], aStr, 3) == 0) {
            char val[len - 2];
            for (int j = 3; j < len; j++) {
                val[j - 3] = argv[i][j];
            }
            val[len - 3] = '\0';
            ACCURACY = strtod(val, (char**)NULL);
        }
        if (len > 10 && strncmp(argv[i], accStr, 10) == 0) {
            char val[len - 9];
            for (int j = 10; j < len; j++) {
                val[j - 10] = argv[i][j];
            }
            val[len - 10] = '\0';
            ACCURACY = strtod(val, (char**)NULL);
        }
    }

    if (myrank == 0) printf("Size:%d,Accuracy:%f\n", SIZE, ACCURACY);
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

// Sends a line of results to the next process, dir determines the direction,
// 0 is up, 1 is down
//void sendResults(double **arr, int height, int myRank, int nproc, int loop)
//{
//    MPI_Request topReq;
//    MPI_Request bottomReq;
//    MPI_Request notificationReq;
//
//    if (myRank > 0) MPI_Isend(arr[1], SIZE, MPI_DOUBLE, myRank - 1, loop, MPI_COMM_WORLD, &topReq);
//    if (myRank < nproc - 1) MPI_Isend(arr[height - 2], SIZE, MPI_DOUBLE, myRank + 1, loop, MPI_COMM_WORLD, &bottomReq);
//
//    MPI_Status stat;
//    if (myRank > 0) MPI_Recv(arr[0], SIZE, MPI_DOUBLE, myRank - 1, loop, MPI_COMM_WORLD, &stat);
//    if (myRank < nproc - 1) MPI_Recv(arr[height - 1], SIZE, MPI_DOUBLE, myRank + 1, loop, MPI_COMM_WORLD, &stat);
//}

//int sendFinal(double **arr, int height, int myRank, int nproc, int loop, int changed)
//{
//    if (myRank <= 0) return -1;
//
//    MPI_Request topReq;
//    MPI_Request bottomReq;
//    MPI_Request notificationReq;
//
//    MPI_Isend(arr[1], SIZE, MPI_DOUBLE, myRank - 1, loop, MPI_COMM_WORLD, &topReq);
//    if (myRank < nproc - 1) MPI_Isend(arr[height - 2], SIZE, MPI_DOUBLE, myRank + 1, loop, MPI_COMM_WORLD, &bottomReq);
//    MPI_Isend(&changed, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &notificationReq);
//
//    MPI_Request reqs[height - 2];
//
//    for (int i = 1; i < height - 1; i++) {
//        MPI_Isend(arr[i], SIZE, MPI_DOUBLE, 0, loop, MPI_COMM_WORLD, &reqs[i - 1]);
//    }
//
//    int *cont = (int*) malloc(sizeof(int));
//
//    MPI_Bcast(cont, 1, MPI_INT, 0, MPI_COMM_WORLD);
//
//    int ret = cont[0];
//    free(cont);
//
//    if (ret) {
//        for (int i = 1; i < height - 1; i++) {
//            MPI_Cancel(&reqs[i - 1]);
//            MPI_Request_free(&reqs[i - 1]);
//        }
//    }
//
//    MPI_Status stat;
//    MPI_Wait(&topReq, &stat);
//
//    if (myRank < SIZE - 1) {
//        MPI_Status stat;
//        MPI_Wait(&bottomReq, &stat);
//    }
//
//    return ret;
//}

int sendAndReceive(double **arr, int height, int myRank, int nproc, int loop, int changed)
{
    MPI_Request topReq;
    MPI_Request bottomReq;

    int cont[1];
    cont[0] = changed;

    if (myRank > 0) {
        MPI_Send(cont, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Isend(arr[1], SIZE, MPI_DOUBLE, myRank - 1, loop, MPI_COMM_WORLD, &topReq);
    }
    if (myRank < nproc - 1) MPI_Isend(arr[height - 2], SIZE, MPI_DOUBLE, myRank + 1, loop, MPI_COMM_WORLD, &bottomReq);

    if (myRank == 0) {
        for (int i = 1; i < nproc; i++) {
            int j[1];
            MPI_Status stat;
            MPI_Recv(j, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &stat);
            if (j[0] == 1) {
                cont[0] = 1;
            }
        }
    }

    MPI_Bcast(cont, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (cont[0]) {
        MPI_Status stat;
        if (myRank > 0) MPI_Recv(arr[0], SIZE, MPI_DOUBLE, myRank - 1, loop, MPI_COMM_WORLD, &stat);
        if (myRank < nproc - 1) MPI_Recv(arr[height - 1], SIZE, MPI_DOUBLE, myRank + 1, loop, MPI_COMM_WORLD, &stat);

        if (myRank > 0) MPI_Wait(&topReq, &stat);
        if (myRank < nproc - 1) MPI_Wait(&bottomReq, &stat);
        return 1;
    } else {
        for (int i = 1; i < height - 1; i++) {
            MPI_Send(arr[i], SIZE, MPI_DOUBLE, 0, loop, MPI_COMM_WORLD);
        }
        if (myRank == nproc - 1) {
            MPI_Send(arr[height - 1], SIZE, MPI_DOUBLE, 0, loop, MPI_COMM_WORLD);
        }
        return 0;
    }
}

//int checkFinal(double *line, int changed, int myRank, int nproc)
//{
//    int *cont = (int*) malloc(sizeof(int));
//
//    cont[0] = changed;
//    for (int i = 1; i < nproc; i++) {
//        int *j = (int*) malloc(sizeof(int));
//        MPI_Status stat;
//        MPI_Recv(j, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &stat);
//        if (j[0] == 1) {
//            cont[0] = 1;
//        }
//        free(j);
//    }
//
//    MPI_Bcast(cont, 1, MPI_INT, 0, MPI_COMM_WORLD);
//
//    int ret = cont[0];
//    free(cont);
//
//    return ret;
//}

double **getFinalArray(double **arr, int height, int width, int nproc, int myRank, int loop)
{
    if (myRank != 0) return NULL;

    double **j = createArray(SIZE, SIZE);

    for (int i = 0; i < height - 1; i++) {
        free(j[i]);
        j[i] = arr[i];
    }

    for (int i = 1; i < nproc; i++) {
        int startRow = (width / nproc);
        if (width % nproc) startRow++;
        startRow = startRow * i;

        int endRow;
        if (i == (nproc - 1)) {
            endRow = width + 1;
        } else {
            endRow = (width / nproc);
            if (width % nproc) endRow++;
            endRow = endRow  * (i + 1);
        }

        int iHeight = endRow - startRow;

        for (int k = startRow; k < endRow; k++) {
            MPI_Status stat;
            MPI_Recv(j[startRow], width, MPI_DOUBLE, i, loop, MPI_COMM_WORLD, &stat);
        }
    }

    return j;
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

    processArgs(argc, argv, myrank);

    int width = SIZE;
    double accuracy = ACCURACY;

    int lineCount = (width - 2) / nproc;
    int rowsCovered = lineCount * nproc;
    int rowsLeft = (width - 2) - rowsCovered;

    int height = lineCount + 2;
    if (myrank < rowsLeft) height++;

    // Can try to improve later to add per cell rather than per line
    // How to deal with border?
//    int round = SIZE % nproc;
//
//    int startRow;
//    if (myrank == 0) {
//        startRow = 1;
//    } else {
//        startRow = (SIZE / nproc);
//        if (round) startRow++;
//        startRow = startRow * myrank;
//    }
//
//    int endRow;
//    if (myrank == (nproc - 1)) {
//        endRow = SIZE;
//    } else {
//        endRow = (SIZE / nproc);
//        if (round) endRow++;
//        endRow = endRow  * (myrank + 1);
//    }
//
//    int height = endRow - startRow + 2;

    // Arrays for each process will store top, left, right, bottom and the internal values

    // Instead of the following could look at using the old cell code with calc x and y
    // If y index == 0 then
    // take top value
    // else if x index % SIZE == 0 (this is same as == 0) then
    // take left value
    // else if x index % SIZE == (SIZE - 1) then
    // take right value
    // else if y index % SIZE == (SIZE - 1) then
    // take bottom value

    double **in = createArray(height, SIZE);
    double **out = createArray(height, SIZE);

    if (myrank == 0) {
        for (int i = 0; i < SIZE; i++) {
            in[0][i] = 1.0;
        }
        for (int i = 1; i < height; i++) {
            in[i][0] = 1.0;
            in[i][SIZE - 1] = 1.0;
            for (int j = 1; j < SIZE - 1; j++) {
                in[i][j] = 0.0;
            }
        }
    } else if (myrank == nproc - 1) {
        for (int i = 0; i < SIZE; i++) {
            in[height - 1][i] = 1.0;
        }
        for (int i = 0; i < height - 1; i++) {
            in[i][0] = 1.0;
            in[i][SIZE - 1] = 1.0;
            for (int j = 1; j < SIZE - 1; j++) {
                in[i][j] = 0.0;
            }
        }
    } else {
        for (int i = 0; i < height - 1; i++) {
            in[i][0] = 1.0;
            in[i][SIZE - 1] = 1.0;
            for (int j = 1; j < SIZE - 1; j++) {
                in[i][j] = 0.0;
            }
        }
    }

    int changed = 1;
    int loop = 0;
    while (changed) {
        loop++;
        changed = calculate(in, out, height, SIZE, accuracy);
        changed = sendAndReceive(out, height, myrank, nproc, loop, changed);
        double **tmp = out;
        out = in;
        in = tmp;
    }

    if (myrank == 0) {
        double **final;
        final = getFinalArray(in, height, width, nproc, myrank, loop);
        printArray(final, width, width);
    }

//    double **final;
//    if (myrank == 0) {
//        final = getFinalArray(out, height, nproc, myrank, loop);
//        // Do stuff with final array
//        freeArray(final, height);
//        freeArray(in + 1, height - 2);
//    } else {
//        freeArray(in, height);
//        freeArray(out + 1, height - 2);
//    }
//    freeArray(in, height);
//    freeArray(out + 1, height - 2);

    // if significant change then
    //   sendResults(); - Should also signal to continue
    //   receiveResults();
    //   ensure sendResults() has finished then next loop
    // else
    //   sendFinal(); - This sends to 0 and its partner like sendResults does
    //   receiveFinal();
    //   if receiveFinal == 1 then
    //     wait for sendFinal(); to finish then next loop

    // Might not need to actually get the array at end

    MPI_Finalize();
    return 0;
}
