#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

/**
 * The structure used for carrying the program arguments back to the main
 * function. This avoids the use of global variables.
 */
typedef struct arguments {
    int size;
    double accuracy;
    int randSeed;
} ARGUMENTS;

typedef struct result {
    int loop;
    double **arr;
} RESULT;

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

double sumArray(double **arr, int height, int width)
{
    double sum = 0.0;

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            sum += arr[i][j];
        }
    }
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

void freeInnerArrays(double **arr, int height)
{
    for (int i = 0; i < height; i++) {
        free(arr[i]);
    }
}

void freeArray(double **arr, int height)
{
    freeInnerArrays(arr, height);
    free(arr);
}

int compareArrays(double **arr1, double **arr2, int height, int width)
{
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (arr1[i][j] != arr2[i][j]) {
                printf("Array check failed on: %d, %d with values %1.2f, %1.2f\n", i, j, arr1[i][j], arr2[i][j]);
                return 1;
            }
        }
    }
    return 0;
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
    ar->randSeed = time(NULL);

    char *sizeStr = "-size=";
    char *sStr = "-s=";
    char *accStr = "-accuracy=";
    char *aStr = "-a=";
    char *seedStr = "-seed=";
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
        if (len > 6 && strncmp(argv[i], seedStr, 6) == 0) {
            char val[len - 5];
            for (int j = 6; j < len; j++) {
                val[j - 6] = argv[i][j];
            }
            val[len - 6] = '\0';
            ar->randSeed = (int) strtoll(val, (char**)NULL, 10);
        }
    }

    return ar;
}

/**
 * The sequential version of the program.
 *
 * @param tmp1 : The input array.
 * @param height : The height of the input array.
 * @param width : The width of the input array.
 * @param accuracy : The accuracy to work to.
 * @return : The number of loops completed.
 */
RESULT *seq(double **in, int height, int width, double accuracy) {
    double **arr = in;
    double **arr2 = createArray(height, width);

    free(arr2[0]);
    free(arr2[height - 1]);
    arr2[0] = arr[0];
    arr2[height - 1] = arr[height - 1];

    for (int i = 1; i < height - 1; i++) {
        arr2[i][0] = arr[i][0];
        arr2[i][width - 1] = arr[i][width - 1];
    }

    int changed = 1;
    int loop = 0;

    while (changed) {
        changed = 0;
        loop++;

        for (int i = 1; i < height - 1; i++) {
            for (int j = 1; j < width - 1; j++) {
                double val = arr[i - 1][j];
                val += arr[i + 1][j];
                val += arr[i][j - 1];
                val += arr[i][j + 1];
                val /= 4;
                arr2[i][j] = val;

                if (!changed && (fabs((val - arr[i][j])) > accuracy)) changed = 1;
            }
        }

        double **tmp = arr;
        arr = arr2;
        arr2 = tmp;
    }

    freeInnerArrays(arr2 + 1, height - 2);
    free(arr2);

    RESULT *res = (RESULT *) malloc(sizeof(RESULT));
    res->loop = loop;
    res->arr = arr;
    return res;
}

/**
 * Calculates the values for the current iteration from the input array.
 *
 * @param in : The input array.
 * @param height : The height of the input array.
 * @param width : The width of the input array.
 * @param accuracy : The accuracy to work to.
 * @return : Whether any values changed by more than the accuracy.
 */
int calculate(double **in, int height, int width, double accuracy)
{
    double *line = (double *) calloc(width,sizeof(double));
    int changed = 0;

    line[0] = in[1][0];
    line[width - 1] = in[1][width - 1];
    for (int i = 1; i < width - 1; i++) {
        double val = (in[0][i] + in[2][i] + in[1][i - 1] + in[1][i + 1]) / 4;
        line[i] = val;

        if (!changed && (fabs((val - in[1][i])) > accuracy)) changed = 1;
    }

    double *tmp = in[1];
    in[1] = line;
    line = tmp;

    for (int i = 2; i < height - 1; i++) {
        line[0] = in[i][0];
        line[width - 1] = in[i][width - 1];

        for (int j = 1; j < width - 1; j++) {
            double val = (line[j] + in[i + 1][j] + in[i][j - 1] + in[i][j + 1]) / 4;
            line[j] = val;

            if (!changed && (fabs((val - in[i][j])) > accuracy)) changed = 1;
        }

        tmp = in[i];
        in[i] = line;
        line = tmp;
    }

    free(line);

    return changed;
}

int sendRow(double *row, int width, int dest, int tag)
{
    int w = width;
    int acc = 0;
    int rc;

    while (w > 0) {
        if (w > 30000) {
            rc = MPI_Send(row + acc, 30000, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
            if (rc != MPI_SUCCESS) return rc;
            w -= 30000;
            acc += 30000;
        } else {
            rc = MPI_Send(row + acc, w, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
            if (rc != MPI_SUCCESS) return rc;
            w -= w;
        }
    }

    return MPI_SUCCESS;
}

int receiveRow(double *row, int width, int dest, int tag)
{
    int w = width;
    int acc = 0;
    int rc;

    while (w > 0) {
        MPI_Status stat;
        if (w > 30000) {
            rc = MPI_Recv(row + acc, 30000, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &stat);
            if (rc != MPI_SUCCESS) return rc;
            w -= 30000;
            acc += 30000;
        } else {
            rc = MPI_Recv(row + acc, w, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &stat);
            if (rc != MPI_SUCCESS) return rc;
            w -= w;
        }
    }

    return MPI_SUCCESS;
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
    int w = width;
    int acc = 0;

    int reqCount = width / 30000;
    if (width % 30000) reqCount++;
    reqCount *= 2;
    if (myrank > 0 && myrank < nproc - 1) reqCount *= 2;

    MPI_Request reqs[reqCount];

    reqCount = 0;

    while (w > 0) {
        if (w > 30000) {
            if (myrank > 0) MPI_Isend(arr[1] + acc, 30000, MPI_DOUBLE, myrank - 1, loop, MPI_COMM_WORLD, &reqs[reqCount++]);
            if (myrank < nproc - 1) MPI_Isend(arr[height - 2] + acc, 30000, MPI_DOUBLE, myrank + 1, loop, MPI_COMM_WORLD, &reqs[reqCount++]);
            w -= 30000;
            acc += 30000;
        } else {
            if (myrank > 0) MPI_Isend(arr[1] + acc, w, MPI_DOUBLE, myrank - 1, loop, MPI_COMM_WORLD, &reqs[reqCount++]);
            if (myrank < nproc - 1) MPI_Isend(arr[height - 2] + acc, w, MPI_DOUBLE, myrank + 1, loop, MPI_COMM_WORLD, &reqs[reqCount++]);
            w -= w;
        }
    }

    w = width;
    acc = 0;

    while (w > 0) {
        if (w > 30000) {
            if (myrank > 0) MPI_Irecv(arr[0] + acc, 30000, MPI_DOUBLE, myrank - 1, loop, MPI_COMM_WORLD, &reqs[reqCount++]);
            if (myrank < nproc - 1) MPI_Irecv(arr[height - 1] + acc, 30000, MPI_DOUBLE, myrank + 1, loop, MPI_COMM_WORLD, &reqs[reqCount++]);
            w -= 30000;
            acc += 30000;
        } else {
            if (myrank > 0) MPI_Irecv(arr[0] + acc, w, MPI_DOUBLE, myrank - 1, loop, MPI_COMM_WORLD, &reqs[reqCount++]);
            if (myrank < nproc - 1) MPI_Irecv(arr[height - 1] + acc, w, MPI_DOUBLE, myrank + 1, loop, MPI_COMM_WORLD, &reqs[reqCount++]);
            w -= w;
        }
    }

    for (int i = 0; i < reqCount; i++) {
        MPI_Status stat;
        MPI_Wait(&reqs[i], &stat);
    }

    return 0;
}

double **distributeArray(double **arr, int width, int nproc)
{
    int lineCount = (width - 2) / nproc;
    int rowsCovered = lineCount * nproc;
    int rowsLeft = (width - 2) - rowsCovered;
    int height = lineCount + 2;
    int currHeight = height;
    if (0 < rowsLeft) currHeight++;

    double **retArr = createArray(currHeight, width);

    for (int i = 0; i < currHeight; i++) {
        for (int j = 0; j < width; j++) {
            retArr[i][j] = arr[i][j];
        }
    }

    int curr = currHeight - 1;
    for (int i = 1; i < nproc; i++) {
        currHeight = height;
        if (i < rowsLeft) currHeight++;

        for (int j = -1; j < currHeight - 1; j++) {
            int sendStat = sendRow(arr[curr + j], width, i, 0);
            if (sendStat != MPI_SUCCESS) {
                printf("Error receiving array.\n");
                MPI_Abort(MPI_COMM_WORLD, sendStat);
            }
        }

        curr += (currHeight - 2);
    }

    return retArr;
}

double **createAndDistributeArray(int width, int nproc, double accuracy)
{
    int lineCount = (width - 2) / nproc;
    int rowsCovered = lineCount * nproc;
    int rowsLeft = (width - 2) - rowsCovered;
    int currHeight = lineCount + 1;
    if (0 < rowsLeft) currHeight++;
    int nextHeight = currHeight;

    double **retArr = createArray(currHeight + 1, width);
    double *line = (double *) calloc(width, sizeof(double));
    for (int i = 1; i < width - 1; i++) {
        line[i] = 0.0;
    }

    // For scalability testing, this was set as +-10, this reduces the
    // number of loops required, so we can see the effect of the
    // parallelism.

    // These default values should result in the values settling
    // roughly as the centre of the array becomes non-zero.
    double max = (width / 2) * accuracy;
    double min = 0 - max;
    double range = (max - min);
    double div = RAND_MAX / range;

    for (int i = 0; i < width; i++) {
        retArr[0][i] = min + rand() / div;
    }

    for (int i = 1; i < currHeight; i++) {
        retArr[i][0] = min + rand() / div;
        retArr[i][width - 1] = min + rand() / div;
        for (int j = 1; j < width - 1; j++) {
            retArr[i][j] = 0.0;
        }
    }

    if (nproc > 1) {
        int sendStat = sendRow(retArr[currHeight - 1], width, 1, 0);
        if (sendStat != MPI_SUCCESS) {
            printf("Error sending array.\n");
            MPI_Abort(MPI_COMM_WORLD, sendStat);
        }
    }

    for (int i = 1; i < nproc; i++) {
        currHeight = lineCount;
        if (i < rowsLeft) currHeight++;

        for (int j = 0; j < currHeight; j++) {
            line[0] = min + rand() / div;
            line[width - 1] = min + rand() / div;

            if (j == 0) {
                if (i == 1) {
                    for (int k = 0; k < width; k++) {
                        retArr[nextHeight][k] = line[k];
                    }
                } else {
                    int sendStat = sendRow(line, width, i - 1, 0);
                    if (sendStat != MPI_SUCCESS) {
                        printf("Error sending array.\n");
                        MPI_Abort(MPI_COMM_WORLD, sendStat);
                    }
                }
            }

            if (j == currHeight - 1 && i != nproc - 1) {
                int sendStat = sendRow(line, width, i + 1, 0);
                if (sendStat != MPI_SUCCESS) {
                    printf("Error sending array.\n");
                    MPI_Abort(MPI_COMM_WORLD, sendStat);
                }
            }

            int sendStat = sendRow(line, width, i, 0);
            if (sendStat != MPI_SUCCESS) {
                printf("Error sending array.\n");
                MPI_Abort(MPI_COMM_WORLD, sendStat);
            }
        }
    }

    for (int i = 0; i < width; i++) {
        line[i] = min + rand() / div;
    }

    int sendStat = sendRow(line, width, nproc - 1, 0);
    if (sendStat != MPI_SUCCESS) {
        printf("Error sending array.\n");
        MPI_Abort(MPI_COMM_WORLD, sendStat);
    }

    return retArr;
}

double **gatherArray(double **arr, int height, int width, int myrank, int nproc)
{
    double **retArr;

    if (nproc == 1) {
        retArr = createArray(width, width);

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                retArr[i][j] = arr[i][j];
            }
        }

        return retArr;
    }

    if (myrank == 0) {
        retArr = createArray(width, width);

        for (int i = 0; i < height - 1; i++) {
            for (int j = 0; j < width; j++) {
                retArr[i][j] = arr[i][j];
            }
        }

        int curr = height - 1;
        for (int i = 1; i < nproc; i++) {
            int lineCount = (width - 2) / nproc;
            int rowsCovered = lineCount * nproc;
            int rowsLeft = (width - 2) - rowsCovered;
            if (i < rowsLeft) lineCount++;

            for (int j = 0; j < lineCount; j++) {
                int recvStat = receiveRow(retArr[curr], width, i, 0);
                if (recvStat != MPI_SUCCESS) {
                    printf("Error receiving array.\n");
                    MPI_Abort(MPI_COMM_WORLD, recvStat);
                }
                curr++;
            }
            if (i == nproc - 1) {
                int recvStat = receiveRow(retArr[curr], width, i, 0);
                if (recvStat != MPI_SUCCESS) {
                    printf("Error receiving array.\n");
                    MPI_Abort(MPI_COMM_WORLD, recvStat);
                }
                curr++;
            }
        }

        return retArr;
    }

    for (int i = 1; i < height - 1; i++) {
        int sendStat = sendRow(arr[i], width, 0, 0);
        if (sendStat != MPI_SUCCESS) {
            printf("Error receiving array.\n");
            MPI_Abort(MPI_COMM_WORLD, sendStat);
        }
    }
    if (myrank == nproc - 1) {
        int sendStat = sendRow(arr[height - 1], width, 0, 0);
        if (sendStat != MPI_SUCCESS) {
            printf("Error receiving array.\n");
            MPI_Abort(MPI_COMM_WORLD, sendStat);
        }
    }

    return NULL;
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

    srand(ar->randSeed);

    int width = ar->size;
    double accuracy = ar->accuracy;

    free(ar);

    if (myrank == 0) printf("Size:%d,Accuracy:%f\n", width, accuracy);

    int lineCount = (width - 2) / nproc;
    int rowsCovered = lineCount * nproc;
    int rowsLeft = (width - 2) - rowsCovered;

    int height = lineCount + 2;
    if (myrank < rowsLeft) height++;

    double **in;
    double **original;

    if (myrank == 0) {
        if (width <= 10000) {
            original = createArray(width, width);

            // For scalability testing, this was set as +-10, this reduces the
            // number of loops required, so we can see the effect of the
            // parallelism.

            // These default values should result in the values settling
            // roughly as the centre of the array becomes non-zero.
            double max = (width / 2) * accuracy;
            double min = 0 - max;
            double range = (max - min);
            double div = RAND_MAX / range;

            for (int i = 0; i < width; i++) {
                original[0][i] = min + rand() / div;
                original[width - 1][i] = min + rand() / div;
            }
            for (int i = 1; i < width - 1; i++) {
                original[i][0] = min + rand() / div;
                original[i][width - 1] = min + rand() / div;
                for (int j = 1; j < width - 1; j++) {
                    original[i][j] = 0.0;
                }
            }

            in = distributeArray(original, width, nproc);
        } else {
            in = createAndDistributeArray(width, nproc, accuracy);
        }
    } else {
        in = createArray(height, width);

        for (int i = 0; i < height; i++) {
            int recvStat = receiveRow(in[i], width, 0, 0);
            if (recvStat != MPI_SUCCESS) {
                printf("Error receiving array.\n");
                MPI_Abort(MPI_COMM_WORLD, recvStat);
            }
        }
    }

    int cont[1];
    cont[0] = 1;
    int loop = 0;

    struct timespec begin, end;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    while (cont[0]) {
        loop++;
        cont[0] = calculate(in, height, width, accuracy);
        MPI_Request statReq;
        if (myrank > 0) MPI_Isend(cont, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &statReq);
        int sendRes = sendAndReceiveResults(in, height, width, myrank, nproc, loop);
        if (sendRes != 0) {
            if (myrank == 0) printf("Error sending results.\n");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_UNKNOWN);
        }
        if (myrank == 0) {
            int tmpVal = cont[0];
            for (int i = 1; i < nproc; i++) {
                MPI_Status stat;
                MPI_Recv(cont, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &stat);
                if (!tmpVal && cont[0]) tmpVal = 1;
            }
            cont[0] = tmpVal;
        } else {
            MPI_Status stat;
            MPI_Wait(&statReq, &stat);
        }
        MPI_Bcast(cont, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    clock_gettime(CLOCK_MONOTONIC, &end);

    if (myrank == 0) {
        printf("Loops: %d\n", loop);
        printf("Microseconds:%lu\n",
               (unsigned long) (((end.tv_sec - begin.tv_sec) * 1e6) +
                                ((end.tv_nsec - begin.tv_nsec) / 1e3)));
    }

    // Used for correctness testing
    if (width <= 10000) {
        double **final = gatherArray(in, height, width, myrank, nproc);
        if (myrank == 0) {
            if (final == NULL) {
                printf("Error gathering final array.\n");
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_UNKNOWN);
            }
            RESULT *res = seq(original, width, width, accuracy);

            if (loop == res->loop) printf("Loop check succeeded.\n");
            else printf("Loop check failed.\n");

            int cmpRes = compareArrays(final, res->arr, width, width);

            if (cmpRes) printf("Array comparison failed.\n");
            else printf("Array comparison succeeded.\n");

            printf("Sum: %1.5f\n", sumArray(final, width, width));

            freeArray(final, width);
            freeArray(res->arr, width);
            free(res);
        }
    } else if (myrank == 0) {
        printf("Sum: %1.5f\n", sumArray(in, height, width));
    }

    freeArray(in, height);

    MPI_Finalize();
    return 0;
}
