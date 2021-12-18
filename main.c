#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>

// The size of the array to compute (excluding edges)
unsigned long long int SIZE = 50;
// The accuracy to compute the array to
double ACCURACY = 0.01;

/**
 * Takes the input arguments and set the global variables appropriately.
 *
 * The available command line arguments are:
 *      -size (or -s)     : The size of the array to compute e.g. -s=1000 would
 *                          be a 1000 x 1000 array.
 *      -threads (or -t)  : The number of threads to use for the computation.
 *      -accuracy (or -a) : The accuracy the array must reach.
 *      -barriers (or -b) : Turns on the use of barriers.
 *      -start            : TThe number of threads to start with.
 *      -end              : TThe number of threads to end with.
 *
 * E.g. ./main -s=1000 -t=10 -a=0.001
 *
 * @param argc : The argc value for the program.
 * @param argv : The argv value for the program.
 */
void processArgs(int argc, char *argv[]) {
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
            SIZE = (unsigned long long int) strtoll(val, (char **)NULL, 10);
        }
        if (len > 6 && strncmp(argv[i], sizeStr, 6) == 0) {
            char val[len - 5];
            for (int j = 6; j < len; j++) {
                val[j - 6] = argv[i][j];
            }
            val[len - 6] = '\0';
            SIZE = (unsigned long long int) strtoll(val, (char **)NULL, 10);
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

    printf("Size:%d,Accuracy:%f\n", (int) (SIZE + 2), ACCURACY);
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

    // How to deal with border?
    unsigned long long startCell = (SIZE * SIZE / (unsigned long long) nproc) * myrank;
    if (myrank == (nproc - 1)) {
        unsigned long long endCell = SIZE + 1;
    } else {
        unsigned long long endCell = (SIZE * SIZE / (unsigned long long) nproc) * (myrank + 1);
    }

    // Instead of the following could look at using the old cell code with calc x and y
    // If y index == 0 then
    // take top value
    // else if x index % SIZE == 0 (this is same as == 0) then
    // take left value
    // else if x index % SIZE == (SIZE - 1) then
    // take right value
    // else if y index % SIZE == (SIZE - 1) then
    // take bottom value

    MPI_Finalize();
    return 0;
}
