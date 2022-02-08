#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/**
 * The structure used for carrying the program arguments back to the main
 * function. This avoids the use of global variables.
 */
typedef struct arguments {
    int size;
    double accuracy;
    int randSeed;
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
    ARGUMENTS *ar = processArgs(argc, argv);

    srand(ar->randSeed);

    int width = ar->size;
    int height = width;
    double accuracy = ar->accuracy;

    free(ar);

    printf("Size:%d,Accuracy:%f\n", width, accuracy);

    // For scalability testing, this was set as +-10, this reduces the number
    // of loops required, so we can see the effect of the parallelism.

    // These default values should result in the values settling roughly as the
    // centre of the array becomes non-zero.
    double randMax = 10.0;
    double randMin = -10.0;
    double randRange = (randMax - randMin);
    double randDiv = RAND_MAX / randRange;

    double **arr = createArray(height, width);
    double **arr2 = createArray(height, width);

    for (int i = 0; i < width; i++) {
        arr[0][i] = randMin + rand() / randDiv;
    }
    for (int i = 1; i < width - 1; i++) {
        arr[i][0] = randMin + rand() / randDiv;
        arr[i][width - 1] = randMin + rand() / randDiv;
        for (int j = 1; j < width - 1; j++) {
            arr[i][j] = 0.0;
        }
    }
    for (int i = 0; i < width; i++) {
        arr[width - 1][i] = randMin + rand() / randDiv;
    }

    free(arr2[0]);
    free(arr2[height - 1]);
    arr2[0] = arr[0];
    arr2[height - 1] = arr[height - 1];

    for (int i = 1; i < height - 1; i++) {
        arr2[i][0] = arr[i][0];
        arr2[i][width - 1] = arr[i][width - 1];
    }

    double *line = (double *) calloc(width,sizeof(double));
    int changed = 1;
    int loop = 0;

    struct timespec begin, end;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    while (changed) {
        changed = 0;
        loop++;

        line[0] = arr[1][0];
        line[width - 1] = arr[1][width - 1];
        for (int i = 1; i < width - 1; i++) {
            double val = (arr[0][i] + arr[2][i] + arr[1][i - 1] + arr[1][i + 1]) / 4;
            line[i] = val;

            if (!changed && (fabs((val - arr[1][i])) > accuracy)) changed = 1;
        }

        double *tmp = arr[1];
        arr[1] = line;
        line = tmp;

        for (int i = 2; i < height - 1; i++) {
            line[0] = arr[i][0];
            line[width - 1] = arr[i][width - 1];

            for (int j = 1; j < width - 1; j++) {
                double val = (line[j] + arr[i + 1][j] + arr[i][j - 1] + arr[i][j + 1]) / 4;
                line[j] = val;

                if (!changed && (fabs((val - arr[i][j])) > accuracy)) changed = 1;
            }

            tmp = arr[i];
            arr[i] = line;
            line = tmp;
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &end);

    free(line);

    freeInnerArrays(arr2 + 1, height - 2);
    free(arr2);

    printf("Loops: %d\n", loop);
    printf("Microseconds:%lu\n",
           (unsigned long) (((end.tv_sec - begin.tv_sec) * 1e6) +
                            ((end.tv_nsec - begin.tv_nsec) / 1e3)));
    printf("Sum: %1.5f\n", sumArray(arr, height, width));

    freeArray(arr, height);

    return 0;
}
