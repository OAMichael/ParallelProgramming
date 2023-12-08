#include <iostream>
#include <cmath>
#include <cstring>

#include <omp.h>

#include "Common.h"


int main(int argc, char **argv) {

    My2DArray<double> a(ISIZE, JSIZE);
    My2DArray<double> b(ISIZE, JSIZE);

    for (int i = 0; i < ISIZE; ++i) {
        for (int j = 0; j < JSIZE; ++j) {
            a[i][j] = 10 * i + j;
            b[i][j] = 0;
        }
    }

    int num_threads = 1;
    if (argc >= 2) {
        num_threads = std::atoi(argv[1]);
        if (num_threads <= 0) {
            std::cerr << "Number of threads must be > 0!" << std::endl;
            return -1;
        }
    }

    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);

    std::string outFilename = "task3_OMP_" + std::to_string(num_threads) + ".dat";

    double tmpArr[JSIZE];

    // Start time
    double start = omp_get_wtime();

    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < ISIZE; ++i) {
            for (int j = 0; j < JSIZE; ++j) {
                a[i][j] = std::sin(0.005 * a[i][j]);
            }
        }
        #pragma omp barrier

        #pragma omp for
        for (int i = 5; i < ISIZE; ++i) {
            for (int j = 0; j < JSIZE - 2; ++j) {
                b[i][j] = a[i - 5][j + 2] * 1.5;
            }
        }
    }

    // End time
    double end = omp_get_wtime();

    std::cout << "Number of executors: " << num_threads << std::endl;
    std::cout << "Elapsed time: " << end - start << " seconds" << std::endl;

    serializeToFile(outFilename, b);
    return 0;
}
