#include <iostream>
#include <cmath>

#include <omp.h>

#include "Common.h"


int main(int argc, char **argv) {

    My2DArray<double> a(ISIZE, JSIZE);

    for (int i = 0; i < ISIZE; ++i) {
        for (int j = 0; j < JSIZE; ++j) {
            a[i][j] = 10 * i + j;
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

    std::string outFilename = "baseline_OMP_" + std::to_string(num_threads) + ".dat";

    // Start time
    double start = omp_get_wtime();

    #pragma omp parallel for
    for (int i = 0; i < ISIZE; ++i) {
        for (int j = 0; j < JSIZE; ++j) {
            a[i][j] = std::sin(2 * a[i][j]);
        }
    }

    // End time
    double end = omp_get_wtime();

    std::cout << "Number of executors: " << num_threads << std::endl;
    std::cout << "Elapsed time: " << end - start << " seconds" << std::endl;
    
    serializeToFile(outFilename, a);
    return 0;
}
