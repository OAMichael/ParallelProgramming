#include <chrono>
#include <random>


#include "5.1_StrassenImpl.h"

#define MATRIX_SIZE 1024

//#define OUTPUT_MATRIX

using MatrixType = float;

int main() {

    Matrix<MatrixType> A{MATRIX_SIZE};
    Matrix<MatrixType> B{MATRIX_SIZE};

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<> random_distribution(0, 1);


    for (int i = 0; i < A.getSize(); ++i) {
        for (int j = 0; j < A.getSize(); ++j) {
            A[i][j] = random_distribution(rng);
        }
    }

    for (int i = 0; i < B.getSize(); ++i) {
        for (int j = 0; j < B.getSize(); ++j) {
            B[i][j] = random_distribution(rng);
        }
    }

#ifdef OUTPUT_MATRIX
    std::ofstream outFileA;
    outFileA.open("./5.1_data/A.dat", std::ios::out | std::ios::trunc);
    printMatrix(A, outFileA);    
    outFileA.close();

    std::ofstream outFileB;
    outFileB.open("./5.1_data/B.dat", std::ios::out | std::ios::trunc);
    printMatrix(B, outFileB);    
    outFileB.close();
#endif

    omp_set_nested(1);
    omp_set_num_threads(NUM_THREADS);
    
    {
        Matrix<MatrixType> C{MATRIX_SIZE};

        auto start = std::chrono::high_resolution_clock::now();

        MatrixMul_Strassen(A, B, C);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();

        printf("Elapsed Strassen (MatSize = %d, NumThreads = %d):  %lf seconds\n", MATRIX_SIZE, NUM_THREADS, microseconds / 1000000.0);

#ifdef OUTPUT_MATRIX
        std::ofstream outFileC;
        outFileC.open("./5.1_data/C.dat", std::ios::out | std::ios::trunc);
        printMatrix(C, outFileC);    
        outFileC.close();
#endif
    }

    return 0;
}
