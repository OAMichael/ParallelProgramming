#include <chrono>

#include "3.1_MatrixMul.h"

#define MATRIX_SIZE 1000
//#define OUTPUT_MATRIX
//#define MAKE_TEST

#if defined(MATRIX_SIZE) && (MATRIX_SIZE < 10) && !defined(MAKE_TEST)
#define PRINT_MATRIX
#endif

using MatrixType = double;

int main() {

#ifndef MAKE_TEST
    Matrix<MatrixType, MATRIX_SIZE> A;
    Matrix<MatrixType, MATRIX_SIZE> B;

    for (int i = 0; i < A.getSize(); ++i) {
        for (int j = 0; j < A.getSize(); ++j) {
            A[i][j] = i + j;
        }
    }

    for (int i = 0; i < B.getSize(); ++i) {
        for (int j = 0; j < B.getSize(); ++j) {
            B[i][j] = i + j;
        }
    }

    Matrix<MatrixType, MATRIX_SIZE> C;

#ifdef PRINT_MATRIX
    printf("A:\n");
    printMatrix(A);
    printf("\nB:\n");
    printMatrix(B);
    printf("\n");
#endif  // PRINT_MATRIX

#ifdef OUTPUT_MATRIX
    std::ofstream outFileA;
    outFileA.open("./3.1_data/A.dat", std::ios::out | std::ios::trunc);
    printMatrix(A, outFileA);    
    outFileA.close();

    std::ofstream outFileB;
    outFileB.open("./3.1_data/B.dat", std::ios::out | std::ios::trunc);
    printMatrix(B, outFileB);    
    outFileB.close();
#endif  // OUTPUT_MATRIX


    {
        auto start = std::chrono::high_resolution_clock::now();

        naiveMatrixMultiplication(A, B, C);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();
        
        printf("Elapsed naive:                      %lf seconds\n", microseconds / 1000000.0);

#ifdef PRINT_MATRIX
        printf("C:\n");
        printMatrix(C);
        printf("\n");
#endif  // PRINT_MATRIX
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        transposeMatrixMultiplication(A, B, C);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();

        printf("Elapsed transpose:                  %lf seconds\n", microseconds / 1000000.0);

#ifdef PRINT_MATRIX
        printf("C:\n");
        printMatrix(C);
        printf("\n");
#endif  // PRINT_MATRIX
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        sumAndPrecopyMatrixMultiplication(A, B, C);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();
        
        printf("Elapsed sum and precopy:            %lf seconds\n", microseconds / 1000000.0);

#ifdef PRINT_MATRIX
        printf("C:\n");
        printMatrix(C);
        printf("\n");
#endif  // PRINT_MATRIX
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        naiveMatrixMultiplication_Parallel(A, B, C);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();
        
        printf("Elapsed naive parallel:             %lf seconds\n", microseconds / 1000000.0);

#ifdef PRINT_MATRIX
        printf("C:\n");
        printMatrix(C);
        printf("\n");
#endif  // PRINT_MATRIX
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        transposeMatrixMultiplication_Parallel(A, B, C);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();

        printf("Elapsed transpose parallel:         %lf seconds\n", microseconds / 1000000.0);

#ifdef PRINT_MATRIX
        printf("C:\n");
        printMatrix(C);
        printf("\n");
#endif  // PRINT_MATRIX
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        sumAndPrecopyMatrixMultiplication_Parallel(A, B, C);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();
        
        printf("Elapsed sum and precopy parallel:   %lf seconds\n", microseconds / 1000000.0);

#ifdef PRINT_MATRIX
        printf("C:\n");
        printMatrix(C);
        printf("\n");
#endif  // PRINT_MATRIX
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        naiveMatrixMultiplication_ParallelCollapse(A, B, C);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();
        
        printf("Elapsed naive parallel collapse:    %lf seconds\n", microseconds / 1000000.0);

#ifdef PRINT_MATRIX
        printf("C:\n");
        printMatrix(C);
        printf("\n");
#endif  // PRINT_MATRIX

#ifdef OUTPUT_MATRIX
        std::ofstream outFileC;
        outFileC.open("./3.1_data/C.dat", std::ios::out | std::ios::trunc);
        printMatrix(C, outFileC);    
        outFileC.close();
#endif  // OUTPUT_MATRIX
    }


#else   // ^^^ !MAKE_TEST ^^^  vvv MAKE_TEST vvv

    printf("n m1 m2 m3 m4 m5 m6 m7\n");

    testFunction<2,    double>();
    testFunction<3,    double>();
    testFunction<4,    double>();
    testFunction<5,    double>();
    testFunction<10,   double>();
    testFunction<15,   double>();
    testFunction<20,   double>();
    testFunction<30,   double>();
    testFunction<40,   double>();
    testFunction<50,   double>();
    testFunction<75,   double>();
    testFunction<100,  double>();
    testFunction<150,  double>();
    testFunction<200,  double>();
    testFunction<300,  double>();
    testFunction<400,  double>();
    testFunction<500,  double>();
    testFunction<750,  double>();
    testFunction<1000, double>();
    testFunction<1250, double>();
    testFunction<1600, double>();

#endif  

    return 0;
}
