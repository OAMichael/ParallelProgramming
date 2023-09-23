#ifndef MATRIX_MUL_H
#define MATRIX_MUL_H

#include <fstream>
#include <iostream>
#include <cstring>

#include <omp.h>


template<typename T, int n = 1>
struct Matrix {

    T* data;            // data stored in row-major order

    inline int getSize()  const { return n; }
    inline void clear() { std::memset(data, 0, n * n * sizeof(T)); }

    inline const T* operator[](const int i) const { return &data[n * i]; }
    inline T* operator[](const int i) { return data + n * i; }

    Matrix() {
        data = new T[n * n];
    }

    ~Matrix() {
        delete[] data;
    }
};



template<typename T, int n>
std::ostream& printMatrix(const Matrix<T, n>& mat, std::ostream& outFile = std::cout) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            outFile << mat[i][j] << " ";
        }
        outFile << std::endl;
    }

    return outFile;
}


template<typename T, int n>
void naiveMatrixMultiplication(const Matrix<T, n>& A, const Matrix<T, n>& B, Matrix<T, n>& C) {

    C.clear();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}


template<typename T, int n>
void transposeMatrixMultiplication(const Matrix<T, n>& A, const Matrix<T, n>& B, Matrix<T, n>& C) {

    C.clear();
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            T A_ik = A[i][k];
            for (int j = 0; j < n; ++j)
                C[i][j] += A_ik * B[k][j];
        }
    }
}


template<typename T, int n>
void sumAndPrecopyMatrixMultiplication(const Matrix<T, n>& A, const Matrix<T, n>& B, Matrix<T, n>& C) {

    C.clear();
    T B_col_j[n] = {};
    for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k)
            B_col_j[k] = B[k][j];

        for (int i = 0; i < n; ++i) {
            T S = 0;
            for (int k = 0; k < n; ++k) {
                S += A[i][k] * B_col_j[k];
            }
            C[i][j] = S;
        }
    }
}


template<typename T, int n>
void naiveMatrixMultiplication_Parallel(const Matrix<T, n>& A, const Matrix<T, n>& B, Matrix<T, n>& C) {

    C.clear();
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}


template<typename T, int n>
void transposeMatrixMultiplication_Parallel(const Matrix<T, n>& A, const Matrix<T, n>& B, Matrix<T, n>& C) {

    C.clear();
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            T A_ik = A[i][k];
            for (int j = 0; j < n; ++j)
                C[i][j] += A_ik * B[k][j];
        }
    }
}


template<typename T, int n>
void sumAndPrecopyMatrixMultiplication_Parallel(const Matrix<T, n>& A, const Matrix<T, n>& B, Matrix<T, n>& C) {

    C.clear();
    #pragma omp parallel for
    for (int j = 0; j < n; ++j) {
        T B_col_j[n] = {};
        for (int k = 0; k < n; ++k)
            B_col_j[k] = B[k][j];

        for (int i = 0; i < n; ++i) {
            T S = 0;
            for (int k = 0; k < n; ++k) {
                S += A[i][k] * B_col_j[k];
            }
            C[i][j] = S;
        }
    }
}


template<typename T, int n>
void naiveMatrixMultiplication_ParallelCollapse(const Matrix<T, n>& A, const Matrix<T, n>& B, Matrix<T, n>& C) {

    C.clear();
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}



template<int n, typename T = int>
void testFunction(const int threadsNum = 4) {

    omp_set_num_threads(threadsNum);

    Matrix<T, n> A;
    Matrix<T, n> B;
    Matrix<T, n> C;

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

    printf("%d ", n);

    {
        auto start = std::chrono::high_resolution_clock::now();

        naiveMatrixMultiplication(A, B, C);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();

        printf("%lf ", microseconds / 1000000.0);
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        transposeMatrixMultiplication(A, B, C);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();

        printf("%lf ", microseconds / 1000000.0);
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        sumAndPrecopyMatrixMultiplication(A, B, C);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();
        
        printf("%lf ", microseconds / 1000000.0);
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        naiveMatrixMultiplication_Parallel(A, B, C);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();
        
        printf("%lf ", microseconds / 1000000.0);
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        transposeMatrixMultiplication_Parallel(A, B, C);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();

        printf("%lf ", microseconds / 1000000.0);
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        sumAndPrecopyMatrixMultiplication_Parallel(A, B, C);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();
        
        printf("%lf ", microseconds / 1000000.0);
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        naiveMatrixMultiplication_ParallelCollapse(A, B, C);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();
        
        printf("%lf\n", microseconds / 1000000.0);
    }
}

#endif
