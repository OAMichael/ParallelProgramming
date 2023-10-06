#ifndef STRASSEN_H
#define STRASSEN_H

#include <fstream>
#include <iostream>

#include <omp.h>

#define NUM_THREADS 8

#define MATRIX_BLOCK_SIZE               64
#define STRASSEN_SWITCH_MATRIX_SIZE     4 * MATRIX_BLOCK_SIZE


template<typename T>
struct genericVector2;

using vec2 = genericVector2<float>;
using ivec2 = genericVector2<int>;
using uvec2 = genericVector2<unsigned int>;


template<typename T> 
class MatrixView;

template<typename T>
class Matrix;


template<typename T>
std::ostream& printMatrix(const Matrix<T>& mat, std::ostream& outFile = std::cout);

template<typename T> 
std::ostream& printMatrix(const MatrixView<T>& matView, std::ostream& outFile = std::cout);

// C = A + B
template<typename T>
void MatrixAdd(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C);

// C = A + B
template<typename T>
void MatrixAdd(const MatrixView<T>& A, const MatrixView<T>& B, Matrix<T>& C);

// C = A - B
template<typename T>
void MatrixSub(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C);

// C = A - B
template<typename T>
void MatrixSub(const MatrixView<T>& A, const MatrixView<T>& B, Matrix<T>& C);

// C = A * B
template<typename T>
void MatrixMul_Strassen(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C);

// C = A * B
template<typename T>
void MatrixMul_Strassen(const MatrixView<T>& A, const MatrixView<T>& B, Matrix<T>& C, int depth);

// C = A * B
template<typename T>
void MatrixMul_Naive(const MatrixView<T>& A, const MatrixView<T>& B, Matrix<T>& C);

// C = A * B
template<typename T>
void MatrixMul_Naive(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C);

// C = A * B
template<typename T>
void MatrixMul_SumPrecopy(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C);

// C = A * B
template<typename T>
void MatrixMul_SumPrecopy_Vector(const MatrixView<T>& A, const MatrixView<T>& B, Matrix<T>& C);

// C = C + A * B
template<typename T>
void MatrixMulAdd_SumPrecopy_Vector(const MatrixView<T>& A, const MatrixView<T>& B, Matrix<T>& C);

// C = A * B
template<typename T>
void MatrixMul_Blocks(const MatrixView<T>& A, const MatrixView<T>& B, Matrix<T>& C);


#endif  // STRASSEN_H
