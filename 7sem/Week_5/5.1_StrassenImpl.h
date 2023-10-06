#ifndef STRASSEN_IMPL_H
#define STRASSEN_IMPL_H

#include <cstring>
#include <x86intrin.h>


#include "5.1_Strassen.h"


template<typename T>
struct genericVector2 {
    T x, y;

    genericVector2() {};
    genericVector2(const T& newX, const T& newY) : x{newX}, y{newY} {};

    inline genericVector2& operator+=(const genericVector2& other) {
        x += other.x;
        y += other.y;

        return *this;
    }

    inline genericVector2 operator+(const genericVector2& other) {
        return genericVector2(x + other.x, y + other.y);
    }

    inline genericVector2& operator=(const genericVector2& other) {
        x = other.x;
        y = other.y;

        return *this;
    }

};


template<typename T>
class Matrix {

private:
    T* data_;            // data stored in row-major order

    unsigned int size_;

public:
    inline unsigned int getSize()  const { return size_; }
    inline void clear() { std::memset(data_, 0, size_ * size_ * sizeof(T)); }

    inline const T* operator[](const int i) const { return data_ + size_ * i; }
    inline T* operator[](const int i) { return data_ + size_ * i; }

    inline MatrixView<T> getView() const {
        return MatrixView<T>(*this, {0, 0}, {size_, size_});
    }

    Matrix(const int size) {
        data_ = new T[size * size];
        size_ = size;
    }

    ~Matrix() {
        delete[] data_;
    }
};


template<typename T> 
class MatrixView {

private:
    const Matrix<T>* matrix_;

public:
    uvec2 start = uvec2{0, 0};
    uvec2 end = uvec2{0, 0};

    MatrixView() = delete;
    MatrixView(const Matrix<T>& mat) : matrix_{&mat} {};
    MatrixView(const Matrix<T>& mat, const uvec2& startView, const uvec2& endView)
    : matrix_{&mat},
    start{startView},
    end{endView}
    {};


    // To make view inside a view, recursive
    MatrixView(const MatrixView<T>& mat) {
        matrix_ = mat.getMatrix();
        start = mat.start;
        end = mat.end;
    };
    
    // To make view inside a view, recursive
    MatrixView(const MatrixView<T>& mat, const uvec2& startView, const uvec2& endView) {
        matrix_ = mat.getMatrix();

        start = mat.start;
        end = mat.start;

        start += startView;
        end += endView;
    };
    
    inline void makeView(const uvec2& startView, const uvec2& endView) {
        end = start + endView;
        start = start + startView;
    }

    inline const T* operator[](const int i) const { return matrix_->operator[](start.y + i) + start.x; }

    inline const Matrix<T>* getMatrix() const { return matrix_; }
};



template<typename T>
std::ostream& printMatrix(const Matrix<T>& mat, std::ostream& outFile) {
    int n = mat.getSize();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            outFile << mat[i][j] << " ";
        }
        outFile << std::endl;
    }

    return outFile;
}

template<typename T> 
std::ostream& printMatrix(const MatrixView<T>& matView, std::ostream& outFile) {
    int xSize = matView.end.x - matView.start.x;
    int ySize = matView.end.y - matView.start.y;
    
    for (int i = 0; i < ySize; ++i) {
        for (int j = 0; j < xSize; ++j) {
            outFile << matView[i][j] << " ";
        }
        outFile << std::endl;
    }

    return outFile;
}


// C = A + B
template<typename T>
void MatrixAdd(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C) {
    int n = C.getSize();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

// C = A + B
template<typename T>
void MatrixAdd(const MatrixView<T>& A, const MatrixView<T>& B, Matrix<T>& C) {
    int xSize = A.end.x - A.start.x;
    int ySize = A.end.y - A.start.y;

    for (int i = 0; i < ySize; ++i) {
        for (int j = 0; j < xSize; ++j) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

// C = A - B
template<typename T>
void MatrixSub(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C) {
    int n = C.getSize();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

// C = A - B
template<typename T>
void MatrixSub(const MatrixView<T>& A, const MatrixView<T>& B, Matrix<T>& C) {
    int xSize = A.end.x - A.start.x;
    int ySize = A.end.y - A.start.y;
    
    for (int i = 0; i < xSize; ++i) {
        for (int j = 0; j < ySize; ++j) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

// C = A * B
template<typename T>
void MatrixMul_Strassen(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C) {
	MatrixMul_Strassen(A.getView(), B.getView(), C, 0);
}


// C = A * B
template<typename T>
void MatrixMul_Strassen(const MatrixView<T>& A, const MatrixView<T>& B, Matrix<T>& C, int depth) {

    C.clear();
    const int n = C.getSize();
    if (n <= STRASSEN_SWITCH_MATRIX_SIZE) {
        MatrixMul_Blocks(A, B, C);
        return;
    }

    const int halfSize = n / 2;

    uvec2 viewVecStart{0, 0};
    uvec2 viewVecEnd{0, 0};

    
    MatrixView<T> A_11(A);
    MatrixView<T> A_12(A);
    MatrixView<T> A_21(A);
    MatrixView<T> A_22(A);
    MatrixView<T> B_11(B);
    MatrixView<T> B_12(B);
    MatrixView<T> B_21(B);
    MatrixView<T> B_22(B);


    viewVecStart.x = 0;
    viewVecStart.y = 0;
    viewVecEnd.x = halfSize;
    viewVecEnd.y = halfSize;

    A_11.makeView(viewVecStart, viewVecEnd);
    B_11.makeView(viewVecStart, viewVecEnd);


    viewVecStart.x = halfSize;
    viewVecStart.y = 0;
    viewVecEnd.x = n;
    viewVecEnd.y = halfSize;

    A_12.makeView(viewVecStart, viewVecEnd);
    B_12.makeView(viewVecStart, viewVecEnd);


    viewVecStart.x = 0;
    viewVecStart.y = halfSize;
    viewVecEnd.x = halfSize;
    viewVecEnd.y = n;

    A_21.makeView(viewVecStart, viewVecEnd);
    B_21.makeView(viewVecStart, viewVecEnd);


    viewVecStart.x = halfSize;
    viewVecStart.y = halfSize;
    viewVecEnd.x = n;
    viewVecEnd.y = n;

    A_22.makeView(viewVecStart, viewVecEnd);
    B_22.makeView(viewVecStart, viewVecEnd);


    Matrix<T> M_1{halfSize};
    Matrix<T> M_2{halfSize};
    Matrix<T> M_3{halfSize};
    Matrix<T> M_4{halfSize};
    Matrix<T> M_5{halfSize};
    Matrix<T> M_6{halfSize};
    Matrix<T> M_7{halfSize};

    // Little hack to force parallelism only once
    #pragma omp parallel if (depth < 1)
    {
        #pragma omp single
        {
            // ============= M_1 ============= //
            #pragma omp task
            {
                Matrix<T> tmpMat1{halfSize};
                Matrix<T> tmpMat2{halfSize};

                MatrixAdd(A_11, A_22, tmpMat1);
                MatrixAdd(B_11, B_22, tmpMat2);

                MatrixMul_Strassen(tmpMat1.getView(), tmpMat2.getView(), M_1, depth + 1);
            }
            

            // ============= M_2 ============= //
            #pragma omp task
            {
                Matrix<T> tmpMat{halfSize};

                MatrixAdd(A_21, A_22, tmpMat);
                
                MatrixMul_Strassen(tmpMat.getView(), B_11, M_2, depth + 1);
            }


            // ============= M_3 ============= //
            #pragma omp task
            {
                Matrix<T> tmpMat{halfSize};
                
                MatrixSub(B_12, B_22, tmpMat);

                MatrixMul_Strassen(A_11, tmpMat.getView(), M_3, depth + 1);
            }
            

            // ============= M_4 ============= //
            #pragma omp task
            {
                Matrix<T> tmpMat{halfSize};
                
                MatrixSub(B_21, B_11, tmpMat);
            
                MatrixMul_Strassen(A_22, tmpMat.getView(), M_4, depth + 1);
            }
            

            // ============= M_5 ============= //
            #pragma omp task
            {
                Matrix<T> tmpMat{halfSize};
                
                MatrixAdd(A_11, A_12, tmpMat);
            
                MatrixMul_Strassen(tmpMat.getView(), B_22, M_5, depth + 1);
            }


            // ============= M_6 ============= //
            #pragma omp task
            {
                Matrix<T> tmpMat1{halfSize};
                Matrix<T> tmpMat2{halfSize};
                
                MatrixSub(A_21, A_11, tmpMat1);
                MatrixAdd(B_11, B_12, tmpMat2);
            
                MatrixMul_Strassen(tmpMat1.getView(), tmpMat2.getView(), M_6, depth + 1);
            }


            // ============= M_7 ============= //
            #pragma omp task
            {
                Matrix<T> tmpMat1{halfSize};
                Matrix<T> tmpMat2{halfSize};
                
                MatrixSub(A_12, A_22, tmpMat1);
                MatrixAdd(B_21, B_22, tmpMat2);
            
                MatrixMul_Strassen(tmpMat1.getView(), tmpMat2.getView(), M_7, depth + 1);
            }
        }
    }

    // Merging results into final matrix
    for (int i = 0; i < halfSize; ++i) {
        for (int j = 0; j < halfSize; ++j) {
            int idxLUX = j + 0;
            int idxLUY = i + 0;

            int idxRUX = j + halfSize;
            int idxRUY = i + 0;

            int idxLDX = j + 0;
            int idxLDY = i + halfSize;

            int idxRDX = j + halfSize;
            int idxRDY = i + halfSize;


            C[idxLUY][idxLUX] = M_1[i][j] + M_4[i][j] - M_5[i][j] + M_7[i][j];
            C[idxRUY][idxRUX] = M_3[i][j] + M_5[i][j];
            C[idxLDY][idxLDX] = M_2[i][j] + M_4[i][j];
            C[idxRDY][idxRDX] = M_1[i][j] - M_2[i][j] + M_3[i][j] + M_6[i][j];
        }
    }
}


// C = A * B
template<typename T>
void MatrixMul_Naive(const MatrixView<T>& A, const MatrixView<T>& B, Matrix<T>& C) {
    int n = C.getSize();

    C.clear();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}


// C = A * B
template<typename T>
void MatrixMul_Naive(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C) {
    int n = C.getSize();

    C.clear();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}


// C = A * B
template<typename T>
void MatrixMul_SumPrecopy(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C) {
    int n = C.getSize();

    C.clear();
    T* B_col_j = new T[n];

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

    delete[] B_col_j;
}


// C = A * B
template<typename T>
void MatrixMul_SumPrecopy_Vector(const MatrixView<T>& A, const MatrixView<T>& B, Matrix<T>& C) {
    int n = C.getSize();
    C.clear();

    __m256 va, vb, vtemp;
    __m128 vlow, vhigh, vresult;

    T* B_col_j = new T[n];

    for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k)
            B_col_j[k] = B[k][j];

        for (int i = 0; i < n; ++i) {
            T S = 0;
            for (int k = 0; k < n; k += 8) {

                va = _mm256_loadu_ps(A[i] + k);
                vb = _mm256_loadu_ps(B_col_j + k);

                vtemp = _mm256_mul_ps(va, vb);

                vhigh = _mm256_extractf128_ps(vtemp, 1);

                vresult = _mm_add_ps(_mm256_castps256_ps128(vtemp), vhigh);
                vresult = _mm_hadd_ps(vresult, vresult);
                vresult = _mm_hadd_ps(vresult, vresult);

                S += _mm_cvtss_f32(vresult);
            }
            C[i][j] = S;
        }
    }

    delete[] B_col_j;
}


// C = C + A * B
template<typename T>
void MatrixMulAdd_SumPrecopy_Vector(const MatrixView<T>& A, const MatrixView<T>& B, Matrix<T>& C) {
    int n = C.getSize();

    __m256 va, vb, vtemp;
    __m128 vlow, vhigh, vresult;

    T* B_col_j = new T[n];

    for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k)
            B_col_j[k] = B[k][j];

        for (int i = 0; i < n; ++i) {
            T S = 0;
            for (int k = 0; k < n; k += 8) {

                va = _mm256_loadu_ps(A[i] + k);
                vb = _mm256_loadu_ps(B_col_j + k);

                vtemp = _mm256_mul_ps(va, vb);

                vhigh = _mm256_extractf128_ps(vtemp, 1);

                vresult = _mm_add_ps(_mm256_castps256_ps128(vtemp), vhigh);
                vresult = _mm_hadd_ps(vresult, vresult);
                vresult = _mm_hadd_ps(vresult, vresult);

                S += _mm_cvtss_f32(vresult);
            }
            C[i][j] += S;
        }
    }

    delete[] B_col_j;
}


// C = A * B
template<typename T>
void MatrixMul_Blocks(const MatrixView<T>& A, const MatrixView<T>& B, Matrix<T>& C) {
    const int N = A.end.x - A.start.x;
    const int effectiveBlock = std::min(N, MATRIX_BLOCK_SIZE);

    C.clear();
    #pragma omp parallel for
    for (unsigned int i = 0; i < N; i += effectiveBlock) {
        for (unsigned int j = 0; j < N; j += effectiveBlock) {
            Matrix<T> C_ij{effectiveBlock};
            C_ij.clear();

            for (unsigned int k = 0; k < N; k += effectiveBlock) {
                uvec2 A_view_start{k, i};
                uvec2 A_view_end{k + effectiveBlock, i + effectiveBlock};
                MatrixView<T> block_A_ik{A, A_view_start, A_view_end};
                                    
                uvec2 B_view_start{j, k};
                uvec2 B_view_end{j + effectiveBlock, k + effectiveBlock};
                MatrixView<T> block_B_kj{B, B_view_start, B_view_end};

                MatrixMulAdd_SumPrecopy_Vector(block_A_ik, block_B_kj, C_ij);
            }

            for (unsigned int l = 0; l < effectiveBlock; ++l) {
                for (unsigned int m = 0; m < effectiveBlock; ++m) {
                    C[i + l][j + m] = C_ij[l][m];
                }
            }
        }
    }
}

#endif 	// STRASSEN_IMPL_H
