#include <iostream>
#include <complex>
#include <vector>
#include <fstream>
#include <chrono>

#include <x86intrin.h>

#include <omp.h>

//#define MAKE_TEST

#define NUM_THREADS 2
#define FFT_SWITCH_SIZE 16384


static constexpr const double PI = 3.141592653589793;


// ==================================== Simple sequential algorithm ==================================== //
void FFT_Sequential(std::vector<std::complex<double>>& a) {
    const size_t n = a.size();
    if (n == 1) {
        return;
    }
    const size_t halfSize = n / 2;

    std::vector<std::complex<double>> even(halfSize);
    std::vector<std::complex<double>> odd(halfSize);
    for (size_t i = 0; i < halfSize; ++i) {
        even[i] = a[2 * i];
        odd[i]  = a[2 * i + 1];
    }

    FFT_Sequential(even);
    FFT_Sequential(odd);

    const double ang = -2 * PI / n;
    std::complex<double> w(1.0);
    std::complex<double> wn(std::cos(ang), std::sin(ang));
    for (size_t i = 0; i < halfSize; ++i) {
        a[i] = even[i] + w * odd[i];
        a[i + halfSize] = even[i] - w * odd[i];
        w *= wn;
    }
}
// ===================================================================================================== //

// =========================== Sequential algorithm using vector instuctions =========================== //
void FFT_Sequential_Vector(std::vector<std::complex<double>>& a) {
    const static std::complex<double> imaginaryUnit = std::complex<double>(0.0, 1.0);

    const size_t n = a.size();
    if (n == 1) {
        return;
    }
    
    if (n == 2) {
        const std::complex<double> tmp0 = a[0];
        const std::complex<double> tmp1 = a[1];
        a[0] = tmp0 + tmp1;
        a[1] = tmp0 - tmp1;
        return;
    }
    
    if (n == 4) {
        const std::complex<double> tmp0 = a[0] + a[2];
        const std::complex<double> tmp1 = a[0] - a[2];
        const std::complex<double> tmp2 = a[1] + a[3];
        const std::complex<double> tmp3 = a[1] - a[3];

        a[0] = tmp0 + tmp2;
        a[2] = tmp0 - tmp2;
        a[1] = tmp1 - tmp3 * imaginaryUnit;
        a[3] = tmp1 + tmp3 * imaginaryUnit;
        return;
    }

    const size_t halfSize = n / 2;

    std::vector<std::complex<double>> even(halfSize);
    std::vector<std::complex<double>> odd(halfSize);
    for (size_t i = 0; i < halfSize; ++i) {
        even[i] = a[2 * i];
        odd[i]  = a[2 * i + 1];
    }

    FFT_Sequential_Vector(even);
    FFT_Sequential_Vector(odd);

    const double ang = -2 * PI / n;
    std::complex<double> exponent(std::cos(ang), std::sin(ang));

    // Some baddass is going to happen because of complex multiplication and vector instructions
    double* evenReal = new double[halfSize];
    double* evenImag = new double[halfSize];
    double* oddReal  = new double[halfSize];
    double* oddImag  = new double[halfSize];
    double* wnReal   = new double[halfSize];
    double* wnImag   = new double[halfSize];

    std::complex<double> tmp(1.0);
    for (size_t i = 0; i < halfSize; ++i) {
        wnReal[i] = tmp.real();
        wnImag[i] = tmp.imag();

        evenReal[i] = even[i].real();
        evenImag[i] = even[i].imag();
        
        oddReal[i] = odd[i].real();
        oddImag[i] = odd[i].imag();

        tmp *= exponent;
    } 

    for (size_t i = 0; i < halfSize; i += 4) {
        __m256d eR  = _mm256_loadu_pd(&evenReal[i]);
        __m256d eI  = _mm256_loadu_pd(&evenImag[i]);
        __m256d oR  = _mm256_loadu_pd(&oddReal[i]);
        __m256d oI  = _mm256_loadu_pd(&oddImag[i]);
        __m256d wnR = _mm256_loadu_pd(&wnReal[i]);
        __m256d wnI = _mm256_loadu_pd(&wnImag[i]);

        __m256d wn_o_Real = _mm256_sub_pd(_mm256_mul_pd(oR, wnR), _mm256_mul_pd(oI, wnI));
        __m256d wn_o_Imag = _mm256_add_pd(_mm256_mul_pd(oR, wnI), _mm256_mul_pd(oI, wnR));

        __m256d a_i_Real = _mm256_add_pd(eR, wn_o_Real);
        __m256d a_i_Imag = _mm256_add_pd(eI, wn_o_Imag);

        __m256d a_i_plus_half_Real = _mm256_sub_pd(eR, wn_o_Real);
        __m256d a_i_plus_half_Imag = _mm256_sub_pd(eI, wn_o_Imag);

        alignas(32) double a_i_reals[4];
        alignas(32) double a_i_imags[4];

        alignas(32) double a_i_plus_half_reals[4];
        alignas(32) double a_i_plus_half_imags[4];

        _mm256_storeu_pd(&a_i_reals[0], a_i_Real);
        _mm256_storeu_pd(&a_i_imags[0], a_i_Imag);
        _mm256_storeu_pd(&a_i_plus_half_reals[0], a_i_plus_half_Real);
        _mm256_storeu_pd(&a_i_plus_half_imags[0], a_i_plus_half_Imag);

        for (size_t k = 0; k < 4; ++k) {
            a[i + k] = std::complex<double>(a_i_reals[k], a_i_imags[k]);
            a[i + halfSize + k] = std::complex<double>(a_i_plus_half_reals[k], a_i_plus_half_imags[k]);
        }
    }

    delete[] evenReal;
    delete[] evenImag;
    delete[] oddReal;
    delete[] oddImag;
    delete[] wnReal;
    delete[] wnImag;
}
// ===================================================================================================== //

// ======================================== Parallel algorithm ========================================= //
void FFT_Parallel(std::vector<std::complex<double>>& a, int depth) {
    const size_t n = a.size();

    if (n < FFT_SWITCH_SIZE) {
        FFT_Sequential(a);
        return;
    }

    const size_t halfSize = n / 2;

    std::vector<std::complex<double>> even(halfSize);
    std::vector<std::complex<double>> odd(halfSize);

    #pragma omp parallel for
    for (size_t i = 0; i < halfSize; ++i) {
        even[i] = a[2 * i];
        odd[i]  = a[2 * i + 1];
    }

    // Little hack to force parallelism only once
    #pragma omp parallel if (depth < 1)
    {
        #pragma omp single
        {
            #pragma omp task
            {
                FFT_Parallel(even, depth + 1);
            }
            #pragma omp task
            {
                FFT_Parallel(odd, depth + 1);
            }
        }
    }

    const double ang = -2 * PI / n;
    std::complex<double> w(1.0);
    std::complex<double> wn(std::cos(ang), std::sin(ang));
    for (size_t i = 0; i < halfSize; ++i) {
        a[i] = even[i] + w * odd[i];
        a[i + halfSize] = even[i] - w * odd[i];
        w *= wn;
    }
}

void FFT_Parallel(std::vector<std::complex<double>>& a) {
    FFT_Parallel(a, 0);
}
// ===================================================================================================== //

// =========================== Parallel algorithm using vector instructions ============================ //
void FFT_Parallel_Vector(std::vector<std::complex<double>>& a, int depth) {
    const size_t n = a.size();

    if (n < FFT_SWITCH_SIZE) {
        FFT_Sequential_Vector(a);
        return;
    }

    const size_t halfSize = n / 2;

    std::vector<std::complex<double>> even(halfSize);
    std::vector<std::complex<double>> odd(halfSize);

    #pragma omp parallel for
    for (size_t i = 0; i < halfSize; ++i) {
        even[i] = a[2 * i];
        odd[i]  = a[2 * i + 1];
    }

    // Little hack to force parallelism only once
    #pragma omp parallel if (depth < 1)
    {
        #pragma omp single
        {
            #pragma omp task
            {
                FFT_Parallel_Vector(even, depth + 1);
            }
            #pragma omp task
            {
                FFT_Parallel_Vector(odd, depth + 1);
            }
        }
    }

    const double ang = -2 * PI / n;
    std::complex<double> w(1.0);
    std::complex<double> wn(std::cos(ang), std::sin(ang));
    for (size_t i = 0; i < halfSize; ++i) {
        a[i] = even[i] + w * odd[i];
        a[i + halfSize] = even[i] - w * odd[i];
        w *= wn;
    }
}

void FFT_Parallel_Vector(std::vector<std::complex<double>>& a) {
    FFT_Parallel_Vector(a, 0);
}
// ===================================================================================================== //





int main(int argc, char* argv[]) {
    if (argc < 3) {
        printf("Usage: %s <in_filename> <out_filename>\n", argv[0]);
        return 0;
    }
    
    std::ifstream inFile;
    inFile.open(argv[1]);
    if (!inFile.is_open()) {
        printf("Could not open file: %s\n", argv[1]);
        return -1;
    }

    int N = 0;
    inFile >> N;

    std::vector<std::complex<double>> signal;
    for (int i = 0; i < N; ++i) {
        double real, imag;

        inFile >> real;
        inFile >> imag;

        signal.emplace_back(real, imag);
    }
    inFile.close();

    omp_set_nested(1);
    omp_set_num_threads(NUM_THREADS);

#ifdef MAKE_TEST

    std::vector<std::complex<double>> signalCopies[4];
    for (auto& copy : signalCopies) {
        std::copy(signal.begin(), signal.end(), std::back_inserter(copy));
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        FFT_Sequential(signalCopies[0]);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();

        printf("Elapsed time for sequential:        %lf seconds\n", microseconds / 1000000.0);
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        FFT_Sequential_Vector(signalCopies[1]);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();

        printf("Elapsed time for sequential vector: %lf seconds\n", microseconds / 1000000.0);
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        FFT_Parallel(signalCopies[2]);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();

        printf("Elapsed time for parallel:          %lf seconds\n", microseconds / 1000000.0);
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        FFT_Parallel_Vector(signalCopies[3]);

        uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start
        ).count();

        printf("Elapsed time for parallel vector:   %lf seconds\n", microseconds / 1000000.0);
    }

#else
    auto start = std::chrono::high_resolution_clock::now();

    FFT_Parallel_Vector(signal);

    uint64_t microseconds = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::high_resolution_clock::now() - start
    ).count();

    printf("Elapsed time:        %lf seconds\n", microseconds / 1000000.0);


    std::ofstream outFile;
    outFile.open(argv[2]);
    if (!outFile.is_open()) {
        printf("Could not open file: %s\n", argv[2]);
        return -1;        
    }

    outFile << N << std::endl;
    for (int i = 0; i < N; ++i) {
        outFile << signal[i].real() << std::endl;
        outFile << signal[i].imag() << std::endl;
    }
    outFile.close();

#endif

    return 0;
}
