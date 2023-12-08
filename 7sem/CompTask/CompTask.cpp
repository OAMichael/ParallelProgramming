#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <omp.h>


constexpr const size_t m = 12;
constexpr const size_t N = (1 << m) - 1;
constexpr const size_t numIterations = 5000;
constexpr const size_t bVarianceNum = 11;
constexpr const double h = 1.0 / N;

constexpr const char *defaultOutFilename = "./test.dat";

/*
 * Solve system of linear algebraic equations:
 * 
 * Ax = f
 * 
 * where A is threediagonal matrix with diagonals a, b and c
 * 
 */
void solveThreediagonalSLAE(double a[N], double b[N], double c[N], double f[N], double x[N]) {

    #pragma omp parallel
    {
        // Forward pass
        uint64_t stride = 1;
        for(uint64_t nn = N, low = 2; nn > 1; nn /= 2, low *= 2, stride *= 2) {
            #pragma omp for
            for(int i = low - 1; i < N; i += stride * 2) {
                float alpha = -a[i] / b[i - stride];
                float gamma = -c[i] / b[i + stride];
                a[i] = alpha * a[i - stride];
                b[i] = alpha * c[i - stride] + b[i] + gamma * a[i + stride];
                c[i] = gamma * c[i + stride];
                f[i] = alpha * f[i - stride] + f[i] + gamma * f[i + stride];
            } 
        } 

        #pragma omp barrier

        // Reverse pass
        x[N / 2] = f[N / 2] / b[N / 2];
        for(stride /= 2; stride >= 1; stride /= 2) {
            #pragma omp for
            for(uint64_t i = stride - 1; i < N; i += stride * 2) {
                x[i] = (f[i]
                - (i - stride > 0 ? a[i] * x[i - stride] : 0.0)
                - (i + stride < N ? c[i] * x[i + stride] : 0.0)
                ) / b[i];
            }
        }
    }
}


double calculateResidue(const double yPrev[N], const double yCurr[N]) {
    double res = 0;
    for (size_t i = 0; i < N; ++i) {
        res += std::fabs(yPrev[i] - yCurr[i]);
    }
    return res;
}



int main(int argc, char *argv[]) {

    int num_threads = 1;
    double epsilon = 1e-32;

    std::ofstream outFile;
    std::string outFilename = defaultOutFilename;

    for (size_t i = 0; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-w" || arg == "--workers") {
            if (i + 1 < argc) {
                num_threads = std::atoi(argv[i + 1]);
                ++i;
                continue;
            }
        }
        if (arg == "-e" || arg == "--epsilon") {
            if (i + 1 < argc) {
                epsilon = std::atof(argv[i + 1]);
                ++i;
                continue;
            }
        }
        if (arg == "-o") {
            if (i + 1 < argc) {
                outFilename = argv[i + 1];
                ++i;
                continue;
            }
        }
    }

    omp_set_num_threads(num_threads);
    omp_set_dynamic(0);

    double *x = new double[N]();
    double *yPrev = new double[N]();
    double *yCurr = new double[N]();
    double *a = new double[N]();
    double *b = new double[N]();
    double *c = new double[N]();
    double *f = new double[N]();

    double **solution_of_b = new double*[bVarianceNum]();
    for (int sol = 0; sol < bVarianceNum; ++sol) {
        solution_of_b[sol] = new double[N]();
    }

    for (int sol = 0; sol < bVarianceNum; ++sol) {
        double rightBound = sol * 1.0 / (bVarianceNum - 1);

        // Initial approximation is: y = 1 + (b - 2) * x + x^2 
        for (size_t i = 0; i < N; ++i) {
            const double currX = i * h; 
            x[i] = currX;
            yPrev[i] = 1 + (rightBound - 2) * currX + currX * currX;
        }

        // Start time
        double start = omp_get_wtime();

        // Main iterative loop
        for (size_t iters = 0; iters < numIterations; ++iters) {
            constexpr const double h2_over_12 = h * h / 12.0;

            a[0] = 0;
            b[0] = 1;
            c[0] = 0;

            f[0] = 1;

            #pragma omp parallel for schedule(static, 3)
            for (size_t k = 1; k < N - 1; ++k) {
                const double exp_y_k_minus_1 = std::exp(yPrev[k - 1]);
                const double exp_y_k         = std::exp(yPrev[k]);
                const double exp_y_k_plus_1  = std::exp(yPrev[k + 1]);
            
                a[k] = 1.0 - h2_over_12 * exp_y_k_plus_1;
                b[k] = -2.0 - 10 * h2_over_12 * exp_y_k;
                c[k] = 1.0 - h2_over_12 * exp_y_k_minus_1;

                f[k] = h2_over_12 * (exp_y_k_plus_1 * (1.0 - yPrev[k + 1]) + 10.0 * exp_y_k * (1.0 - yPrev[k]) + exp_y_k_minus_1 * (1.0 - yPrev[k - 1]));
            }

            a[N - 1] = 0;
            b[N - 1] = 1;
            c[N - 1] = 0;

            f[N - 1] = rightBound;


            solveThreediagonalSLAE(a, b, c, f, yCurr);

            if (calculateResidue(yPrev, yCurr) < epsilon) {
                break;
            }

            std::memcpy(yPrev, yCurr, N * sizeof(double));
        }

        // End time
        double end = omp_get_wtime();

        std::memcpy(solution_of_b[sol], yCurr, N * sizeof(double));

        std::cout << "b = " << rightBound << std::endl;
        std::cout << "epsilon = " << std::scientific << epsilon << std::defaultfloat << std::endl;
        std::cout << "Number of executors: " << num_threads << std::endl;
        std::cout << "Elapsed time: " << end - start << " seconds\n" << std::endl;
    }


    outFile.open(outFilename);
    if (!outFile.is_open()) {
        std::cerr << "Could not open file: " << outFilename << std::endl;
    }
    else {
        for(int i = 0; i < N; ++i) {
            outFile << x[i];
            for (int sol = 0; sol < bVarianceNum; ++sol) {
                outFile << " " << solution_of_b[sol][i];
            }
            outFile << std::endl;
        }
        outFile.close();
    }
    
    delete[] x;
    delete[] yPrev;
    delete[] yCurr;
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] f;

    for (int sol = 0; sol < bVarianceNum; ++sol) {
        delete[] solution_of_b[sol];
    } 
    delete[] solution_of_b;

    return 0;
}
