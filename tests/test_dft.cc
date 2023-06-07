#include <iostream>
#include <complex>
#include <cmath>
#include <vector>

#include <random>

#include <tests/benchmark_timer.h>
#include <core/dft.h>

#include <core/parallel_dft.h>

#include <unistd.h>
#include <string>

template <typename T>
void PrintArray(T* array, size_t N) {
    for (size_t i=0; i<N; i++)
        std::cout << array[i] << " ";
    std::cout << std::endl;
}

void truncate(Complex *array, size_t N, FloatType limit) {
    for(size_t i=0; i<N; i++) {
        FloatType real, imag;
        if (std::abs(array[i].real()) < limit) {
            real = 0;
        } else {
            real = array[i].real();
        }
        if (std::abs(array[i].imag()) < limit) {
            imag = 0;
        } else {
            imag = array[i].imag();
        }
        array[i] = Complex(real, imag);
    }
}

bool checkIsClose(Complex *a, Complex *b, size_t N, FloatType eps=1E-3) {
    for (size_t i=0; i<N; i++) {
        if(norm(a[i] - b[i]) > eps)
            return false;
    }
    return true;
}

void TestSequentialFFT(size_t N) {
    constexpr FloatType max_val = 1000;

    std::vector<Complex> x(N);
    for (size_t i=0; i<N; i++) {
        x[i] = (rand() % 2*max_val) - max_val;
    }

    std::vector<Complex> d1(N), d2(N), d3(N);

    auto func1 = [&](){ base_dft::DFT(x.begin(), x.end(), d1.begin()); };
    auto func2 = [&](){ recursive_fft::DFT(x.begin(), x.end(), d2.begin());};
    auto func3 = [&](){ iterative_fft::DFT(x.begin(), x.end(), d3.begin());};

    timeFunction(func1, "Base");
    timeFunction(func2, "Recursive");
    timeFunction(func3, "Iterative");


    assert(checkIsClose(d1.data(), d2.data(), N));
    assert(checkIsClose(d2.data(), d3.data(), N));
    assert(checkIsClose(d3.data(), d1.data(), N));


    auto func4 = [&](){ base_dft::IDFT(d1.begin(), d1.end(), d1.begin()); };
    auto func5 = [&](){ recursive_fft::IDFT(d2.begin(), d2.end(), d2.begin());};
    auto func6 = [&](){ iterative_fft::IDFT(d3.begin(), d3.end(), d3.begin());};

    timeFunction(func4, "Inv Base");
    timeFunction(func5, "Inv Recursive");
    timeFunction(func6, "Inv Iterative");


    assert(checkIsClose(d1.data(), x.data(), N));
    assert(checkIsClose(d2.data(), x.data(), N));
    assert(checkIsClose(d3.data(), x.data(), N));
}

void CompareSequentialFFT(size_t N) {
    constexpr FloatType max_val = 1000;

    std::vector<Complex> x(N);
    for (size_t i=0; i<N; i++) {
        x[i] = (rand() % 2*max_val) - max_val;
    }

    std::vector<Complex> d1(N), d2(N), d3(N);

    auto func2 = [&](){ recursive_fft::DFT(x.begin(), x.end(), d2.begin());};
    auto func3 = [&](){ iterative_fft::DFT(x.begin(), x.end(), d3.begin());};

    timeFunction(func2, "Recursive");
    timeFunction(func3, "Iterative");
}

void TestParallelDFT(size_t N) {
    constexpr FloatType max_val = 1000;

    std::vector<Complex> x(N);
    for (size_t i=0; i<N; i++) {
        x[i] = (rand() % 2*max_val) - max_val;
    }

    std::vector<Complex> d01(N), d02(N), d03(N), d11(N), d12(N), d13(N), d_seq1(N), d_seq2(N), d_seq3(N);

    auto func_seq1 = [&](){ base_dft::DFT(x.begin(), x.end(), d_seq1.begin()); };
    auto func_seq2 = [&](){ recursive_fft::DFT(x.begin(), x.end(), d_seq2.begin()); };
    auto func_seq3 = [&](){ iterative_fft::DFT(x.begin(), x.end(), d_seq3.begin()); };

    auto func01 = [&](){ base_dft::ParallelDFT(x.begin(), x.end(), d01.begin(), FixedThreadsParallelizer()); };
    auto func02 = [&](){ recursive_fft::ParallelDFT(x.begin(), x.end(), d02.begin(), FixedThreadsParallelizer());};
    auto func03 = [&](){ iterative_fft::ParallelDFT(x.begin(), x.end(), d03.begin(), FixedThreadsParallelizer());};

    auto func11 = [&](){ base_dft::ParallelDFT(x.begin(), x.end(), d11.begin(), OmpParallelizer()); };
    auto func12 = [&](){ recursive_fft::ParallelDFT(x.begin(), x.end(), d12.begin(), OmpParallelizer());};
    auto func13 = [&](){ iterative_fft::ParallelDFT(x.begin(), x.end(), d13.begin(), OmpParallelizer());};

    timeFunction(func_seq1, "Base - sequential");
    timeFunction(func_seq2, "Recursive - sequential");
    timeFunction(func_seq3, "Iterative - sequential");

    timeFunction(func01, "Fixed Threads Base - parallel");
    timeFunction(func02, "Fixed Threads Recursive - parallel  ");
    timeFunction(func03, "Fixed Threads Iterative - parallel  ");

    timeFunction(func11, "Omp Base - parallel");
    timeFunction(func12, "Omp Recursive - parallel  ");
    timeFunction(func13, "Omp Iterative - parallel  ");


    assert(checkIsClose(d01.data(), d_seq3.data(), N));
    assert(checkIsClose(d02.data(), d_seq3.data(), N));
    assert(checkIsClose(d03.data(), d_seq3.data(), N));

    assert(checkIsClose(d11.data(), d_seq3.data(), N));
    assert(checkIsClose(d12.data(), d_seq3.data(), N));
    assert(checkIsClose(d13.data(), d_seq3.data(), N));
}

void CompareParallelDFT(size_t N) {
    constexpr FloatType max_val = 1000;

    std::vector<Complex> x(N);
    for (size_t i=0; i<N; i++) {
        x[i] = (rand() % 2*max_val) - max_val;
    }

    std::vector<Complex> d01(N), d02(N), d03(N), d11(N), d12(N), d13(N), d_seq1(N), d_seq2(N), d_seq3(N);

    auto func_seq2 = [&](){ recursive_fft::DFT(x.begin(), x.end(), d_seq2.begin()); };
    auto func_seq3 = [&](){ iterative_fft::DFT(x.begin(), x.end(), d_seq3.begin()); };

    auto func02 = [&](){ recursive_fft::ParallelDFT(x.begin(), x.end(), d02.begin(), FixedThreadsParallelizer());};
    auto func03 = [&](){ iterative_fft::ParallelDFT(x.begin(), x.end(), d03.begin(), FixedThreadsParallelizer());};

    auto func12 = [&](){ recursive_fft::ParallelDFT(x.begin(), x.end(), d12.begin(), OmpParallelizer());};
    auto func13 = [&](){ iterative_fft::ParallelDFT(x.begin(), x.end(), d13.begin(), OmpParallelizer());};

    timeFunction(func_seq2, "Recursive - sequential");
    timeFunction(func_seq3, "Iterative - sequential");

    timeFunction(func02, "Fixed Threads Recursive - parallel  ");
    timeFunction(func03, "Fixed Threads Iterative - parallel  ");

    timeFunction(func12, "Omp Recursive - parallel  ");
    timeFunction(func13, "Omp Iterative - parallel  ");
}

int main()
{

    std::string line(50, '-');
    std::cout << line << std::endl;

    std::cout << ">>> Input Size 2^6\n";
    TestParallelDFT(1 << 6);
    std::cout << line << std::endl;

    std::cout << ">>> Input Size 2^12\n";
    TestParallelDFT(1 << 12);
    std::cout << line << std::endl;

    std::cout << ">>> Input Size 2^18\n";
    CompareParallelDFT(1 << 18);
    std::cout << line << std::endl;
    
    std::cout << ">>> Input Size 2^24\n";
    CompareParallelDFT(1 << 24);
    std::cout << line << std::endl;

    return 0;
}