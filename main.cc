#include <iostream>
#include <complex>
#include <cmath>
#include <vector>

#include <random>

#include <tests/benchmark_timer.h>
#include <core/dft.h>

#include <core/parallel_dft.h>

#include <unistd.h>

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

void TestSequentialFFT() {
    constexpr FloatType max_val = 1000;
    constexpr size_t N = 1 << 12;

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

void TestParallelFor() {

    const size_t N = 1 << 18;

    std::atomic<int> par_sum(0);
    int seq_sum = 0;

    std::vector<int> array(N);

    for (size_t i=0; i < N; i++) {
        array[i] = rand() % 2;
    }

    auto parallel_sum = [&]() {
        auto foo = [&](int i){
            par_sum.fetch_add(array[i]);
        };

        FixedThreadsParallelizer parallelizer{6};
        parallelizer.parallel_for(0, N, foo);
    };

    auto sequential_sum = [&]() {
        for (auto x : array) {
            seq_sum += x;
        }
    };

    timeFunction(parallel_sum, "Parallel sum");
    timeFunction(sequential_sum, "Sequential sum");

    std::cout << "Parallel sum: " << par_sum << std::endl;
    std::cout << "Sequential sum: " << seq_sum << std::endl;
}

void TestParallelCalls() {
    const size_t N = 1 << 18;
    const size_t num_funcs = 100;

    std::atomic<int> par_sum(0);
    int seq_sum = 0;

    std::vector<int> array(N);

    for (size_t i=0; i < N; i++) {
        array[i] = rand() % 2;
    }

    std::vector<std::function<void(void)>> funcs(num_funcs);

    size_t idx_begin = 0;
    for (size_t i = 0; i < num_funcs; i++) {
        const size_t idx_end = std::min(N, idx_begin + N/num_funcs);
        funcs[i] = [&array, &par_sum, idx_begin, idx_end](){
            for (size_t j = idx_begin; j < idx_end; j++) {
                par_sum.fetch_add(array[j]);
            }
        };
        idx_begin = idx_end;
    }

    auto parallel_sum = [&]() {
        FixedThreadsParallelizer parallelizer{6};
        parallelizer.parallel_calls(funcs);
    };

    auto sequential_sum = [&]() {
        for (auto x : array) {
            seq_sum += x;
        }
    };

    timeFunction(parallel_sum, "Parallel sum");
    timeFunction(sequential_sum, "Sequential sum");

    std::cout << "Parallel sum: " << par_sum << std::endl;
    std::cout << "Sequential sum: " << seq_sum << std::endl;
}

template <typename Parallelizer>
void TestParallelDFT(size_t N, const Parallelizer &parallelizer) {
    constexpr FloatType max_val = 1000;

    std::vector<Complex> x(N);
    for (size_t i=0; i<N; i++) {
        x[i] = (rand() % 2*max_val) - max_val;
    }

    std::vector<Complex> d1(N), d2(N), d3(N), d_seq1(N), d_seq2(N), d_seq3(N);

    // auto func_seq1 = [&](){ base_dft::DFT(x.begin(), x.end(), d_seq1.begin()); };
    auto func_seq2 = [&](){ recursive_fft::DFT(x.begin(), x.end(), d_seq2.begin()); };
    auto func_seq3 = [&](){ iterative_fft::DFT(x.begin(), x.end(), d_seq3.begin()); };

    // auto func1 = [&](){ base_dft::ParallelDFT(x.begin(), x.end(), d1.begin(), parallelizer); };
    auto func2 = [&](){ recursive_fft::ParallelDFT(x.begin(), x.end(), d2.begin(), parallelizer);};
    auto func3 = [&](){ iterative_fft::ParallelDFT(x.begin(), x.end(), d3.begin(), parallelizer);};

    // timeFunction(func_seq1, "Base - sequential");
    timeFunction(func_seq2, "Recursive - sequential");
    timeFunction(func_seq3, "Iterative - sequential");

    // timeFunction(func1, "Base - parallel");
    timeFunction(func2, "Recursive - parallel  ");
    timeFunction(func3, "Iterative - parallel  ");


    // assert(checkIsClose(d1.data(), d_seq3.data(), N));
    assert(checkIsClose(d2.data(), d_seq3.data(), N));
    assert(checkIsClose(d3.data(), d_seq3.data(), N));
}

int main()
{
    size_t N = 1 << 24;
    size_t num_threads = 4;

    std::cout << "Fixed Threads Parallelizer.\n";
    FixedThreadsParallelizer parallelizer{num_threads};
    TestParallelDFT(N, parallelizer);

    std::cout << "OMP Parallelizer.\n";
    OmpParallelizer omp_parallelizer{};
    TestParallelDFT(N, omp_parallelizer);

    // std::cout << "Parallel For" << std::endl;
    // TestParallelFor();

    // std::cout << "Parallel Calls" << std::endl;
    // TestParallelCalls();
    // TestSequentialFFT();


    return 0;
}