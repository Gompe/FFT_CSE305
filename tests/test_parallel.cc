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

int main() {
    TestParallelFor();
    TestParallelCalls();
}