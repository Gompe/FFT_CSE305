#include <core/parallel.h>

#include <iostream>
#include <vector>
#include <thread>

FixedThreadsParallelizer::FixedThreadsParallelizer() 
    : m_limit_thread_count(std::thread::hardware_concurrency()) {}

FixedThreadsParallelizer::FixedThreadsParallelizer(const size_t limit_thread_count)
    : m_limit_thread_count(limit_thread_count) {}


FixedThreadsParallelizer::ThreadGuard::ThreadGuard(const FixedThreadsParallelizer *instance) 
    : instance(instance) {
    
    // Acquire the maximum number of threads for the loop
    size_t curr_thread_count = instance->m_thread_count.load();;
    do {
        this->n_claimed_threads = instance->m_limit_thread_count - curr_thread_count;
    } while (!instance->m_thread_count.compare_exchange_weak(curr_thread_count, instance->m_limit_thread_count));
}

FixedThreadsParallelizer::ThreadGuard::~ThreadGuard() {
    // Acquire the maximum number of threads for the loop
    size_t curr_thread_count = instance->m_thread_count.load(); 
    size_t new_thread_count = curr_thread_count - this->n_claimed_threads;
    
    while (!instance->m_thread_count.compare_exchange_weak(curr_thread_count, new_thread_count)) {
        new_thread_count = curr_thread_count - this->n_claimed_threads;
    }
}

void FixedThreadsParallelizer::parallel_for(const int first, const int last, const std::function<void(int)> &func) const {
    const size_t length = std::max(last - first, 0);

    const ThreadGuard thread_guard(this); 

    const size_t num_threads = 1 + thread_guard.n_claimed_threads;

    const size_t chunk_size = length / num_threads;
    const size_t remainder = length % num_threads;

    std::vector<std::thread> workers(num_threads - 1);

    auto do_work = [&func](const int work_first, const int work_last) {
        for (int i = work_first; i < work_last; i++) {
            func(i);
        }
    };

    int work_first = first;
    for (size_t i = 0; i < num_threads - 1; i++) {
        const int extra_work = (remainder > i) ? 1 : 0;
        const int work_last = work_first + chunk_size + extra_work;

        workers[i] = std::thread(do_work, work_first, work_last);
        work_first = work_last;
    }

    // Last thread never does extra work
    do_work(work_first, last);

    for (size_t i = 0; i < num_threads - 1; i++) {
        workers[i].join();
    }
}

void FixedThreadsParallelizer::parallel_calls(std::vector< std::function<void(void)> > funcs) const {
    
    AtomicFIFO< std::function<void(void)> > fifo(funcs);
    const auto work_loop = [&](){
        while (true) {
            auto maybe_func = fifo.pop();
            if (maybe_func.has_value()) {
                auto func = maybe_func.value();
                func();
            }
            else {
                break;
            }
        }
    };

    const ThreadGuard thread_guard(this);
    const size_t num_threads = 1 + thread_guard.n_claimed_threads;

    std::vector<std::thread> workers(num_threads - 1);
    for (auto& worker : workers) {
        worker = std::thread(work_loop);
    }

    work_loop();

    for (auto& worker : workers) {
        worker.join();
    }
}

void OmpParallelizer::parallel_for(const int first, const int last, const std::function<void(int)> &func) const {
    #pragma omp parallel for
    for (int k = first; k < last; k++) {
        func(k);
    }
}

void OmpParallelizer::parallel_calls(std::vector< std::function<void(void)> > funcs) const {
    #pragma omp parallel for
    for (size_t k = 0; k < funcs.size(); k++) {
        funcs[k]();
    }
}