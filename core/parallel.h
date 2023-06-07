#pragma once

#ifndef CORE_PARALLEL_H
#define CORE_PARALLEL_H

#include <core/dft.h>

#include <functional>
#include <thread>
#include <atomic>
#include <mutex>
#include <queue>
#include <optional>

/// parallel_for
/// Implements a parallel for loop. 
/// 
/// Calls the function func with all integers in [first, last) in parallel.
// void parallel_for(const int first, const int last, const std::function<void(int)> &func, const size_t num_threads);
// void parallel_for(const int first, const int last, const std::function<void(int)> &func);

// Very Simple Thread Safe Queue with mutex locks.
template <typename T>
class AtomicFIFO {
public:
    AtomicFIFO(const std::vector<T> &list) {
        for (auto value : list) {
            m_fifo.push(std::move(value));
        }
    }

    void push(const T& value) {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_fifo.push(value);
    }

    std::optional<T> pop() {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (!m_fifo.empty()) {
            T output = m_fifo.front();
            m_fifo.pop();
            return output;
        }
        return {};
    }

private:
    std::queue<T> m_fifo;
    mutable std::mutex m_mutex;
};

class FixedThreadsParallelizer {
public:
    FixedThreadsParallelizer();
    FixedThreadsParallelizer(const size_t limit_thread_count);

    void parallel_for(const int first, const int last, const std::function<void(int)> &func) const;

    void parallel_calls(std::vector< std::function<void(void)> > funcs) const;

private:
    size_t m_limit_thread_count;
    mutable std::atomic<size_t> m_thread_count{1};

    // RAII -
    // [1] Claim threads on constructor 
    // [2] Return threads to the Parallelizer on the destructor
    struct ThreadGuard {
        size_t n_claimed_threads = 0;
        const FixedThreadsParallelizer *instance;

        ThreadGuard() = delete;
        ThreadGuard(const ThreadGuard&) = delete;
        ThreadGuard(ThreadGuard&&) = delete;
        ThreadGuard& operator=(const ThreadGuard&) = delete;

        ThreadGuard(const FixedThreadsParallelizer *instance, const size_t max_num_threads = std::numeric_limits<size_t>::max());
        ~ThreadGuard();
    };    
};

class OmpParallelizer {
public:
    void parallel_for(const int first, const int last, const std::function<void(int)> &func) const;
    void parallel_calls(std::vector< std::function<void(void)> > funcs) const;
};

#endif