#pragma once

#include <iostream>
#include <cmath>
#include <algorithm>
#include <cassert>

#include <core/fft_types.h>

namespace fft_utils {
    inline bool IsPowerOfTwo(size_t N) {
        return ((bool) N) && !(bool) (N & (N-1));
    }

    inline int IntLog2(size_t N) {
        int output = -1;
        while (N != 0) {
            output++;
            N = N >> 1; 
        }
        return output;
    }

    inline int PowerOfTwo(size_t N) {
        return (1 << N);
    }

    // Outputs a number that corresponds to the first (from least significant to
    // most significant) n_bits of n reversed in binary.
    inline int ReverseBits(int n, int n_bits) {
        constexpr int max_n_bits = 31;

        int n_reversed = 0;
        for (int bit_idx=0; bit_idx < max_n_bits; bit_idx++) {
            n_reversed <<= 1;
            n_reversed |= (n & 1);
            n >>= 1;
        }

        n_reversed >>= (max_n_bits - n_bits);
        return n_reversed;
    }

    // BitReversalPermutation permutes the values in [first...last] using the
    // bit reversal permutation. It stores the output starting at d_first.
    //  It is allowed that d_first == first.
    // It is required that first, last, d_first are Random Access Iterators with
    // std::distance(first, last) being a power of 2.
    template < class InputIt, class OutputIt >
    inline void BitReversalPermutation(InputIt first, InputIt last, OutputIt d_first) {
        const int N = std::distance(first, last);
        const int logN = IntLog2(N);

        assert(IsPowerOfTwo(N));

        for (int i = 0; i < N; i++) {
            int j = fft_utils::ReverseBits(i, logN);

            // Here we want to make d_first[i] = first[j] and d_first[j] = first[i]
            // The following code ensures that there are no problems if first==d_first

            if (j < i) {
                // Already swapped i and j
                continue; 
            }

            if (d_first == first) {
                std::swap(d_first[i], d_first[j]);
            }
            else {
                d_first[i] = first[j];
                d_first[j] = first[i];
            }
        }
    }

    inline std::complex<double> RootOfUnity(int N, int k) {
        // Returns exp(i*2*pi*k/N)
        double theta = 2 * M_PI * k / (double) N;
        return std::complex<double> (cos(theta), sin(theta));
    }
};