#pragma once

#include <iostream>
#include <cmath>
#include <algorithm>

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

    inline Complex RootOfUnity(int N, int k) {
        // Returns exp(i*2*pi*k/N)
        FloatType theta = 2 * M_PI * k / (FloatType) N;
        return Complex(cos(theta), sin(theta));
    }
};