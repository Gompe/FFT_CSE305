#pragma once

#include <number_theory/number_theory.h>

#include <core/fft_types.h>
#include <core/fft_utils.h>

#include <cassert>


template < class InputIt, class OutputIt >
static void ImplModularFft(InputIt first, InputIt last, OutputIt d_first,
                            nt::Integer p, nt::Integer g, bool is_inverse_transform) {
    
    const int N = std::distance(first, last);
    const int logN = fft_utils::IntLog2(N);
    const int k = (p - 1) / N;

    // For debugging
    assert(N == (1 << logN));
    assert(p % N == 1);
    assert(nt::IsPrime(p));

    if (is_inverse_transform) {
        // Multiplicative Inverse of g mod p
        g = nt::ModularExponentiation(g, p-2, p);
    }

    // Stores the bit reversal permutation of [first...last] in d_first
    fft_utils::BitReversalPermutation(first, last, d_first);

    // Take mod p for all elements in d_first
    for (int i = 1; i < N; i++) {
        d_first[i] = (d_first[i]) % p;
    }

    // omega^N === 1 (mod p) since k*N == p-1
    nt::Integer omega = nt::ModularExponentiation(g, k, p);

    for (int s = 1; s <= logN; s++) {
        // twiddle is a 2^s root of unity mod p.
        nt::Integer twiddle = nt::ModularExponentiation(omega, fft_utils::PowerOfTwo(logN - s), p);

        for (int k = 0; k < N; k += fft_utils::PowerOfTwo(s)) {
            nt::Integer twiddle_factor = 1;

            // Set both halves of the out array at the same time
            // j = 1, 4, 8, 16, ..., N / 2
            for (int j = 0; j < fft_utils::PowerOfTwo(s-1); j++) {
                nt::Integer a = d_first[k + j];
                nt::Integer b = (twiddle_factor * d_first[k + j + fft_utils::PowerOfTwo(s-1)]) % p;

                // Compute pow(twiddle, j)
                twiddle_factor = (twiddle_factor * twiddle) % p;

                d_first[k + j] = (a + b) % p;
                d_first[k + j + fft_utils::PowerOfTwo(s-1)] = (a - b) % p;
            }
        }
    }

    // C++ may use negative remainders. This brings all results to the range
    // [0...p-1]
    for (int i = 1; i < N; i++) {
        d_first[i] = (d_first[i] >= 0) ? d_first[i] : (d_first[i] + p);
    }
}

template < class InputIt, class OutputIt >
void ModularFftTransform(InputIt first, InputIt last, OutputIt d_first, nt::Integer p,
                         nt::Integer g) {
    ImplModularFft(first, last, d_first, p, g, false);
}

template < class InputIt, class OutputIt >
void ModularFftInverseTransform(InputIt first, InputIt last, OutputIt d_first, nt::Integer p,
                         nt::Integer g) {

    ImplModularFft(first, last, d_first, p, g, true);

    // Divide Output by N (modulo p)
    const int N = std::distance(first, last);
    // Multiplicative inverse of N modulo p
    const int inv_N = nt::ModularExponentiation(N, p-2, p);

    for (int i=0; i<N; i++) {
        d_first[i] = (d_first[i] * inv_N) % p;
    }
}