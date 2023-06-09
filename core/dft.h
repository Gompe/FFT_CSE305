#ifndef DFT_H_
#define DFT_H_

#include <iostream>
#include <algorithm>
#include <cassert>

// Multithreading
#include <thread>
#include <future>

#include <core/fft_types.h>
#include <core/fft_utils.h>

#define RECURSIVE_FFT_BASE_CASE_SIZE 1

namespace dft_detail {
    template < class InputIt, class OutputIt >
    bool IsMemEqual (InputIt first, OutputIt d_first) {
        const void *mem_addr_src = (void *) (& (*first));
        const void *mem_addr_dst = (void *) (& (*d_first));

        return (mem_addr_src == mem_addr_dst);
    }

    template < class InputIt, class OutputIt, class Impl >
    struct Wrapper {

        void call(InputIt first, InputIt last, OutputIt d_first, bool is_inverse_transform) {

            const bool condition = IsMemEqual(first, d_first);
            const size_t N = std::distance(first, last);
            using ComplexType = typename std::iterator_traits<OutputIt>::value_type;

            if (!condition) {
                Impl{}.template operator()<InputIt, OutputIt>(first, last, d_first, 1, is_inverse_transform);
            }
            else {
                std::vector<ComplexType> storage(N);
                Impl{}.template operator()<InputIt, typename std::vector<ComplexType>::iterator>(first, last, storage.begin(), 1, is_inverse_transform);
                std::copy(storage.begin(), storage.end(), d_first);
            }

            if (is_inverse_transform) {
                std::for_each(d_first, d_first + N, [N](ComplexType& value){ value /= (ComplexType) N; });
            }
        }
    };
}

// Implements O(N^2) dft
namespace base_dft {

    struct ImplDFT {
        template < class InputIt, class OutputIt >
        void operator()(InputIt first, InputIt last, OutputIt d_first, size_t stride,
                        bool is_inverse_transform) {
            using ComplexType = typename std::iterator_traits<OutputIt>::value_type; 

            const size_t N = std::distance(first, last);
            const size_t n = 1 + (N - 1)/stride;

            for (size_t k = 0; k < n; k++) {
                ComplexType twiddle = (ComplexType) fft_utils::RootOfUnity(n, is_inverse_transform ? k : -((int) k));

                d_first[k] = (ComplexType) 0;
                ComplexType twiddle_factor = (ComplexType) 1; 

                for (size_t m = 0; m < N; m += stride) {
                    d_first[k] += first[m] * twiddle_factor;
                    
                    // twiddle_factor = root_of_unity^(k * (m+1)) afte the following line
                    twiddle_factor *= twiddle;
                }
            }
        }
    };

    template < class InputIt, class OutputIt >
    void DFT(InputIt first, InputIt last, OutputIt d_first) {
        dft_detail::Wrapper<InputIt, OutputIt, ImplDFT> wrapper;
        wrapper.call(first, last, d_first, false);
    }

    template < class InputIt, class OutputIt >
    void IDFT(InputIt first, InputIt last, OutputIt d_first) {
        dft_detail::Wrapper<InputIt, OutputIt, ImplDFT> wrapper;
        wrapper.call(first, last, d_first, true);
    } 
}; // namespace base_dft

namespace recursive_fft {
    // Implementation inspired by pseudocode in
    // https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm

    struct ImplDFT {
        template < class InputIt, class OutputIt >
        void operator()(InputIt first, InputIt last, OutputIt d_first, size_t stride,
                        bool is_inverse_transform) {
            
            using ComplexType = typename std::iterator_traits<OutputIt>::value_type; 

            const size_t N = std::distance(first, last);
            const size_t n = 1 + (N - 1)/stride;

            if (n <= RECURSIVE_FFT_BASE_CASE_SIZE) {
                base_dft::ImplDFT{}.template operator()<InputIt, OutputIt>(first, last, d_first, stride, is_inverse_transform);
                return;
            }

            if (n % 2 != 0) {
                std::cout << n << std::endl;
            }
            assert(n % 2 == 0);

            // Even terms 0, 2*stride, 4*stride,...
            this->template operator()<InputIt, OutputIt>(first, last, d_first, 2*stride, is_inverse_transform);

            // Odd terms 1*stride, 3*stride, ...
            this->template operator()<InputIt, OutputIt>(first + stride, last, d_first + n/2, 2*stride, is_inverse_transform);

            // exp(2PI i/n) if IDFT
            // else exp(-2PI i/n)
            ComplexType root_of_unity = (ComplexType) fft_utils::RootOfUnity(n, is_inverse_transform ? 1 : -1);
            ComplexType twiddle = (ComplexType) 1;

            for (size_t k = 0; k < n/2; k++) {
                ComplexType p = d_first[k];
                ComplexType q = twiddle * d_first[k + n/2];
                d_first[k] = p + q;
                d_first[k + n/2] = p - q;

                // twiddle = root_of_unity^(k+1) after next line
                twiddle *= root_of_unity;
            }
        } 
    };

    template < class InputIt, class OutputIt >
    void DFT(InputIt first, InputIt last, OutputIt d_first) {
        dft_detail::Wrapper<InputIt, OutputIt, ImplDFT> wrapper;
        wrapper.call(first, last, d_first, false);
    }

    template < class InputIt, class OutputIt >
    void IDFT(InputIt first, InputIt last, OutputIt d_first) {
        dft_detail::Wrapper<InputIt, OutputIt, ImplDFT> wrapper;
        wrapper.call(first, last, d_first, true);
    } 
}; // namespace recursive_fft

namespace iterative_fft {
    // Implementation inspired by the code at 
    // https://github.com/roguh/cuda-fft/blob/069554b979d9ce82257bf3d3efa2a386d3abc2a1/main.cu#L247
    
    struct ImplDFT {
        template < class InputIt, class OutputIt >
        void operator()(InputIt first, InputIt last, OutputIt d_first, size_t stride,
                        bool is_inverse_transform) {
            
            using ComplexType = typename std::iterator_traits<OutputIt>::value_type; 

            const size_t N = std::distance(first, last);

            const size_t n = 1 + (N - 1)/stride;
            const size_t logn = fft_utils::IntLog2(n);

            for (size_t i = 0; i < n; i++) {
                // Reverse the logN bits in the index.
                size_t rev = fft_utils::ReverseBits(i, logn);

                // Base case: set the output to the bit-reversed input.
                d_first[i] = first[stride * rev];
            }

            int transform_sign = is_inverse_transform ? -1 : 1;
            // Loop over 2, 4, 8, 16, ..., n
            for (size_t s = 1; s <= logn; s++) {
                // twiddle = exp(-2 PI i / 2^s)
                ComplexType twiddle = (ComplexType) fft_utils::RootOfUnity(fft_utils::PowerOfTwo(s), -transform_sign);

                // Iterate through out in strides of length m=2**s
                // Set k to 0, 2^s, 2 * 2^s, 3 * 2^s, ..., N-2^s
                for (size_t k = 0; k < n; k += fft_utils::PowerOfTwo(s)) {
                    ComplexType twiddle_factor = (ComplexType) 1;

                    for (size_t j = 0; j < (size_t) fft_utils::PowerOfTwo(s-1); j++) {
                        ComplexType a = d_first[k + j];
                        ComplexType b = twiddle_factor * d_first[k + j + fft_utils::PowerOfTwo(s-1)];
                        d_first[k + j] = a + b;
                        d_first[k + j + fft_utils::PowerOfTwo(s-1)] = a - b;

                        twiddle_factor *= twiddle;
                    }
                }
            }
        }
    }; 

    template < class InputIt, class OutputIt >
    void DFT(InputIt first, InputIt last, OutputIt d_first) {
        dft_detail::Wrapper<InputIt, OutputIt, ImplDFT> wrapper;
        wrapper.call(first, last, d_first, false);
    }

    template < class InputIt, class OutputIt >
    void IDFT(InputIt first, InputIt last, OutputIt d_first) {
        dft_detail::Wrapper<InputIt, OutputIt, ImplDFT> wrapper;
        wrapper.call(first, last, d_first, true);
    } 
}; // namespace iterative_fft

#endif