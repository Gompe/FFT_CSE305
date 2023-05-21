/**
 * Implementation of the (non-concurrent/parallel) DFT algorithms. 
 * 
 * 
 * Remark. Our implementation supports values of N that fits in signed 32 bit
 * integers. That is: N <= 2,147,483,647.
*/

#include <core/dft.h>

#define DFT_CC_DEBUG 1

static constexpr int RECURSIVE_FFT_BASE_CASE_SIZE = 32;

///////////////////////////////////////////////////////////////////////////////
/// BaseDft Implementation O(N^2)

// Computes the Fourier Transform of x[0], x[s], x[2s], ..., x[(N-1)s]
// and stores it in out[0], ... , out[N-1]
void BaseDft::Transform(const Complex *x, Complex *out, int N, int stride){
    BaseDft::ImplementationTransforms(x, out, N, stride, false);
}

void BaseDft::InverseTransform(const Complex *x, Complex *out, int N, int stride){
    BaseDft::ImplementationTransforms(x, out, N, stride, true);
    // Divide each term by N
    std::for_each(out, out+N, [N](Complex &elem){elem /= (FloatType) N;});
}

void BaseDft::ImplementationTransforms(
        const Complex *x,
        Complex *out,
        int N,
        int stride,
        bool is_inverse) {
    
    int transform_sign = is_inverse ? -1 : 1;

    #ifdef OPEN_MP
        #pragma omp parallel for
    #endif
    for(int k=0; k<N; k++){
        out[k] = 0;
        for (int n=0; n<N; n++){
            out[k] += x[stride*n] * fft_utils::RootOfUnity(N, -transform_sign*k*n);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
/// Recursive FFT implementation O(N log N)

int RecursiveFft::m_base_case_size = RECURSIVE_FFT_BASE_CASE_SIZE;

void RecursiveFft::Transform(const Complex *x, Complex *out, int N, int stride) {
    if (!fft_utils::IsPowerOfTwo(N)) {
        throw std::logic_error("ERROR: In RecursiveFft::Transform: "
                                "radix-2 can only be used if N is a power of 2.\n");
    }

    RecursiveFft::ImplementationTransforms(x, out, N, stride, false);
}

void RecursiveFft::InverseTransform(const Complex *x, Complex *out, int N, int stride) {
    if (!fft_utils::IsPowerOfTwo(N)) {
        throw std::logic_error("ERROR: In RecursiveFft::InverseTransform: "
                                "radix-2 can only be used if N is a power of 2.\n");
    }

    RecursiveFft::ImplementationTransforms(x, out, N, stride, true);
    // Divide each term by N
    std::for_each(out, out+N, [N](Complex &elem){elem /= (FloatType) N;});
}

void RecursiveFft::ImplementationTransforms(
        const Complex *x,
        Complex *out,
        int N,
        int stride,
        bool is_inverse) {
    
    // Implementation inspired by pseudocode in
    // https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm

    if(N <= RecursiveFft::m_base_case_size) {
        BaseDft::ImplementationTransforms(x, out, N, stride, is_inverse);
        return;
    }

    #if DFT_CC_DEBUG
        assert(N%2 == 0);
    #endif

    RecursiveFft::ImplementationTransforms(x, out, N/2, 2*stride, is_inverse);
    RecursiveFft::ImplementationTransforms(x + stride, out + N/2, N/2, 2*stride, is_inverse);

    int transform_sign = is_inverse ? -1 : 1;
    for (int k=0; k < N/2; k++) {
        Complex p = out[k];
        Complex q = fft_utils::RootOfUnity(N, -transform_sign * k) * out[k + N/2];
        out[k] = p + q;
        out[k + N/2] = p-q;
    }
}

///////////////////////////////////////////////////////////////////////////////
/// Iterative FFT implementation O(N log N)

void IterativeFft::Transform(const Complex *x, Complex *out, int N, int stride)
{
    if (!fft_utils::IsPowerOfTwo(N)) {
        throw std::logic_error("ERROR: In IterativeFft::Transform: "
                                "this algorithm can only be used if N is a power of 2.\n");
    }

    IterativeFft::ImplementationTransforms(x, out, N, stride, false);
}

void IterativeFft::InverseTransform(const Complex *x, Complex *out, int N, int stride)
{
    if (!fft_utils::IsPowerOfTwo(N)) {
        throw std::logic_error("ERROR: In IterativeFft::Transform: "
                                "this algorithm can only be used if N is a power of 2.\n");
    }

    IterativeFft::ImplementationTransforms(x, out, N, stride, true);
    // Divide each term by N
    std::for_each(out, out+N, [N](Complex &elem){elem /= (FloatType) N;});
}

// Implementation similar to 
// https://github.com/roguh/cuda-fft/blob/069554b979d9ce82257bf3d3efa2a386d3abc2a1/main.cu#L247
void IterativeFft::ImplementationTransforms(
        const Complex *x, 
        Complex *out, 
        int N, 
        int stride,  
        bool is_inverse) {
    
    int logN = fft_utils::IntLog2(N); 

    for (int i = 0; i < N; i++) {
        // Reverse the logN bits in the index.
        int rev = fft_utils::ReverseBits(i, logN);

        // Base case: set the output to the bit-reversed input.
        out[i] = x[stride * rev];
    }

    int transform_sign = is_inverse ? -1 : 1;
    // Loop over 2, 4, 8, 16, ..., N
    for (int s = 1; s <= logN; s++) {
        // twiddle = exp(-2 PI i / 2^s)
        Complex twiddle = fft_utils::RootOfUnity(fft_utils::PowerOfTwo(s), -transform_sign);

        // Iterate through out in strides of length m=2**s
        // Set k to 0, 2^s, 2 * 2^s, 3 * 2^s, ..., N-2^s

        #ifdef OPEN_MP
            #pragma omp parallel for
        #endif
        for (int k = 0; k < N; k += fft_utils::PowerOfTwo(s)) {
            Complex twiddle_factor = 1;

            // Set both halves of the out array at the same time
            // j = 1, 4, 8, 16, ..., N / 2
            for (int j = 0; j < fft_utils::PowerOfTwo(s-1); j++) {
                Complex a = out[k + j];
                Complex b = twiddle_factor * out[k + j + fft_utils::PowerOfTwo(s-1)];

                // Compute pow(twiddle, j)
                twiddle_factor *= twiddle;

                out[k + j] = a + b;
                out[k + j + fft_utils::PowerOfTwo(s-1)] = a - b;
            }
        }
    }
}