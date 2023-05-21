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

/// O(N^2) DFT calculator
class BaseDft {
public:
    static void Transform(const Complex *x, Complex *out, int N, int stride=1);
    static void InverseTransform(const Complex *x, Complex *out, int N, int stride=1);

    static void ImplementationTransforms(
        const Complex *x, 
        Complex *out, 
        int N, 
        int stride,  
        bool is_inverse);
};

/// Recursive O(N log N) DFT calculator
// Uses the Radix-2 Algorithm. In particular, it only works for powers of 2.
class RecursiveFft {
    static int m_base_case_size;
public:
    static void Transform(const Complex *x, Complex *out, int N, int stride=1);
    static void InverseTransform(const Complex *x, Complex *out, int N, int stride=1);

    static void ImplementationTransforms(
        const Complex *x, 
        Complex *out, 
        int N, 
        int stride,  
        bool is_inverse);
};

/// Iterative O(N log N) DFT Calculator
class IterativeFft {
public:
    static void Transform(const Complex *x, Complex *out, int N, int stride=1);
    static void InverseTransform(const Complex *x, Complex *out, int N, int stride=1);

    static void ImplementationTransforms(
        const Complex *x, 
        Complex *out, 
        int N, 
        int stride,  
        bool is_inverse);
};

#endif