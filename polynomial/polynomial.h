#pragma once

#ifndef POLYNOMIAL_POLYNOMIAL_H
#define POLYNOMIAL_POLYNOMIAL_H

#include <core/dft.h>
#include <core/fft_types.h>

#include <number_theory/number_theory.h>
#include <core/modular_fft.h>

#include <core/parallel.h>
#include <numeric>

template < class T >
using PolynomialCoefficients = std::vector<T>;

template < class T >
static std::vector<T> MakeCoefficients(T value) {
    std::vector<T> out = {value};
    return out;
}

constexpr size_t LIMIT_NAIVE_MULTIPLY = 8;

// Immutable Polynomial
template < class T >
class Polynomial {
private:
    std::vector<T> m_coefs = {(T) 0};

    void FixSize() {
        size_t ptr = m_coefs.size() - 1;
        while (ptr > 0 && m_coefs[ptr] == static_cast<T>(0)) {
            ptr--;
        }

        m_coefs.erase(m_coefs.begin() + ptr + 1, m_coefs.end());

        // while ( (m_coefs.size() >= 2) && (m_coefs[m_coefs.size() - 1] == (T) 0) ) {
        //     m_coefs.pop_back();
        // }
    }

public:
    Polynomial() = default;
    Polynomial(const std::vector<T> &m_coefs) : m_coefs(m_coefs) {
        FixSize();
    }
    Polynomial(T value) {
        m_coefs = {(T) 0};
    }

    // Const iterators
    typename std::vector<T>::const_iterator ConstBegin() const {
        return m_coefs.begin();
    }

    typename std::vector<T>::const_iterator ConstEnd() const {
        return m_coefs.end();
    }

    // Getter
    T operator[](size_t i) const {
        return (i >= m_coefs.size()) ? (T) 0 : m_coefs[i];
    }

    // Setter
    void SetCoefficient(size_t i, T value) {
        m_coefs.resize(std::max(m_coefs.size(), i + 1));
        m_coefs[i] = value;
        FixSize();
    }

    inline size_t Degree() const {
        return m_coefs.size() - 1;
    }

    Polynomial operator+(const Polynomial &other) const {
        std::vector<T> new_coefs(std::max(this->Degree(), other.Degree()));

        // Adds two polynomials
        for (size_t i=0; i < new_coefs.size(); i++) {
            new_coefs[i] = this->operator[](i) + other[i];
        }

        // Pops trailling zeros
        while (new_coefs.size() > 1 && new_coefs[new_coefs.size() - 1] == 0) {
            new_coefs.pop_back();
        }

        return Polynomial(new_coefs);
    }

    Polynomial operator-(const Polynomial &other) const {
        // Identical to previous code but with minus operations
        std::vector<T> new_coefs(std::max(this->Degree(), other.Degree()));

        // Adds two polynomials
        for (size_t i=0; i < new_coefs.size(); i++) {
            new_coefs[i] = this->operator[](i) - other[i];
        }

        // Pops trailing zeros
        while (new_coefs.size() > 1 && new_coefs[new_coefs.size() - 1] == 0) {
            new_coefs.pop_back();
        }

        return Polynomial(new_coefs);
    }

    // Unary minus
    Polynomial operator-() const {
        std::vector<T> new_coefs(m_coefs);
        std::for_each(new_coefs.begin(), new_coefs.end(), [](T& x){ x = -x; });
        return Polynomial(new_coefs);
    }

    // Polynomial times scalar
    Polynomial operator*(T scalar) const {
        if (scalar == (T) 0) {
            return Polynomial();
        }
        std::vector<T> new_coefs(m_coefs);
        std::for_each(new_coefs.begin(), new_coefs.end(), [scalar](T& x){ x *= scalar; });
        return Polynomial(new_coefs);
    }

    // Polynomial over scalar
    Polynomial operator/(T scalar) const {
        if (scalar == (T) 0) {
            return Polynomial();
        }
        std::vector<T> new_coefs(m_coefs);
        std::for_each(new_coefs.begin(), new_coefs.end(), [scalar](T& x){ x /= scalar; });
        return Polynomial(new_coefs);
    }

    template < class T1, class T2 >
    friend Polynomial<Complex> ComplexMultiply(const Polynomial<T1> &A, const Polynomial<T2> &B);
};

template <typename TypeFrom, typename TypeTo>
Polynomial<TypeTo> CastPolynomial(const Polynomial<TypeFrom> &P) {
    return Polynomial<TypeTo>(std::vector<TypeTo>(P.ConstBegin(), P.ConstEnd()));

} 

template < class T >
Polynomial<T> operator*(T value, const Polynomial<T> &other) {
    if (value == (T) 0) {
        return Polynomial<T>();
    }
    std::vector<T> new_coefs(other.Degree() + 1);
    for (int i = 0; i < other.Degree(); i++) {
        new_coefs[i] = value * other[i];
    }

    return Polynomial<T>(new_coefs);
}

// Multiplication
template < class T1, class T2 >
Polynomial<Complex> NaiveMultiply(const Polynomial<T1> &A, const Polynomial<T2> &B) {
    const size_t degree_A = A.Degree();
    const size_t degree_B = B.Degree();
    
    // A * B has degree = degree_A + degree_B or 0 if one of the polynomials is 0.
    const size_t degree_product = degree_A + degree_B;

    std::vector<Complex> coefs_AB(degree_product + 1);

    for (size_t k=0; k <= degree_product; k++) {
        coefs_AB[k] = 0;
        for (size_t l = 0; l <= k; l++) {
            coefs_AB[k] += A[l] * B[k-l];
        }
    }

    return Polynomial<Complex>(coefs_AB);
}

template <class T>
Polynomial<T> NaiveMultiply(const Polynomial<T> &A, const Polynomial<T> &B) {
    const size_t degree_A = A.Degree();
    const size_t degree_B = B.Degree();
    
    // A * B has degree = degree_A + degree_B or 0 if one of the polynomials is 0.
    const size_t degree_product = degree_A + degree_B;

    std::vector<T> coefs_AB(degree_product + 1);

    for (size_t k=0; k <= degree_product; k++) {
        coefs_AB[k] = 0;
        for (size_t l = 0; l <= k; l++) {
            coefs_AB[k] += A[l] * B[k-l];
        }
    }

    return Polynomial<T>(coefs_AB);
}

template < class T1, class T2 >
Polynomial<Complex> ComplexMultiply(const Polynomial<T1> &A, const Polynomial<T2> &B) {

    const size_t degree_A = A.Degree();
    const size_t degree_B = B.Degree();

    if (degree_A <= LIMIT_NAIVE_MULTIPLY || degree_B <= LIMIT_NAIVE_MULTIPLY) {
        return NaiveMultiply<T1, T2>(A, B);
    }

    // A * B has degree = degree_A + degree_B or 0 if one of the polynomials is 0.
    const size_t degree_product = degree_A + degree_B;
    
    // Next power of 2 after degree_product
    size_t N = (1 << (fft_utils::IntLog2(degree_product) + 1));

    std::vector<Complex> rep_A(N);
    std::vector<Complex> rep_B(N);

    // Sequential Version of the code:
    // std::fill(rep_A.begin(), rep_A.end(), (Complex) 0);
    // std::fill(rep_B.begin(), rep_B.end(), (Complex) 0);

    // std::copy(A.ConstBegin(), A.ConstEnd(), rep_A.begin());
    // std::copy(B.ConstBegin(), B.ConstEnd(), rep_B.begin());

    // // Perform DFT inplace -> rep_A now is the value representation of polynomial
    // iterative_fft::DFT(rep_A.begin(), rep_A.end(), rep_A.begin());
    // iterative_fft::DFT(rep_B.begin(), rep_B.end(), rep_B.begin());

    // Perform the 2 FFTs in parallel (~2x Faster)
    FixedThreadsParallelizer parallelizer(2);

    auto TransformA = [&](){
        std::fill(rep_A.begin(), rep_A.end(), (Complex) 0);
        std::copy(A.ConstBegin(), A.ConstEnd(), rep_A.begin());
        iterative_fft::DFT(rep_A.begin(), rep_A.end(), rep_A.begin());
    };

    auto TransformB = [&](){
        std::fill(rep_B.begin(), rep_B.end(), (Complex) 0);
        std::copy(B.ConstBegin(), B.ConstEnd(), rep_B.begin());
        iterative_fft::DFT(rep_B.begin(), rep_B.end(), rep_B.begin());
    };

    std::vector<std::function<void(void)>> tasks = {TransformA, TransformB};
    parallelizer.parallel_calls(tasks);

    // Multiply A * B in values domain
    std::vector<Complex> rep_AB(N);
    std::transform(rep_A.begin(), rep_A.end(), rep_B.begin(), rep_AB.begin(), 
                    [](Complex a, Complex b){ return a * b; });
    
    // Inverse transform
    iterative_fft::IDFT(rep_AB.begin(), rep_AB.end(), rep_AB.begin());
    
    // Only keep the first deg_A + deg_B coefficients
    rep_AB.erase(rep_AB.begin() + degree_product + 1, rep_AB.end());

    return Polynomial<Complex>(rep_AB);
}

template < class T1, class T2 >
Polynomial<FloatType> RealMultiply(const Polynomial<T1> &A, const Polynomial<T2> &B) {
    Polynomial<Complex> AB = ComplexMultiply(A, B);
    PolynomialCoefficients<FloatType> coefs_AB(AB.Degree() + 1);
    for (size_t k = 0; k <= AB.Degree(); k++) {
        coefs_AB[k] = AB[k].real();
    }

    return Polynomial<FloatType>(coefs_AB);
}

/// Multiplies A*B (mod p).
/// p must be a prime such that p === 1 (mod N) for N being the smallest power
/// of 2 larger than degree(A) + degree(B)
Polynomial<nt::Integer> ModularMultiply(const Polynomial<nt::Integer> &A, const Polynomial<nt::Integer> &B, const nt::Integer p) {

    if (A.Degree() <= LIMIT_NAIVE_MULTIPLY || B.Degree() <= LIMIT_NAIVE_MULTIPLY) {
        auto AB = NaiveMultiply<nt::Integer>(A, B);
        for (size_t k = 0; k <= AB.Degree(); k++) {
            AB.SetCoefficient(k, nt::SafeMod(AB[k], p));
        }
        return AB;
    }


    const size_t degree_output = A.Degree() + B.Degree();
    // Smallest power of 2 larger than the degree of the output
    const size_t N = fft_utils::PowerOfTwo(1 + fft_utils::IntLog2(degree_output));
    assert(p % N == 1);

    // Resize the coefficient vectors to have length N
    auto pad_coefs = [N](const auto &polynomial) {
        std::vector<nt::Integer> coefs(polynomial.ConstBegin(), polynomial.ConstEnd());
        coefs.resize(N);
        return coefs;
    };

    const std::vector<nt::Integer> coefs_A = pad_coefs(A);
    const std::vector<nt::Integer> coefs_B = pad_coefs(B);

    std::vector<nt::Integer> values_A(N), values_B(N);

    const nt::Integer g = nt::PrimitiveRootModPrime(p);

    // Perform the 2 FFTs in parallel
    FixedThreadsParallelizer parallelizer(2);

    auto TransformA = [&](){
        ModularFftTransform(coefs_A.begin(), coefs_A.end(), values_A.begin(), p, g);
    };

    auto TransformB = [&](){
        ModularFftTransform(coefs_B.begin(), coefs_B.end(), values_B.begin(), p, g);
    };

    // Evaluate Polynomials A and B at Nth roots of unity mod p
    std::vector<std::function<void(void)>> tasks = {TransformA, TransformB};
    parallelizer.parallel_calls(tasks);


    const auto mul = [p](nt::Integer a, nt::Integer b) {
        return (a * b) % p;
    };

    // Evaluate Polynomial AB at the same points (point-wise multiplication of values_A, values_B)
    std::vector<nt::Integer> values_AB(N);
    std::transform(values_A.begin(), values_A.end(), values_B.begin(), values_AB.begin(), mul);

    // Do Langrange interpolation to recover the coefficients of AB
    std::vector<nt::Integer> coefs_AB(N);
    ModularFftInverseTransform(values_AB.begin(), values_AB.end(), coefs_AB.begin(), p, g);

    return Polynomial<nt::Integer>(coefs_AB);
}


Polynomial<nt::Integer> IntegerMultiply(const Polynomial<nt::Integer> &A, const Polynomial<nt::Integer> &B) {
    
    if (A.Degree() <= LIMIT_NAIVE_MULTIPLY || B.Degree() <= LIMIT_NAIVE_MULTIPLY) {
        return NaiveMultiply<nt::Integer>(A, B);
    }

    const size_t degree_output = A.Degree() + B.Degree();

    // Large power of 2. It must be larger than the degree of the output.
    const size_t exponent = std::max(14, 1 + fft_utils::IntLog2(degree_output));
    const nt::Integer N = fft_utils::PowerOfTwo(exponent);

    // We will compute the product AB modulo 2 different primes
    size_t n_moduli = 2;
    Polynomial<nt::Integer> polynomials[n_moduli];

    const auto primes = nt::FindPrimesInAP(N, n_moduli);

    std::vector<std::function<void(void)>> tasks(n_moduli);

    for (size_t i = 0; i < n_moduli; i++) {
        tasks[i] = [&polynomials, &A, &B, &primes, i](){
            polynomials[i] = ModularMultiply(A, B, primes[i]);
        };
        // polynomials[i] = ModularMultiply(A, B, primes[i]);
    }

    FixedThreadsParallelizer parallelizer{};
    parallelizer.parallel_calls(tasks);

    // Now we recover the int coefficients from the CRT
    std::vector<nt::Integer> out_coefficients(degree_output + 1);
    auto find_coefficient = [&](int k) {
        std::vector<nt::Integer> remainders(n_moduli);
        for (size_t i = 0; i < n_moduli; i++) {
            remainders[i] = polynomials[i][k];
        } 

        out_coefficients[k] = nt::ChineseRemainderTheorem(remainders, primes);
    };

    parallelizer.parallel_for(0, degree_output + 1, find_coefficient);

    // Sequential Version
    // for (size_t k = 0; k <= degree_output; k++) {
    //     std::vector<nt::Integer> remainders(n_moduli);
    //     for (size_t i = 0; i < n_moduli; i++) {
    //         remainders[i] = polynomials[i][k];
    //     } 

    //     out_coefficients[k] = nt::ChineseRemainderTheorem(remainders, primes);
    // }

    // Make all coefficients to have the smallest possible norm while keeping
    // their modulo wrt each prime
    nt::Integer prime_product = std::accumulate(primes.begin(), primes.end(), (nt::Integer) 1,
        [](nt::Integer acc, nt::Integer p) -> nt::Integer { return acc * p; });

    std::for_each(out_coefficients.begin(), out_coefficients.end(), 
        [prime_product](nt::Integer &c) { c = (c <= prime_product/2) ? c : c - prime_product; });

    return Polynomial<nt::Integer>(out_coefficients);
}

Polynomial<int> IntegerMultiply(const Polynomial<int> &A, const Polynomial<int> &B) {    

    if (A.Degree() <= LIMIT_NAIVE_MULTIPLY || B.Degree() <= LIMIT_NAIVE_MULTIPLY) {
        return NaiveMultiply<int>(A, B);
    }

    const std::vector<nt::Integer> coefs_A(A.ConstBegin(), A.ConstEnd());
    const std::vector<nt::Integer> coefs_B(B.ConstBegin(), B.ConstEnd());

    const auto AB = IntegerMultiply(Polynomial<nt::Integer>(coefs_A), Polynomial<nt::Integer>(coefs_B));

    // May lose precision. It is ok.
    const std::vector<int> coefs_AB(AB.ConstBegin(), AB.ConstEnd());

    return Polynomial<int>(coefs_AB);
}

#endif