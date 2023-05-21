#pragma once

#include <core/dft.h>
#include <core/fft_types.h>

#include <number_theory/number_theory.h>
#include <core/modular_fft.h>

template < class T >
using PolynomialCoefficients = std::vector<T>;

template < class T >
static PolynomialCoefficients<T> MakeCoefficients(T value) {
    PolynomialCoefficients<T> out = {value};
    return out;
}

template < class T >
static PolynomialCoefficients<T> PopTraillingZeros(const PolynomialCoefficients<T> &coefs) {

}

// Immutable Polynomial
template < class T >
class Polynomial {
    const PolynomialCoefficients<T> coefs = {(T) 0};

public:
    Polynomial() = default;
    Polynomial(const PolynomialCoefficients<T> &coefs) : coefs(coefs) {}

    Polynomial(T value) : coefs(MakeCoefficients(value)) {}

    T operator[](size_t i) const {
        return (i >= coefs.size()) ? (T) 0 : coefs[i];
    }

    inline size_t Degree() const {
        return coefs.size() - 1;
    }

    Polynomial operator+(const Polynomial &other) const {
        PolynomialCoefficients<T> new_coefs(std::max(this->Degree(), other.Degree()));

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
        PolynomialCoefficients<T> new_coefs(std::max(this->Degree(), other.Degree()));

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
        PolynomialCoefficients<T> new_coefs(coefs);
        std::for_each(new_coefs.begin(), new_coefs.end(), [](T& x){ x = -x; });
        return Polynomial(new_coefs);
    }

    // Polynomial times scalar
    Polynomial operator*(T scalar) const {
        if (scalar == (T) 0) {
            return Polynomial();
        }
        PolynomialCoefficients<T> new_coefs(coefs);
        std::for_each(new_coefs.begin(), new_coefs.end(), [scalar](T& x){ x *= scalar; });
        return Polynomial(new_coefs);
    }

    // Polynomial over scalar
    Polynomial operator/(T scalar) const {
        if (scalar == (T) 0) {
            return Polynomial();
        }
        PolynomialCoefficients<T> new_coefs(coefs);
        std::for_each(new_coefs.begin(), new_coefs.end(), [scalar](T& x){ x /= scalar; });
        return Polynomial(new_coefs);
    }
};

template < class T >
Polynomial<T> operator*(T value, const Polynomial<T> &other) {
    if (value == (T) 0) {
        return Polynomial<T>();
    }
    PolynomialCoefficients<T> new_coefs(other.Degree() + 1);
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

template < class T1, class T2 >
Polynomial<Complex> ComplexMultiply(const Polynomial<T1> &A, const Polynomial<T2> &B) {
    const size_t degree_A = A.Degree();
    const size_t degree_B = B.Degree();
    
    // A * B has degree = degree_A + degree_B or 0 if one of the polynomials is 0.
    const size_t degree_product = degree_A + degree_B;
    
    // Next power of 2 after degree_product
    size_t N = (1 << (fft_utils::IntLog2(degree_product) + 1));

    std::vector<Complex> coefs_A(N);
    std::vector<Complex> coefs_B(N);

    for (size_t k = 0; k < N; k++) {
        coefs_A[k] = A[k];
        coefs_B[k] = B[k];
    }

    std::vector<Complex> evals_A(N);
    std::vector<Complex> evals_B(N);

    // Fourier transform
    IterativeFft::Transform(&coefs_A[0], &evals_A[0], N);
    IterativeFft::Transform(&coefs_B[0], &evals_B[0], N);

    // Multiply in frequency domain
    std::vector<Complex> evals_AB(N);
    std::transform(evals_A.begin(), evals_A.end(), evals_B.begin(), evals_AB.begin(), 
                    [](Complex a, Complex b){ return a * b; });
    
    // Inverse transform
    std::vector<Complex> coefs_AB(N);
    IterativeFft::InverseTransform(&evals_AB[0], &coefs_AB[0], N);
    
    // Only keep the first deg_A + deg_B coefficients
    coefs_AB.erase(std::next(coefs_AB.begin(), degree_product + 1), coefs_AB.end());

    return Polynomial<Complex>(coefs_AB);
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

Polynomial<int> IntegerMultiply(const Polynomial<int> &A, const Polynomial<int> &B) {
    const size_t degree_A = A.Degree();
    const size_t degree_B = B.Degree();
    
    // A * B has degree = degree_A + degree_B or 0 if one of the polynomials is 0.
    const size_t degree_product = degree_A + degree_B;
    
    // Next power of 2 after degree_product
    size_t N = (1 << (fft_utils::IntLog2(degree_product) + 1));

    nt::LL p = nt::FindPrimeInAP(N);
    nt::LL g = nt::PrimitiveRootModPrime(p);

    std::cout << "Using prime: " << p << std::endl;

    std::vector<nt::LL> coefs_A(N);
    std::vector<nt::LL> coefs_B(N);

    for (size_t k = 0; k < N; k++) {
        coefs_A[k] = A[k];
        coefs_B[k] = B[k];
    }

    std::vector<nt::LL> evals_A(N);
    std::vector<nt::LL> evals_B(N);

    ModularFftTransform(coefs_A.begin(), coefs_A.end(), evals_A.begin(), p, g);
    ModularFftTransform(coefs_B.begin(), coefs_B.end(), evals_B.begin(), p, g);


    // Multiply in frequency domain
    std::vector<nt::LL> evals_AB(N);
    std::transform(evals_A.begin(), evals_A.end(), evals_B.begin(), evals_AB.begin(), 
                    [p](nt::LL a, nt::LL b){ return (a * b) % p; });
    
    // Inverse transform
    std::vector<nt::LL> coefs_AB(N);
    ModularFftInverseTransform(evals_AB.begin(), evals_AB.end(), coefs_AB.begin(), p, g);

    // Bring ints to interval [-p/2, +p/2]
    std::vector<int> coefs_output(degree_product + 1);
    for (size_t k = 0; k <= degree_product; k++) {
        coefs_output[k] = (coefs_AB[k] <= p/2) ? coefs_AB[k] : coefs_AB[k] - p;
    }

    return Polynomial<int>(coefs_output);
}