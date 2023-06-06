#include <number_theory/number_theory.h>

#include <assert.h>

#include <iostream>
#include <limits.h>
#include <algorithm>
#include <numeric>
#include <cassert>

namespace nt {
    
Integer SafeMod(const Integer a, const Integer m) {
    return ((a % m) + m) % m;
}

/// Returns a^b (mod n) in time complexity O(log b)
/// The returned value is in the range [0...n-1]
Integer ModularExponentiation(Integer a, Integer b, Integer n) {
    Integer remainder = 1;
    a = a % n;

    while(b != 0) {
        if (b % 2 == 1) {
            remainder = (remainder * a) % n;
        }
        a = (a * a) % n;
        b = b/2;
    }

    return SafeMod(remainder, n);
}

/// Returns if n is Prime - Always returns the correct result
bool IsPrime(Integer n) {
    for (Integer d=2; d*d <= n; d++) {
        if (n % d == 0) {
            return false;
        }
    }
    return true;
}

/// Returns if n is probably prime. 
/// If it returns false, then n is composite for sure.
/// If it returns true, then n is probably prime.
bool IsProbablyPrime(Integer n) {
    if (n % 2 == 0) {
        return false;
    }

    // Checks if 2^(n-1) === 1 (mod n)
    return ModularExponentiation(2, n-1, n) == 1;
}

/// Finds a prime in the arithmetic progression kn + 1. By Dirichlet's theorem
/// such a prime exists for sure.
Integer FindPrimeInAP(const Integer n) {
    Integer maybe_prime = n + 1;

    bool flag_print = false;

    while (true) {
        if (maybe_prime >= INT_MAX && !flag_print) {
            printf("WARNING: FindPrimeInAP is returning a number that does not fit in int.\n");
            flag_print = true;
        } 
        else if (maybe_prime < 0) {
            fprintf(stderr, "ERROR: Overflow in FindPrimeInAP.\n");
            exit(1);
        }
        // Use short circuit evaluation for efficiency 
        if (IsProbablyPrime(maybe_prime) && IsPrime(maybe_prime)) {
            return maybe_prime;
        }

        maybe_prime += n;
    }
}

std::vector<nt::Integer> FindPrimesInAP(const Integer n, const size_t number) {

    std::vector<nt::Integer> out(number);
    size_t count = 0;

    Integer maybe_prime = n + 1;

    bool flag_print = false;

    while (true) {
        if (maybe_prime >= INT_MAX && !flag_print) {
            printf("WARNING: FindPrimeInAP is returning a number that does not fit in int.\n");
            flag_print = true;
        } 
        else if (maybe_prime < 0) {
            fprintf(stderr, "ERROR: Overflow in FindPrimeInAP.\n");
            exit(1);
        }
        // Use short circuit evaluation for efficiency 
        if (IsProbablyPrime(maybe_prime) && IsPrime(maybe_prime)) {
            out[count] = maybe_prime;
            count++;

            if (count == number) {
                return out;
            }
        }

        maybe_prime += n;
    }
}

/// Returns the prime divisors of n in sorted order with multiplicity
std::vector<Integer> PrimeDivisorsWithMultiplicity(Integer n) {
    std::vector<Integer> divisors; 
    Integer p = 2;
    while (n > 1) {
        while (n % p == 0) {
            divisors.push_back(p);
            n = n/p;
        }
        p++;
    }

    return divisors;
}

/// Returns the prime divisors of n in sorted order without multiplicity
std::vector<Integer> PrimeDivisors(Integer n) {
    std::vector<Integer> divisors = PrimeDivisorsWithMultiplicity(n);
    // Erase non-unique prime divisors
    divisors.erase(std::unique(divisors.begin(), divisors.end()), divisors.end());
    return divisors;
}

/// Returns a primitive root modulo the prime p
Integer PrimitiveRootModPrime(Integer p) {
    // prime divisors of phi(p) = p-1
    std::vector<Integer> prime_divisors = PrimeDivisors(p-1);
    for (Integer g=1; g < p; g++) {
        bool flag = false;
        for (auto q : prime_divisors) {
            if (ModularExponentiation(g, (p-1)/q, p) == 1) {
                flag = true;
                break;
            }
        }
        if (!flag) {
            return g;
        }
    }

    fprintf(stderr, "ERROR in PrimitiveRootModPrime: no primitive root found\n.");
    exit(1);
}

// Assumes that r and m are coprime.
Integer MultiplicativeInverse(const Integer value, const Integer m) {

    const Integer r = SafeMod(value, m);
    if (r == 0) {
        std::cerr << "ERROR: Multiplicative Inverse used with non-coprime arguments: " << value << " " << m << "\n";
        exit(1);
    }

    if (r == 1) {
        return 1;
    } 

    const Integer inverse = (-(m/r) * MultiplicativeInverse(m % r, m)) % m;

    return SafeMod(inverse, m);
}

Integer ChineseRemainderTheorem(const std::vector<Integer> &remainders, const std::vector<Integer> &moduli) {
    const size_t k = remainders.size();
    assert(moduli.size() == k);

    const auto mul = [](Integer acc, Integer m) -> Integer {
        return acc * m; 
    };

    const Integer M = std::accumulate(moduli.begin(), moduli.end(), (Integer) 1, mul);

    Integer solution = 0;
    for (size_t i = 0; i < k; i++) {
        const Integer r_i = SafeMod(remainders[i], moduli[i]);
        const Integer m_i = moduli[i];
        const Integer M_i = M / moduli[i];
        const Integer Inv_i = MultiplicativeInverse(M_i, m_i);

        const Integer n_i = (M_i * Inv_i) % M;
        solution = (solution + r_i * n_i) % M; 
    }

    return SafeMod(solution, M);
}


};