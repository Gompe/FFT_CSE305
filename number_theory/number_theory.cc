#include <number_theory/number_theory.h>

#include <limits.h>

namespace nt {
    
/// Returns a^b (mod n) in time complexity O(log b)
/// The returned value is in the range [0...n-1]
LL ModularExponentiation(LL a, LL b, LL n) {
    LL remainder = 1;
    while(b != 0) {
        if (b % 2 == 1) {
            remainder = (remainder * a) % n;
        }
        a = (a * a) % n;
        b = b/2;
    }

    return remainder;
}

/// Returns if n is Prime - Always returns the correct result
bool IsPrime(LL n) {
    for (LL d=2; d*d <= n; d++) {
        if (n % d == 0) {
            return false;
        }
    }
    return true;
}

/// Returns if n is probably prime. 
/// If it returns false, then n is composite for sure.
/// If it returns true, then n is probably prime.
bool IsProbablyPrime(LL n) {
    if (n % 2 == 0) {
        return false;
    }

    // Checks if 2^(n-1) === 1 (mod n)
    return ModularExponentiation(2, n-1, n) == 1;
}

/// Finds a prime in the arithmetic progression kn + 1. By Dirichlet's theorem
/// such a prime exists for sure.
LL FindPrimeInAP(LL n) {
    LL maybe_prime = n + 1;
    bool flag_print = false;

    while (true) {
        if (maybe_prime >= INT_MAX && !flag_print) {
            printf("WARNING: FindPrimeInAP is returning a number that does not fit in int.\n");
            flag_print = true;
        } 
        else if (maybe_prime < 0) {
            fprintf(stderr, "WARNING: Overflow in FindPrimeInAP.\n");
            exit(1);
        }
        // Use short circuit evaluation for efficiency 
        if (IsProbablyPrime(maybe_prime) && IsPrime(maybe_prime)) {
            return maybe_prime;
        }

        maybe_prime += n;
    }
}

/// Returns the prime divisors of n in sorted order with multiplicity
std::vector<LL> PrimeDivisorsWithMultiplicity(LL n) {
    std::vector<LL> divisors; 
    LL p = 2;
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
std::vector<LL> PrimeDivisors(LL n) {
    std::vector<LL> divisors = PrimeDivisorsWithMultiplicity(n);
    // Erase non-unique prime divisors
    divisors.erase(std::unique(divisors.begin(), divisors.end()), divisors.end());
    return divisors;
}

/// Returns a primitive root modulo the prime p
LL PrimitiveRootModPrime(LL p) {
    // prime divisors of phi(p) = p-1
    std::vector<LL> prime_divisors = PrimeDivisors(p-1);
    for (LL g=1; g < p; g++) {
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

};