#pragma once

#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <algorithm>

namespace nt {

using LL = long long;

/// Returns a^b (mod n) in time complexity O(log b)
/// The returned value is in the range [0...n-1]
LL ModularExponentiation(LL a, LL b, LL n);

/// Returns if n is Prime - Always returns the correct result
bool IsPrime(LL n);

/// Returns if n is probably prime. 
/// If it returns false, then n is composite for sure.
/// If it returns true, then n is probably prime.
bool IsProbablyPrime(LL n) ;

/// Finds a prime in the arithmetic progression kn + 1. By Dirichlet's theorem
/// such a prime exists for sure.
LL FindPrimeInAP(LL n);

/// Returns the prime divisors of n in sorted order with multiplicity
std::vector<LL> PrimeDivisorsWithMultiplicity(LL n);

/// Returns the prime divisors of n in sorted order without multiplicity
std::vector<LL> PrimeDivisors(LL n);

/// Returns a primitive root modulo the prime p
LL PrimitiveRootModPrime(LL p);

};