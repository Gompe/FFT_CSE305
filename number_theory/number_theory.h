#pragma once
#pragma GCC optimize "trapv"


#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace nt {

using Integer = long long;

Integer SafeMod(const Integer a, const Integer m);

/// Returns a^b (mod n) in time complexity O(log b)
/// The returned value is in the range [0...n-1]
Integer ModularExponentiation(Integer a, Integer b, Integer n);

/// Returns if n is Prime - Always returns the correct result
bool IsPrime(Integer n);

/// Returns if n is probably prime. 
/// If it returns false, then n is composite for sure.
/// If it returns true, then n is probably prime.
bool IsProbablyPrime(Integer n) ;

/// Finds a prime in the arithmetic progression kn + 1. By Dirichlet's theorem
/// such a prime exists for sure. 
Integer FindPrimeInAP(const Integer n);
std::vector<nt::Integer> FindPrimesInAP(const Integer n, const size_t number);

/// Returns the prime divisors of n in sorted order with multiplicity
std::vector<Integer> PrimeDivisorsWithMultiplicity(Integer n);

/// Returns the prime divisors of n in sorted order without multiplicity
std::vector<Integer> PrimeDivisors(Integer n);

/// Returns a primitive root modulo the prime p
Integer PrimitiveRootModPrime(Integer p);

Integer MultiplicativeInverse(const Integer value, const Integer m);

/// The Chinese Remainder Theorem
/// Given a vector of remainders (r1,r2,...,rk) and pairwise coprime moduli (m1,m2,...,mk)
/// returns the ONLY remainder r (mod. m1m2...mk) such that ri === r (mod. mi) for all i.
Integer ChineseRemainderTheorem(const std::vector<Integer> &remainders, const std::vector<Integer> &moduli);

};