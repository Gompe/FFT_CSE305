#include <iostream>

#include <number_theory/number_theory.h>

int main() {
    // int k = 5;
    // std::cout << "k: " << k << std::endl;
    // std::cout << "Is k probably prime " << nt::IsProbablyPrime(k) << std::endl;
    // std::cout << "Is k prime " << nt::IsPrime(k) << std::endl;

    int n;
    std::cin >> n;

    nt::LL N = 1 << n;
    std::cout << "My number N" << std::endl;
    std::cout << N << std::endl;

    nt::LL p = nt::FindPrimeInAP(N);
    std::cout << "Prime p" << std::endl;
    std::cout << p << std::endl;

    std::vector<nt::LL> divisors = nt::PrimeDivisors(p-1);
    std::cout << "Prime divisors of p-1\n";
    for (auto q : divisors) {
        std::cout << q << " ";
    }
    std::cout << std::endl;

    nt::LL g = nt::PrimitiveRootModPrime(p);
    std::cout << "Prime root g" << std::endl;
    std::cout << g << std::endl;
}
