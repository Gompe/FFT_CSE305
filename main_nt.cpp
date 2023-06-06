#include <iostream>

#include <number_theory/number_theory.h>
#include <core/modular_fft.h>

template <typename T>
void PrintVec(std::vector<T> vec) {
    for (auto x : vec) std::cout << x << " ";
    std::cout << "\n";
}

int main() {

    int n = 5;
    nt::Integer N = 1 << n;
    nt::Integer p = nt::FindPrimeInAP(N);
    nt::Integer g = nt::PrimitiveRootModPrime(p);

    printf("N = %d, p = %d, g = %d \n", N, p);
    
    std::vector<nt::Integer> integers(N);
    std::vector<nt::Integer> out(N);

    for (int i=0; i<N; i++) {
        integers[i] = i;
    }

    ModularFftTransform(integers.begin(), integers.end(), out.begin(), p, g);

    printf("Integers:\n");
    PrintVec(integers);

    printf("\nDFT:\n");
    PrintVec(out);

    ModularFftInverseTransform(out.begin(), out.end(), out.begin(), p, g);

    printf("\nIDFT:\n");
    PrintVec(out);
}
