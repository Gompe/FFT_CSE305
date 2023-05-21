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
    nt::LL N = 1 << n;
    nt::LL p = nt::FindPrimeInAP(N);
    nt::LL g = nt::PrimitiveRootModPrime(p);

    printf("N = %d, p = %d, g = %d \n", N, p);
    
    std::vector<nt::LL> integers(N);
    std::vector<nt::LL> out(N);

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
