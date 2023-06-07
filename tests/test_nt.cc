#include <iostream>

#include <number_theory/number_theory.h>
#include <core/modular_fft.h>

template <typename T>
void PrintVec(std::vector<T> vec) {
    for (auto x : vec) std::cout << x << " ";
    std::cout << "\n";
}

void TestChineseRemainderTheorem();
void TestModularInverse();
void TestModularFFT();

int main() {
    TestChineseRemainderTheorem();
    TestModularInverse();
    TestModularFFT();
}

void TestModularFFT() {
    int n = 5;
    nt::Integer N = 1 << n;
    nt::Integer p = nt::FindPrimeInAP(N);
    nt::Integer g = nt::PrimitiveRootModPrime(p);

    printf("N = %lld, p = %lld, g = %lld \n", N, p, g);
    
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

void TestChineseRemainderTheorem() {
    auto test_instance = [](auto remainders, auto moduli) {
        auto out = nt::ChineseRemainderTheorem(remainders, moduli);
        assert(remainders.size() == moduli.size());
        for (size_t i = 0; i < remainders.size(); i++) {
            if(nt::SafeMod(remainders[i], moduli[i]) != nt::SafeMod(out, moduli[i])) {
                std::cerr << "FAIL: TestChineseRemainderTheorem\n";

                std::cerr << "remainders:\n";
                for(auto r : remainders) {
                    std::cerr << r << " ";
                } 
                std::cerr << "\n";
                std::cerr << "moduli:\n";
                for(auto m : moduli) {
                    std::cerr << m << " ";
                } 
                std::cerr << "\n";
                std::cerr << "output: " << out << "\n";
                for(auto m : moduli) {
                    std::cerr << nt::SafeMod(out, m) << " ";
                } 
                std::cerr << "\n\n";
                return;
            }
        }
    };

    // Test 1
    {
        std::vector<nt::Integer> remainders = {1, 2, 4, 3, 8 , -3   , -3};
        std::vector<nt::Integer> moduli     = {2, 3, 5, 7, 11, 65537, 163841};
        test_instance(remainders, moduli);
    }

    // Test 2
    {
        std::vector<nt::Integer> remainders = {65534, 163838};
        std::vector<nt::Integer> moduli     = {65537, 163841};
        test_instance(remainders, moduli);
    }

}

void TestModularInverse() {
    const std::vector<nt::Integer> primes = {2, 3, 5, 7, 11, 13, 17, 19, 2017, 65537, 163841, 557057};
    const std::vector<nt::Integer> remainders = {1, 2, 4, 7, 8, 25903, 1294,
     19251, 557054, 65537 * 163841LL, 65537 * 557057LL, 163841 * 557057LL, -3, -5, -9};

    for (auto p : primes) {
        for (auto r : remainders) {
            if (nt::SafeMod(r, p) == 0) {
                continue;
            }
            else {
                if(nt::MultiplicativeInverse(r, p) != nt::ModularExponentiation(r, p-2, p)) {
                    std::cout << "FAIL: TestModularInverse with inputs " << r << " " << p << "\n";
                    std::cout << "\tExpected: " << nt::ModularExponentiation(r, p-2, p) << "\n";
                    std::cout << "\tGot: " << nt::MultiplicativeInverse(r, p) << "\n";
                }
            }
        }
    }
}