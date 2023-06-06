#include <polynomial/polynomial.h>

#include <tests/benchmark_timer.h>


template < class T >
static void PrintPolynomial(const Polynomial<T> &P) {
    std::cout << "Polynomial\n";
    for (size_t i = 0; i <= P.Degree(); i++) {
        std::cout << P[i] << " ";
    }
    std::cout << "\n";
}

static int RoundDouble(const double x) {
    return static_cast <int>(x + 0.5 - (x < 0.0));
}

void TestChineseRemainderTheorem();
void TestModularInverse();
void TestIntegerPolynomialMultiplication(const size_t degree = 10000, const int max_coef = 10000);

int main() {
    TestChineseRemainderTheorem();
    TestModularInverse();
    TestIntegerPolynomialMultiplication(40000, 1000);
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

void TestIntegerPolynomialMultiplication(const size_t degree, const int max_coef) {
    std::cout << "Testing Integer Polynomial Multiplication with degree " << 10000 << "\n";

    auto generate_coefficients = [max_coef](auto num_coefs){
        std::vector<int> out(num_coefs);
        for (size_t i = 0; i < num_coefs; i++) {
            out[i] = (random() % (2*max_coef)) - max_coef;
        }
        return out;
    };

    const size_t degree_P = degree;

    // Important to test polynomial with different degrees
    const size_t degree_Q = random() % degree;

    const Polynomial<int> P(generate_coefficients(degree_P));
    const Polynomial<int> Q(generate_coefficients(degree_Q));

    Polynomial<int> PQ_naive, PQ_fft, PQ_round;

    auto make_naive = [&P, &Q, &PQ_naive](){
        PQ_naive = NaiveMultiply(P, Q);
    };

    auto make_fft = [&P, &Q, &PQ_fft](){
        PQ_fft = IntegerMultiply(P, Q);
    };

    auto make_rounding = [&P, &Q, &PQ_round](){
        auto PQ_real = RealMultiply(P, Q);
        std::vector<int> coefs(PQ_real.Degree() + 1);
        for (size_t k = 0; k <= PQ_real.Degree(); k++) {
            PQ_round.SetCoefficient(k, RoundDouble(PQ_real[k]));
        }
    };

    timeFunction(make_naive, "Naive Polynomial Multiply");
    timeFunction(make_fft, "FFT Polynomial Multiply");
    timeFunction(make_rounding, "FFT Real Polynomial Multiply");

    assert(PQ_naive.Degree() == PQ_fft.Degree());
    for (size_t i = 0; i <= PQ_naive.Degree(); i++) {
        if (PQ_naive[i] != PQ_fft[i]) {
            std::cout << "Coefficient " << i << ":\n";
            std::cout << "\tExpected: " << PQ_naive[i] << "\n";
            std::cout << "\tGot: " << PQ_fft[i] << "\n";
            std::cout << "\tReal Multiply: " << RealMultiply(P, Q)[i] << "\n";
        }
        assert(PQ_naive[i] == PQ_fft[i]);
    }

    assert(PQ_naive.Degree() == PQ_round.Degree());
    for (size_t i = 0; i <= PQ_naive.Degree(); i++) {
        if (PQ_naive[i] != PQ_round[i]) {
            std::cout << "Coefficient " << i << ":\n";
            std::cout << "\tExpected: " << PQ_naive[i] << "\n";
            std::cout << "\tGot: " << PQ_round[i] << "\n";
        }
    }
}