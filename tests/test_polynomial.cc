#include <polynomial/polynomial.h>

#include <tests/benchmark_timer.h>

#include <string>

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

void TestIntegerPolynomialMultiplication(const size_t degree = 10000, const int max_coef = 10000);
void ComparePolynomialMultiplication(const size_t degree);

int main() {
    std::string line(50, '-');
    std::cout << line << std::endl;

    std::cout << ">>>input size: 2^6\n";
    TestIntegerPolynomialMultiplication(1 << 6, 100);
    std::cout << line << std::endl;

    std::cout << ">>>input size: 2^12\n";
    TestIntegerPolynomialMultiplication(1 << 12, 100);
    std::cout << line << std::endl;

    std::cout << ">>>input size: 2^18\n";
    ComparePolynomialMultiplication(1 << 18);
    std::cout << line << std::endl;

    std::cout << ">>>input size: 2^24\n";
    ComparePolynomialMultiplication(1 << 24);
    std::cout << line << std::endl;
    
}

void TestIntegerPolynomialMultiplication(const size_t degree, const int max_coef) {
    std::cout << "Testing Integer Polynomial Multiplication with degree " << degree << "\n";

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
        
        // Similar purpose as std::vector::reserve
        PQ_round.SetCoefficient(PQ_real.Degree(), 1);
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

void ComparePolynomialMultiplication(const size_t degree) {
    std::cout << "Testing Integer Polynomial Multiplication with degree " << degree << "\n";

    auto generate_coefficients = [](auto num_coefs){
        std::vector<int> out(num_coefs);
        for (size_t i = 0; i < std::min<size_t>((1<<12), num_coefs); i++) {
            out[i] = rand() % 2;
        }

        out[num_coefs - 1] = 1;
        return out;
    };

    const size_t degree_P = degree;

    // Important to test polynomial with different degrees
    const size_t degree_Q = random() % degree;

    const Polynomial<int> P(generate_coefficients(degree_P));
    const Polynomial<int> Q(generate_coefficients(degree_Q));

    Polynomial<int> PQ_fft, PQ_round;

    auto make_fft = [&P, &Q, &PQ_fft](){
        PQ_fft = IntegerMultiply(P, Q);
    };

    auto make_rounding = [&P, &Q, &PQ_round](){
        auto PQ_real = RealMultiply(P, Q);
        
        // Similar purpose as std::vector::reserve
        PQ_round.SetCoefficient(PQ_real.Degree(), 1);
        for (size_t k = 0; k <= PQ_real.Degree(); k++) {
            PQ_round.SetCoefficient(k, RoundDouble(PQ_real[k]));
        }
    };

    timeFunction(make_fft, "FFT Polynomial Multiply");
    timeFunction(make_rounding, "FFT Real Polynomial Multiply");

    assert(PQ_fft.Degree() == PQ_round.Degree());
    for (size_t i = 0; i <= PQ_fft.Degree(); i++) {
        if (PQ_fft[i] != PQ_round[i]) {
            std::cout << "Coefficient " << i << ":\n";
            std::cout << "\tExpected: " << PQ_fft[i] << "\n";
            std::cout << "\tGot: " << PQ_round[i] << "\n";
        }
    }
}