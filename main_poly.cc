#include <polynomial/polynomial.h>

template < class T >
static void PrintPolynomial(const Polynomial<T> &P) {
    std::cout << "Polynomial\n";
    for (size_t i = 0; i <= P.Degree(); i++) {
        std::cout << P[i] << " ";
    }
    std::cout << "\n";
}

int main() {
    PolynomialCoefficients<int> coeff = {1, 2, 3};
    Polynomial<int> P(coeff);

    PrintPolynomial(P);

    PrintPolynomial(NaiveMultiply(P, P));
    PrintPolynomial(RealMultiply(P, P));
    PrintPolynomial(IntegerMultiply(P, P));
}