#include <polynomial/polynomial.h>
#include <iostream>

#include <stdio.h>
#include <string.h>

int main(int argc, char **argv) {
    if (argc != 2 || strlen(argv[1]) != 1) {
        printf("Usage: %s [MODE]\n", argv[0]);
        printf("\t [MODE] can take the following values:\n");
        printf("\t\t `i` for integer polynomial multiplication\n");
        printf("\t\t `r` for real polynomial multiplication\n");
        exit(1);
    }

    char mode = argv[1][0];

    auto write_polynomial = [](auto &P) {
        for (size_t i = 0; i <= P.Degree(); i++) {
            if (i == 0) {
                std::cout << P[0];
            }
            else if (i == 1) {
                std::cout << P[1] << "X";
            }
            else {
                std::cout << P[i] << "X^" << i;
            }

            if (i != P.Degree()) {
                std::cout << " + ";
            }
        }
        std::cout << "\n";
    };

    auto read_polynomial = [](size_t deg, auto& P){
        for (size_t i = 0; i <= deg; i++) {
            auto coef = P[0];
            std::cin >> coef;
            P.SetCoefficient(i, coef); 
        }
    };

    if (mode == 'i') {
        size_t deg_A, deg_B;
        Polynomial<int> A, B;


        std::cout << "Type the degree of the first polynomial: ";
        std::cin >> deg_A;

        printf("\nEnter the %lu coefficients of the first polynomial:\n", deg_A + 1);
        read_polynomial(deg_A, A);

        std::cout << "Your first polynomial is:\n";
        write_polynomial(A);

        std::cout << "Type the degree of the second polynomial: ";
        std::cin >> deg_B;

        printf("\nEnter the %lu coefficients of the second polynomial:\n", deg_B + 1);
        read_polynomial(deg_B, B);

        std::cout << "Your second polynomial is:\n";
        write_polynomial(B);

        auto AB = IntegerMultiply(A, B);

        std::cout << "The output polynomial is:\n";
        write_polynomial(AB);
    }

    if (mode == 'r') {
        size_t deg_A, deg_B;
        Polynomial<FloatType> A, B;


        std::cout << "Type the degree of the first polynomial: ";
        std::cin >> deg_A;

        printf("\nEnter the %lu coefficients of the first polynomial:\n", deg_A + 1);
        read_polynomial(deg_A, A);

        std::cout << "Your first polynomial is:\n";
        write_polynomial(A);

        std::cout << "Type the degree of the second polynomial: ";
        std::cin >> deg_B;

        printf("\nEnter the %lu coefficients of the second polynomial:\n", deg_B + 1);
        read_polynomial(deg_B, B);

        std::cout << "Your second polynomial is:\n";
        write_polynomial(B);

        auto AB = RealMultiply(A, B);

        std::cout << "The output polynomial is:\n";
        write_polynomial(AB);
    }

}