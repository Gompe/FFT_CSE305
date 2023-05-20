#include <iostream>
#include <complex>
#include <cmath>
#include <vector>

#include <random>

#include "fft_types.h"
#include "my_timer.h"

#include "dft.h"

template <typename T>
void PrintArray(T* array, size_t N) {
    for (size_t i=0; i<N; i++)
        std::cout << array[i] << " ";
    std::cout << std::endl;
}

void truncate(Complex *array, size_t N, FloatType limit) {
    for(size_t i=0; i<N; i++) {
        FloatType real, imag;
        if (std::abs(array[i].real()) < limit) {
            real = 0;
        } else {
            real = array[i].real();
        }
        if (std::abs(array[i].imag()) < limit) {
            imag = 0;
        } else {
            imag = array[i].imag();
        }
        array[i] = Complex(real, imag);
    }
}

bool checkIsClose(Complex *a, Complex *b, size_t N, FloatType eps=1E-3) {
    for (size_t i=0; i<N; i++) {
        if(norm(a[i] - b[i]) > eps)
            return false;
    }
    return true;
}

// const int N = 1 << 20;
const int N = 1 << 20;

Complex x[N], y[N], z[N];

int main()
{   
    for (int i=0; i<N; i++) {
        x[i] = rand() % 2;
    }

    auto func1 = [](){ RecursiveFft::Transform(x, y, N); };
    auto func2 = [](){ IterativeFft::Transform(x, z, N);};

    timeFunction(func1, "Base DFT");
    timeFunction(func2, "Iterative");

    if(checkIsClose(y, z, N)) {
        std::cout << "Code ok" << std::endl;
    } else {
        std::cout << "Code not ok" << std::endl;
    }

    return 0;
}