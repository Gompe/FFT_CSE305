#include <iostream>
#include <random>

#include <compressor/compressor.h>

static void PrintData(const compressor::Data &data) {
    std::cout << "Data " << data.size() << " elements\n";
    for(const auto &elem : data) {
        std::cout << elem << " ";
    }
    std::cout << "\n";
}

static void PrintEncodedData(const compressor::EncodedData &data) {
    std::cout << "Encoded Data " << data.size() << " elements\n";
    for(const auto &elem : data) {
        std::cout << "(" << elem.index  << ", " << elem.value << ") ";
    }
    std::cout << "\n";
}

static void PrintFft(const std::vector<Complex> &sequence) {
    std::cout << "Fft Coefficients " << sequence.size() << " elements\n";
    for(const auto &elem : sequence) {
        std::cout << elem << " ";
    }
    std::cout << "\n";
}

int main(int argc, char **argv) {

    if (argc != 2) {
        fprintf(stderr, "usage: %s num_frequencies", argv[0]);
        exit(1);
    }

    int num_frequencies = std::stoi(argv[1]);

    compressor::Data data = compressor::ReadDataFromStdin();
    compressor::EncodedData compressed_data = compressor::Compress(data, num_frequencies);

    compressor::Data reconstructed_data = compressor::Decompress(compressed_data, data.size());
    compressor::WriteDataToStdout(reconstructed_data);
}