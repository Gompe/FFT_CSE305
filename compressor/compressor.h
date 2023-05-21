#pragma once

#ifndef COMPRESSOR_H_
#define COMPRESSOR_H_

#define DEFAULT_NUM_FREQUENCIES 2 

#include <vector>

#include <core/dft.h>

namespace compressor {

/// Stores one component of the Discrete Fourier Transform of the Data.
struct EncodedItem {
    int index;
    Complex value;

    EncodedItem() = default;
    EncodedItem(int index, Complex value);
};

using Data = std::vector<FloatType>;
using EncodedData = std::vector<EncodedItem>;

EncodedData Compress(const std::vector<FloatType> &data, int num_frequencies=DEFAULT_NUM_FREQUENCIES);
std::vector<FloatType> Decompress(const EncodedData &encoded_data, const int N);


/// utils
std::vector<FloatType> ReadDataFromStdin();
// std::vector<FloatType> ReadDataFromFile(std::string filename);

void WriteDataToStdout(const std::vector<FloatType> &data);

}; // namespace compressor

#endif