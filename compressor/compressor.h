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
Data Decompress(const EncodedData &encoded_data, const int N);


/// utils
Data ReadDataFromStdin();
// TODO: std::vector<FloatType> ReadDataFromFile(std::string filename);

void WriteDataToStdout(const Data &data);

}; // namespace compressor

#endif