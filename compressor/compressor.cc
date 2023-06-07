#include <compressor/compressor.h>
#include <numeric>

namespace compressor {

EncodedItem::EncodedItem(int index, Complex value) 
: index(index), value(value) {}

// Compress / Decompress functions
EncodedData Compress(const Data &data, int num_frequencies) {

    const int N = (data.size() == 1) ? 1 : fft_utils::PowerOfTwo(fft_utils::IntLog2(data.size() - 1) + 1);

    std::vector<Complex> resized_data(N);
    std::copy(data.begin(), data.end(), resized_data.begin());

    const double data_average = std::accumulate(data.begin(), data.end(), 0.) / data.size();
    std::fill(resized_data.begin() + data.size(), resized_data.end(), data_average);


    // At most N frequencies can be captured by the FFT
    num_frequencies = std::min(N, num_frequencies);

    // Data on the frequency domain
    std::vector<Complex> data_freq(N);
    iterative_fft::DFT(resized_data.begin(), resized_data.end(), data_freq.begin());

    // Populate the EncodedData struct with frequency data
    EncodedData compressed_data(N);
    for (int index=0; index < N; index++) {
        compressed_data[index] = EncodedItem(index, data_freq[index]);
    }

    // Partial sort - The first num_frequencies EncodedItem objects correspond
    // to the num_frequencies items of data_freq with largest complex norm.
    std::nth_element(compressed_data.begin(), compressed_data.begin() + num_frequencies,
                    compressed_data.end(), [](const EncodedItem &a, const EncodedItem &b) {
                        return (std::norm(a.value) > std::norm(b.value));
                    });
    
    // Discard all frequencies with small complex norm
    compressed_data.erase(compressed_data.begin() + num_frequencies, compressed_data.end());

    return compressed_data;
}

Data Decompress(const EncodedData &encoded_data, const int output_size) {

    const int N = (output_size == 1) ? 1 : fft_utils::PowerOfTwo(fft_utils::IntLog2(output_size - 1) + 1);

    std::vector<Complex> frequency_data(N);
    std::fill(frequency_data.begin(), frequency_data.end(), (Complex) 0);

    for (const auto &item : encoded_data) {
        assert(item.index >= 0 && item.index < N);
        frequency_data[item.index] = item.value;
    }

    std::vector<Complex> predecoded_data(N);
    iterative_fft::IDFT(frequency_data.begin(), frequency_data.end(), predecoded_data.begin());

    std::vector<FloatType> decoded_data(N);
    std::transform(predecoded_data.begin(), predecoded_data.end(), decoded_data.begin(),
                    [](const Complex &x){ return x.real(); });

    decoded_data.erase(decoded_data.begin() + output_size, decoded_data.end());
    return decoded_data;
}


/// utils
Data ReadDataFromStdin() {
    int size_sequence;
    std::cin >> size_sequence;

    if (size_sequence <= 0) {
        printf("ERROR: the size of the sequence must be positive.\n");
        exit(1);
    }

    Data data(size_sequence);

    for (int i=0; i < size_sequence; i++) {
        std::cin >> data[i];
    }

    return data;
}

void WriteDataToStdout(const Data &data) {
    std::cout << data.size() << "\n";
    for (const auto &elem : data) {
        std::cout << elem << " ";
    }
    std::cout << "\n";
}

}; // namespace compressor