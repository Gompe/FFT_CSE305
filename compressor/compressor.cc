#include <compressor/compressor.h>

namespace compressor {

EncodedItem::EncodedItem(int index, Complex value) 
: index(index), value(value) {}

static std::vector<Complex> FftTransform(const std::vector<FloatType> &data) {
    const int N = data.size();
    // Initialize location of output
    std::vector<Complex> data_dft(N);

    // Cast data to Complex
    std::vector<Complex> data_astype_complex(N);
    for (int i=0; i < N; i++) {
        data_astype_complex[i] = (Complex) data[i];
    }
    
    // Compute FFT
    IterativeFft::Transform(&data_astype_complex[0], &data_dft[0], N);

    return data_dft;
}

// Compress / Decompress functions
EncodedData Compress(const std::vector<FloatType> &data, int num_frequencies) {
    const int N = data.size();

    // At most N frequencies can be captured by the FFT
    num_frequencies = std::min(N, num_frequencies);

    // Data on the frequency domain
    std::vector<Complex> data_freq = FftTransform(data);

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

std::vector<FloatType> Decompress(const EncodedData &encoded_data, const int N) {
    std::vector<Complex> frequency_data(N);
    for (const auto &item : encoded_data) {
        frequency_data[item.index] = item.value;
    }

    std::vector<Complex> predecoded_data(N);
    IterativeFft::InverseTransform(&frequency_data[0], &predecoded_data[0], N);

    std::vector<FloatType> decoded_data(N);
    std::transform(predecoded_data.begin(), predecoded_data.end(), decoded_data.begin(),
                    [](const Complex &x){ return x.real(); });

    return decoded_data;
}


/// utils
std::vector<FloatType> ReadDataFromStdin() {
    int size_sequence;
    std::cin >> size_sequence;

    std::vector<FloatType> data(size_sequence);

    for (int i=0; i < size_sequence; i++) {
        std::cin >> data[i];
    }

    return data;
}

void WriteDataToStdout(const std::vector<FloatType> &data) {
    std::cout << data.size() << "\n";
    for (const auto &elem : data) {
        std::cout << elem << " ";
    }
    std::cout << "\n";
}

}; // namespace compressor