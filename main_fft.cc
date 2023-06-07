#include <core/parallel_dft.h>
#include <vector>
#include <string.h>

int main(int argc, char **argv) {

    char mode = 'd';

    if (argc != 1) {
        if (argc != 2  || strlen(argv[1]) != 1) {
            printf("Usage: %s [MODE]\n", argv[0]);
            printf("\t [MODE] (optional) can take the following values:\n");
            printf("\t\t `d` (default) for the DFT transformation\n");
            printf("\t\t `i` for the inverse dft transformation\n");
            exit(1);
        }

        mode = argv[1][0];
        if (mode != 'd' && mode != 'i') {
            printf("Invalid mode.\n");
            exit(1);
        }
    }

    size_t N;

    printf("Write the size of your sequence (must be a power of 2): ");
    std::cin >> N;

    if (!fft_utils::IsPowerOfTwo(N)) {
        printf("ERROR: Your sequence size is not a power of 2.\n");
        exit(1);
    }

    printf("Write your sequence (Complex Numbers Allowed)\n");

    std::vector<Complex> data(N);
    for (size_t i = 0; i < N; i++) {
        std::cin >> data[i];
    }

    std::cout << "Your data is:\n"; 
    for (size_t i = 0; i < N; i++) {
        std::cout << data[i] << "\n";
    }

    std::cout << "\n";

    if (mode == 'd') {
        iterative_fft::ParallelDFT(data.begin(), data.end(), data.begin(), FixedThreadsParallelizer{});
    }

    if (mode == 'i') {
        iterative_fft::ParallelIDFT(data.begin(), data.end(), data.begin(), FixedThreadsParallelizer{});
    }

    std::string transform;
    if (mode == 'i') {
        transform = "IDFT";
    } else {
        transform = "DFT";
    }

    std::cout << "The " << transform << " of your data is:\n";
    for (size_t i = 0; i < N; i++) {
        std::cout << data[i] << "\n";
    }

}