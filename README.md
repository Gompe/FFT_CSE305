# FFT_CSE305
## Pedro Cabral

### How to run

#### Usage

There are three main (non-test) files: `main_fft.cc`, `main_polynomial.cc`, and `/compressor/main.cc`. In order to compile any of these files, simple type the appropriate make command:
>`$ make fft`
>`$ make polynomial`
>`$ make compressor`

The fft program can take any sequence of size power of 2 from stdin and an option `d` (for DFT) or `i` (for IDFT) and compute the corresponding Parallel Fast Fourier Transform. 

The polynomial program can take any sequence of integer/real polynomials from stdin and compute their product with the corresponding fast polynomial multiplication algorithm.

The compressor program can take any sequence from stdin, compress it internally, and output the data reconstruction to stdout. It is easier to use with in combination to the provided python file `compressor_demo.py`. This file will use matplotlib to plot the original data and the reconstructed data. Moreover, it will provide some information about the reconstruction error.

#### Timing Experiments

For the timing experiments, run the appropriate make commands. In order to verify the timing of the FFT experiments run:

>$ make test_dft

in order to to see the timing of the polynomial experiments, run:

>$ make test_polynomial

