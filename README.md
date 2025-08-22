Private top-k on CKKS
=====================================

## Installation

This code is based on OpenFHE
* [OpenFHE documentation](https://openfhe-development.readthedocs.io/en/latest/)
* [Design paper for OpenFHE](https://eprint.iacr.org/2022/915)
* [OpenFHE website](https://openfhe.org)

Note our implementation is on version 1.1.4

After installing OpenFHE, operate following on "privtopk" directory:

* cd build
* cmake ..
* make
* ./test_topk 
or else, with proper input arguments.

# Tests for Topk

- The order of input is as the follows: ./test_topk size numk iteration 

* example:
* ./test_topk 65536 16 8
Run test on input array size 65536, top-16, 8 iterations.

If you omit an argument, it runs with every prepared values. 
* ./test_topk 65536
Run test on input array size 65536, top-1, 2, 4, 8, 16, 8 iterations at default.

# Tests for Sort

- The order of input is as the follows: ./test_sort size iteration 

* example:
* ./test_sort 256 8
Run test on input array size 256, 8 iterations.

If you omit an argument, it runs with every prepared values. 
* ./test_sort 
Run test on input array size 8 to 256, 8 iterations at default.

# Tests for varying threads

- The order of input is as the follows: ./test_abl number_of_initial_ciphertexts number_of_threads iteration 

* example:
* ./test_abl 8 64 8
Run test with: the number of input ciphertexts 8, number of threads 64, 8 iterations.

If you omit an argument, it runs with every prepared values. 
* ./test_sort 8
Run test with: the number of input ciphertexts 8, number of threads 1~128, 8 iterations at default.



## Note on Structure of the Library

# Files

- algorithms: collection of our algorithms specified in the papers.

- test_topk / test_abl / test_sort: main function for varous test.

- testcode: code containg test pipeline of our algorithms.

- utils: various algorithms for test, including reading/writing of data, random number generation, and precision estimation.
