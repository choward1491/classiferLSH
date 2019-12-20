//
//  unit_tests.hpp
//  classifierLSH
//
//  Created by Christian Howard on 12/14/19.
//  Copyright © 2019 Christian Howard. All rights reserved.
//

#ifndef unit_tests_hpp
#define unit_tests_hpp

#include <cstdint>

namespace unit_test {

// tests for the bitvector class
bool bitvec_and();
bool bitvec_xor();
bool bitvec_dotproduct();
bool bitvec_norm(); // Hamming metric
bool bitvec_equality();
bool bitvec_dist(); // Hamming metric

// tests for the classicalLSH r-near DS
bool classicalLSH_find_1D_approx_rnear(int num_threads, uint32_t U = 8*1024);
bool classicalLSH_find_1D_approx_knearest(int num_threads, int k, uint32_t U = 8*1024);

// tests for the coveringLSH r-near DS
bool coveringLSH_find_1D_approx_rnear(int num_threads, uint32_t U = 8*1024);
bool coveringLSH_find_1D_approx_knearest(int num_threads, int k, uint32_t U = 8*1024);

// tests for classification
void simple_2D_binary_classification_classical(int num_threads, int k, uint32_t U = 256);
void simple_2D_binary_classification_covering(int num_threads, int k, uint32_t U = 256);

// run the tests
void run_tests();

}

#endif /* unit_tests_hpp */
