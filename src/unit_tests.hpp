//
//  unit_tests.hpp
//  classifierLSH
//
//  Created by Christian Howard on 12/14/19.
//  Copyright © 2019 Christian Howard. All rights reserved.
//

#ifndef unit_tests_hpp
#define unit_tests_hpp

namespace unit_test {

// tests for the bitvector class
bool bitvec_and();
bool bitvec_xor();
bool bitvec_dotproduct();
bool bitvec_norm(); // Hamming metric
bool bitvec_equality();
bool bitvec_dist(); // Hamming metric

// tests for the classicalLSH r-near DS
bool classicalLSH_find_1D_approx_rnear();

// tests for the coveringLSH r-near DS
bool coveringLSH_find_1D_approx_rnear();

// run the tests
void run_tests();

}

#endif /* unit_tests_hpp */