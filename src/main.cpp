//
//  main.cpp
//  classifierLSH
//
//  Created by Christian Howard on 12/13/19.
//  Copyright © 2019 Christian Howard. All rights reserved.
//

#include <iostream>
#include "omp.h"
#include "bitvector.hpp"
#include "classicalLSH.hpp"
#include "unit_tests.hpp"

int main(int argc, const char * argv[]) {
    
    // set the number of threads to use
    omp_set_num_threads(8);
    
    unit_test::run_tests();
    unit_test::classicalLSH_find_1D_approx_rnear();
    
    // complete the code
    return 0;
}
