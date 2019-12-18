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
    
    for(int i = 1; i <= 16; i*=2){
        omp_set_num_threads(i);
        unit_test::classicalLSH_find_1D_approx_rnear();
        unit_test::coveringLSH_find_1D_approx_rnear();
        std::cout << std::endl;
    }
    
    // complete the code
    return 0;
}
