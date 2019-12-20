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
    
    /*
    for(int i = 1; i <= 16; i*=2){
        unit_test::classicalLSH_find_1D_approx_rnear(i);
        unit_test::coveringLSH_find_1D_approx_rnear(i);
        std::cout << std::endl;
    }*/
    
    unit_test::classicalLSH_find_1D_approx_knearest(8, 5);
    
    unit_test::coveringLSH_find_1D_approx_knearest(8, 5);
    
    // complete the code
    return 0;
}
