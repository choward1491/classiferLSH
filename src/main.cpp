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
    
    unit_test::simple_2D_binary_classification_classical(8, 5, 512);
    unit_test::simple_2D_binary_classification_covering(8, 5, 512);
    // complete the code
    return 0;
}
