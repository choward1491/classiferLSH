//
//  main.cpp
//  classifierLSH
//
//  Created by Christian Howard on 12/13/19.
//  Copyright © 2019 Christian Howard. All rights reserved.
//

#include <iostream>
#include "omp.h"

int main(int argc, const char * argv[]) {
    
    // set the number of threads to use
    omp_set_num_threads(8);
    
    // run some openmp stuff
    #pragma omp parallel
    #pragma omp critical
        std::cout << "Greetings from thread "<<omp_get_thread_num()<<std::endl;
        return 0;
    
    // complete the code
    return 0;
}
