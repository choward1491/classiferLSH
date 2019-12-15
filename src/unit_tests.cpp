//
//  unit_tests.cpp
//  classifierLSH
//
//  Created by Christian Howard on 12/14/19.
//  Copyright © 2019 Christian Howard. All rights reserved.
//

#include "unit_tests.hpp"
#include "bitvector.hpp"
#include "classicalLSH.hpp"
#include "coveringLSH.hpp"
#include <cmath>
#include <iostream>
#include <map>
#include <string>

namespace unit_test {

// tests for the bitvector class
bool bitvec_and() {
    util::bitvec bv, bv2;
    bv.set_parameters(3, 27);
    bv2.set_parameters(3, 27);
    
    bv.set(0,3);
    bv.set(1,4);
    bv.set(2,27);
    
    bv2.set(1,3);
    
    return (bv2 & bv).norm() == 3;
}
bool bitvec_xor() {
    util::bitvec bv, bv2;
    bv.set_parameters(3, 27);
    bv2.set_parameters(3, 27);
    
    bv.set(0,3);
    bv.set(1,4);
    bv.set(2,27);
    
    bv2.set(1,3);
    return (bv2 ^ bv).norm() == (3 + 1 + 27);
}
bool bitvec_dotproduct() {
    util::bitvec bv, bv2;
    bv.set_parameters(3, 27);
    bv2.set_parameters(3, 27);
    
    bv.set(0,3);
    bv.set(1,4);
    bv.set(2,27);
    
    bv2.set(1,3);
    
    return bv.dot(bv2) == 3;
}
bool bitvec_norm() { // Hamming metric
    util::bitvec bv;
    bv.set_parameters(3, 27);
    bv.set(0,13); bv.set(1,14); bv.set(2,15);
    return bv.norm() == 42;
}
bool bitvec_equality() {
    bool equality_works = true;
    util::bitvec bv, bv2;
    bv.set_parameters(3, 27);
    bv2.set_parameters(3, 27);
    
    equality_works &= (bv == bv2);
    
    bv.set(0,3);
    bv.set(1,4);
    bv.set(2,27);
    bv2.set(1,3);
    
    equality_works &= !(bv == bv2);
    
    bv.set(0,3);
    bv.set(1,4);
    bv.set(2,27);
    bv2.set(0,3);
    bv2.set(1,4);
    bv2.set(2,27);
    
    equality_works &= (bv == bv2);
    
    return equality_works;
}

bool bitvec_dist(){ // Hamming metric
    util::bitvec bv, bv2;
    bv.set_parameters(3, 27);
    bv2.set_parameters(3, 27);
    
    bv.set(0,27); bv.set(1,0); bv.set(2,0);
    bv2.set(0,0); bv2.set(1,27); bv2.set(2,0);
    int difference = util::dist(bv, bv2);
    return (difference == 54);
}

// helper method for L1 distance
uint32_t abs_dist(uint32_t v1, uint32_t v2){
    if( v1 >= v2 ){ return v1 - v2; }
    return v2 - v1;
}
uint32_t l1_dist(const std::vector<uint32_t>& p1,
                 const std::vector<uint32_t>& p2)
{
    assert(p1.size() == p2.size());
    uint32_t dist = 0;
    for(int i = 0; i < p1.size(); ++i){
        dist += abs_dist(p1[i],p2[i]);
    }
    return dist;
}

// tests for the classicalLSH r-near DS
bool classicalLSH_find_1D_approx_rnear() {
    
    using point_t = typename lsh::point_t;
    using tuple_t = typename lsh::tuple_t;
    
    // build the dataset
    uint32_t d = 1, U = 128;
    std::vector<point_t> dataset(U);
    for(uint32_t i = 0; i < U; ++i){
        dataset[i].push_back(i);
    }
    
    // choose the parameters we care about
    uint32_t r = 10, c = 2, k = 15;
    
    // build the LSH for r-near
    classical::near_lsh rLSH;
    rLSH.set_parameters(r,d,U,c);
    rLSH.set_dataset(dataset);
    rLSH.build();
    
    // query
    std::set<tuple_t> result_indices;
    point_t q { 30 };
    int num_found = rLSH.k_near(q, k, result_indices);
    std::cout << "Number matches found are: " << num_found << std::endl;
    
    // print actual errors
    for(tuple_t tuple : result_indices){
        auto idx = tuple.id;
        auto d_  = tuple.dist_val;
        std::cout << "dist(p[" << idx << "], q) = " << l1_dist(dataset[idx], q) << " vs dist = " << d_ << std::endl;
    }
    
    return (num_found >= k);
}

// tests for the coveringLSH r-near DS
bool coveringLSH_find_1D_approx_rnear() {
    using point_t = typename lsh::point_t;
    using tuple_t = typename lsh::tuple_t;
    
    // build the dataset
    uint32_t d = 1, U = 128;
    std::vector<point_t> dataset(U);
    for(uint32_t i = 0; i < U; ++i){
        dataset[i].push_back(i);
    }
    
    // choose the parameters we care about
    uint32_t r = 10, c = 2, k = 15;
    
    // build the LSH for r-near
    covering::near_lsh rLSH;
    rLSH.set_parameters(r,d,U,c);
    rLSH.set_dataset(dataset);
    rLSH.build();
    
    // query
    std::set<tuple_t> result_indices;
    point_t q { 30 };
    int num_found = rLSH.k_near(q, k, result_indices);
    std::cout << "Number matches found are: " << num_found << std::endl;
    
    // print actual errors
    for(tuple_t tuple : result_indices){
        auto idx = tuple.id;
        auto d_  = tuple.dist_val;
        std::cout << "dist(p[" << idx << "], q) = " << l1_dist(dataset[idx], q) << " vs dist = " << d_ << std::endl;
    }
    
    return (num_found >= k);
}

#define TEST(func) test_map[ #func ] = func ;

// run the tests
void run_tests() {
    typedef bool (*test_func_t)();
    
    // test map
    std::map<std::string, test_func_t> test_map;
    
    // add the tests
    TEST(bitvec_and);
    TEST(bitvec_xor);
    TEST(bitvec_dotproduct);
    TEST(bitvec_norm);
    TEST(bitvec_equality);
    TEST(bitvec_dist);
    
    // run the tests
    size_t num_success = 0;
    for(auto& it: test_map){
        bool succeed = (*it.second)();
        if( !succeed ){
            std::cout << "[X] Test (" << it.first <<") has failed." << std::endl;
        }else{
            std::cout << "[✓] Test (" << it.first <<") has succeeded." << std::endl;
            ++num_success;
        }
    }
    
    // print final statistics
    std::cout << "\n\nThere have been:\n" << num_success << " successes\n" << test_map.size() - num_success << " failures" << std::endl;
    
    
}

}
