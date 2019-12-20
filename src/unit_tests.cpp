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
#include "knearest.hpp"
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include "omp.h"
#include <stdio.h>

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
bool classicalLSH_find_1D_approx_rnear(int num_threads, uint32_t U) {
    
    using point_t = typename lsh::point_t;
    using tuple_t = typename lsh::tuple_t;
    
    // set the number of threads
    omp_set_num_threads(num_threads);
    
    // build the dataset
    uint32_t d = 1;
    std::vector<point_t> dataset(U+1);
    for(uint32_t i = 0; i <= U; ++i){
        dataset[i].push_back(i);
    }
    
    // choose the parameters we care about
    uint32_t r = 5, c = 2, k = r;
    
    // build the LSH for r-near
    classical::near_lsh rLSH;
    rLSH.set_parameters(r,d,U,c);
    rLSH.set_dataset(dataset);
    
    double t1 = omp_get_wtime();
    rLSH.build();
    double t2 = omp_get_wtime();
    printf("Time to build classicalLSH with %i threads was %0.5e seconds\n",
           omp_get_max_threads(), t2 - t1);
    
    // query
    bool is_correct = true;
    std::set<tuple_t> result_indices;
    point_t q { 100 };
    
    double t3 = omp_get_wtime();
    int num_found = rLSH.k_near(q, k, result_indices);
    double t4 = omp_get_wtime();
    printf("Time to query classicalLSH once with %i threads was %0.5e seconds\n",
           omp_get_max_threads(), t4 - t3);
    is_correct = (num_found >= k);
    std::cout << "Number matches found are: " << num_found << std::endl;
    
    
    // print actual errors
    for(tuple_t tuple : result_indices){
        auto idx = tuple.id;
        auto d_  = tuple.dist_val;
        is_correct = is_correct && (l1_dist(dataset[idx], q) == d_);
        //std::cout << "dist(p[" << idx << "], q) = " << l1_dist(dataset[idx], q) << " vs dist = " << d_ << std::endl;
    }
    
    return is_correct;
}

bool classicalLSH_find_1D_approx_knearest(int num_threads, int k, uint32_t U) {
    using point_t = typename lsh::point_t;
    using tuple_t = typename lsh::tuple_t;
    using knearest_t = typename approx::knearest<classical::near_lsh>;
    
    // set the number of threads
    omp_set_num_threads(num_threads);
    
    // build the dataset
    uint32_t d = 1;
    std::vector<point_t> dataset(U+1);
    for(uint32_t i = 0; i <= U; ++i){
        dataset[i].push_back(i);
    }
    
    // choose the parameters we care about
    uint32_t c = 2;
    std::vector<uint32_t> radii { 1, 2, 3, 4, 5, 6 };
    
    // build the LSH for r-near
    knearest_t knn_ds;
    knn_ds.set_parameters(radii, d, U, c);
    knn_ds.set_dataset(dataset);
    
    double t1 = omp_get_wtime();
    knn_ds.build();
    double t2 = omp_get_wtime();
    printf("Time to build knn based on classicalLSH with %i threads was %0.5e seconds\n",
           omp_get_max_threads(), t2 - t1);
    
    // query
    bool is_correct = true;
    std::set<tuple_t> result_indices;
    point_t q { 100 };
    
    double t3 = omp_get_wtime();
    int num_found = knn_ds.k_nearest(q, k, result_indices);
    double t4 = omp_get_wtime();
    printf("Time to query knn based on classicalLSH once with %i threads was %0.5e seconds\n",
           omp_get_max_threads(), t4 - t3);
    is_correct = (num_found >= k);
    std::cout << "Number matches found are: " << num_found << std::endl;
    
    
    // print actual errors
    for(tuple_t tuple : result_indices){
        auto idx = tuple.id;
        auto d_  = tuple.dist_val;
        is_correct = is_correct && (l1_dist(dataset[idx], q) == d_);
        std::cout << "dist(p[" << idx << "], q) = " << l1_dist(dataset[idx], q) << " vs dist = " << d_ << std::endl;
    }
    
    return is_correct;
}

// tests for the coveringLSH r-near DS
bool coveringLSH_find_1D_approx_rnear(int num_threads, uint32_t U) {
    using point_t = typename lsh::point_t;
    using tuple_t = typename lsh::tuple_t;
    
    // set the number of threads
    omp_set_num_threads(num_threads);
    
    // build the dataset
    uint32_t d = 1;
    std::vector<point_t> dataset(U+1);
    for(uint32_t i = 0; i <= U; ++i){
        dataset[i].push_back(i);
    }
    
    // choose the parameters we care about
    uint32_t r = 5, c = 2, k = r;
    
    // build the LSH for r-near
    covering::near_lsh rLSH;
    rLSH.set_parameters(r,d,U,c);
    rLSH.set_dataset(dataset);
    double t1 = omp_get_wtime();
    rLSH.build();
    double t2 = omp_get_wtime();
    printf("Time to build coveringLSH with %i threads was %0.5e seconds\n",
           omp_get_max_threads(), t2 - t1);
    
    // query
    std::set<tuple_t> result_indices;
    point_t q { 100 };
    double t3 = omp_get_wtime();
    int num_found = rLSH.k_near(q, k, result_indices);
    double t4 = omp_get_wtime();
    printf("Time to query coveringLSH once with %i threads was %0.5e seconds\n",
           omp_get_max_threads(), t4 - t3);
    bool is_correct = (num_found >= k);
    std::cout << "Number matches found are: " << num_found << std::endl;
    
    // print actual errors
    for(tuple_t tuple : result_indices){
        auto idx = tuple.id;
        auto d_  = tuple.dist_val;
        is_correct = is_correct && (l1_dist(dataset[idx], q) == d_);
        //std::cout << "dist(p[" << idx << "], q) = " << l1_dist(dataset[idx], q) << " vs dist = " << d_ << std::endl;
    }
    
    return is_correct;
}

bool coveringLSH_find_1D_approx_knearest(int num_threads, int k, uint32_t U) {
    using point_t = typename lsh::point_t;
    using tuple_t = typename lsh::tuple_t;
    using knearest_t = typename approx::knearest<covering::near_lsh>;
    
    // set the number of threads
    omp_set_num_threads(num_threads);
    
    // build the dataset
    uint32_t d = 1;
    std::vector<point_t> dataset(U+1);
    for(uint32_t i = 0; i <= U; ++i){
        dataset[i].push_back(i);
    }
    
    // choose the parameters we care about
    uint32_t c = 2;
    std::vector<uint32_t> radii { 1, 2, 3, 4, 5, 6 };
    
    // build the LSH for r-near
    knearest_t knn_ds;
    knn_ds.set_parameters(radii, d, U, c);
    knn_ds.set_dataset(dataset);
    
    double t1 = omp_get_wtime();
    knn_ds.build();
    double t2 = omp_get_wtime();
    printf("Time to build knn based on coveringLSH with %i threads was %0.5e seconds\n",
           omp_get_max_threads(), t2 - t1);
    
    // query
    bool is_correct = true;
    std::set<tuple_t> result_indices;
    point_t q { 100 };
    
    double t3 = omp_get_wtime();
    int num_found = knn_ds.k_nearest(q, k, result_indices);
    double t4 = omp_get_wtime();
    printf("Time to query knn based on coveringLSH once with %i threads was %0.5e seconds\n",
           omp_get_max_threads(), t4 - t3);
    is_correct = (num_found >= k);
    std::cout << "Number matches found are: " << num_found << std::endl;
    
    
    // print actual errors
    for(tuple_t tuple : result_indices){
        auto idx = tuple.id;
        auto d_  = tuple.dist_val;
        is_correct = is_correct && (l1_dist(dataset[idx], q) == d_);
        std::cout << "dist(p[" << idx << "], q) = " << l1_dist(dataset[idx], q) << " vs dist = " << d_ << std::endl;
    }
    
    return is_correct;
}

void simple_2D_binary_classification_classical(int num_threads, int k, uint32_t U) {
    using point_t = typename lsh::point_t;
    using tuple_t = typename lsh::tuple_t;
    using knearest_t = typename approx::knearest<classical::near_lsh>;
    
    // set the number of threads
    omp_set_num_threads(num_threads);
    
    // random number stuff
    uint32_t seed = 123456789;
    uint32_t m = -1, a = 1103515245, b = 12345;
    
    // build the dataset
    uint32_t n = (U*U)/6;
    uint32_t d = 2;
    std::set<point_t> sample_points;
    for(uint32_t i = 0; i < n; ++i){
        seed = ( a * seed + b ) % m;
        point_t pt(2);
        pt[0] = seed % (U+1);
        seed = ( a * seed + b ) % m;
        pt[1] = seed % (U+1);
        sample_points.insert(pt);
    }
    
    std::vector<point_t>    dataset(sample_points.size());
    std::vector<uint32_t>   labels(sample_points.size());
    point_t center(2, U/2);
    {
        uint32_t ii = 0;
        for(auto& pt: sample_points){
            dataset[ii] = pt;
            labels[ii++] = l1_dist(center,pt) > (U / 5);
        }
    }
    
    // choose the parameters we care about
    uint32_t c = 2;
    std::vector<uint32_t> radii { 1, 2, 3, 4, 5, 6 };
    
    // build the LSH for r-near
    knearest_t knn_ds;
    knn_ds.set_parameters(radii, d, U, c);
    knn_ds.set_dataset(dataset);
    
    double t1 = omp_get_wtime();
    knn_ds.build();
    double t2 = omp_get_wtime();
    printf("Time to build knn based on classicalLSH with %i threads was %0.5e seconds\n",
           omp_get_max_threads(), t2 - t1);
    
    // query
    std::set<tuple_t> result_indices;
    point_t q(2);
    FILE* file = fopen("/Users/cjh/Documents/code/classiferLSH/output/classic_circle_class2.csv", "w");
    
    if( file ){
        // try to see how the classifier performs
        for(uint32_t i = 0; i <= U; ++i){
            for(uint32_t j=0; j <= U; ++j){
                q[0] = i; q[1] = j;
                uint32_t label = 0, true_label = l1_dist(center, q) > (U / 5);
                result_indices.clear();
                knn_ds.k_nearest(q, k, result_indices);
                
                // try to figure out what to label the result
                if( result_indices.size()){
                    int num_0 = 0, num_1 = 1;
                    for(auto tuple: result_indices ){
                        if( labels[tuple.id] ){ ++num_1; }
                        else{ ++num_0; }
                    }
                    if( num_1 >= num_0 ){ label = 1; }
                }
                
                // print the result
                fprintf(file, "%u, %u, %u, %u\n", i, j, true_label, label);
                //std::cout << i << " " << j << " " << true_label << " " << label << std::endl;
                
            }// end for j
        }// end for i
        
        fclose(file); file = nullptr;
    }else{
        std::cout << "Could not open the file!" << std::endl;
    }
}
void simple_2D_binary_classification_covering(int num_threads, int k, uint32_t U) {
    using point_t = typename lsh::point_t;
    using tuple_t = typename lsh::tuple_t;
    using knearest_t = typename approx::knearest<covering::near_lsh>;
    
    // set the number of threads
    omp_set_num_threads(num_threads);
    
    // random number stuff
    uint32_t seed = 123456789;
    uint32_t m = -1, a = 1103515245, b = 12345;
    
    // build the dataset
    uint32_t n = (U*U)/6;
    uint32_t d = 2;
    std::set<point_t> sample_points;
    for(uint32_t i = 0; i < n; ++i){
        seed = ( a * seed + b ) % m;
        point_t pt(2);
        pt[0] = seed % (U+1);
        seed = ( a * seed + b ) % m;
        pt[1] = seed % (U+1);
        sample_points.insert(pt);
    }
    
    std::vector<point_t>    dataset(sample_points.size());
    std::vector<uint32_t>   labels(sample_points.size());
    point_t center(2, U/2);
    {
        uint32_t ii = 0;
        for(auto& pt: sample_points){
            dataset[ii] = pt;
            labels[ii++] = l1_dist(center,pt) > (U / 5);
        }
    }
    
    // choose the parameters we care about
    uint32_t c = 2;
    std::vector<uint32_t> radii { 1, 2, 3, 4, 5, 6 };
    
    // build the LSH for r-near
    knearest_t knn_ds;
    knn_ds.set_parameters(radii, d, U, c);
    knn_ds.set_dataset(dataset);
    
    double t1 = omp_get_wtime();
    knn_ds.build();
    double t2 = omp_get_wtime();
    printf("Time to build knn based on coveringLSH with %i threads was %0.5e seconds\n",
           omp_get_max_threads(), t2 - t1);
    
    // query
    std::set<tuple_t> result_indices;
    point_t q(2);
    FILE* file = fopen("/Users/cjh/Documents/code/classiferLSH/output/covering_circle_class2.csv", "w");
    
    if( file ){
        // try to see how the classifier performs
        for(uint32_t i = 0; i <= U; ++i){
            for(uint32_t j=0; j <= U; ++j){
                q[0] = i; q[1] = j;
                uint32_t label = 0, true_label = l1_dist(center, q) > (U / 5);
                result_indices.clear();
                knn_ds.k_nearest(q, k, result_indices);
                
                // try to figure out what to label the result
                if( result_indices.size()){
                    int num_0 = 0, num_1 = 1;
                    for(auto tuple: result_indices ){
                        if( labels[tuple.id] ){ ++num_1; }
                        else{ ++num_0; }
                    }
                    if( num_1 >= num_0 ){ label = 1; }
                }
                
                // print the result
                fprintf(file, "%u, %u, %u, %u\n", i, j, true_label, label);
                //std::cout << i << " " << j << " " << true_label << " " << label << std::endl;
                
            }// end for j
        }// end for i
        
        fclose(file); file = nullptr;
    }else{
        std::cout << "Could not open the file!" << std::endl;
    }
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
    //TEST(classicalLSH_find_1D_approx_rnear);
    //TEST(coveringLSH_find_1D_approx_rnear);
    
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
