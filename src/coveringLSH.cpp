//
//  coveringLSH.cpp
//  classifierLSH
//
//  Created by Christian Howard on 12/14/19.
//  Copyright © 2019 Christian Howard. All rights reserved.
//

#include "coveringLSH.hpp"

namespace covering {

    void near_lsh::set_parameters(uint32_t r_, uint32_t d_, uint32_t U_ , uint32_t c_) {
        r = r_; d = d_; U = U_; c = c_;
    }

    void near_lsh::set_dataset( std::vector<lsh::point_t>& dset ){
        
        // set the pointer to the main dataset
        ds_ref = &dset;
        
        // construct the bitvector version of dataset
        bit_dataset.resize(dset.size());
        for(int i = 0; i < dset.size(); ++i){
            orig2bitvec(dset[i], bit_dataset[i]);
        }// end loop over data vectors
    }

    // build the hash table
    void near_lsh::build() {
        
        // compute nunmber of hash functions and corresponding hash tables
        uint32_t L = 20 * static_cast<uint32_t>(std::pow(ds_ref->size(), 1.0/c));
        
        // resize the hash table list
        hash_tables.clear();
        hash_tables.resize(L);
        
        // choose seed
        uint32_t seed = 123456789;
        uint32_t m = -1, a = 1103515245, b = 12345;
        
        // dataset size
        size_t N = ds_ref->size();
        
        // build each hash table using random bitvector construction
        //#pragma omp parallel for
        for(uint32_t i = 0; i < L; ++i){
            
            // define new random hash function
            lsh::hash_func hf; hf.bv.set_parameters(d, U);
            for(size_t j = 0; j < hf.bv.num_compressed_ints(); ++j){
                //#pragma omp atomic
                seed += (( a * seed + b ) % m - seed);
                hf.bv.set_int(j, seed );
            }
            
            // get the ith hash table
            auto& ht = hash_tables[i];
            ht = lsh::hashtable_t(N, hf);
            
            // hash dataset into the ith hash table
            for(size_t j = 0; j < N; ++j){
                ht[bit_dataset[j]] = j;
            }
        }// loop number of hash functions
        
    }

    int near_lsh::k_near( const lsh::point_t& q, int k, std::set<lsh::tuple_t>& nearest_ind ) const {
        
        // get a bitvec version of query point
        util::bitvec qv;
        orig2bitvec(q, qv);
        
        // counter keeping track of number found
        int num_found = 0;
        nearest_ind.clear();
        
        // loop over hash tables until we find k near points or fail to
        size_t num_ht = hash_tables.size();
        #pragma omp parallel for
        for(size_t i = 0; i < num_ht; ++i){
            
            // only look for new values if we have not
            // yet exceeded the number of values we hope to find
            if(num_found < k){
                
                // get the ith hash table
                const auto& ht = hash_tables[i];
                
                // get the bucket index for the query point
                auto bucket_id = ht.bucket(qv);
                
                // iterate over the bucket elements
                auto it = ht.cbegin(bucket_id);
                auto end= ht.cend(bucket_id);
                for(; it != end; ++it){
                    
                    // compute distance between elements
                    size_t idx = it->second;
                    uint32_t dist_ = util::dist(qv, bit_dataset[idx]);
                    
                    // if the distance is below our approximation
                    // radius, then we will add its index
                    // to our output result
                    if( dist_ <= c*r ){ // satisfying element
                        
                        // only one thread should push to
                        // the vector at once and similarly
                        // for updating the counter
                        #pragma omp critical
                        {
                            if( num_found < k ){
                                lsh::tuple_t output;
                                output.dist_val = dist_;
                                output.id = idx;
                                nearest_ind.insert(output);
                                num_found = nearest_ind.size();
                            }
                        }
                        
                        if( num_found >= k ){ break; }
                        
                    }
                    
                }// end loop over bucket elements
            }
            
        }// end for i
        
        // return the number of found satisfying points
        return num_found;
    }

    void near_lsh::orig2bitvec(const lsh::point_t& p, util::bitvec& bv) const {
        bv.set_parameters(d,U);
        for(int j = 0; j < d; ++j){
            bv.set(j, p[j]);
        }// loop over vector components
    }

}
