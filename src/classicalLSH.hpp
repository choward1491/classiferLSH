//
//  classicalLSH.hpp
//  classifierLSH
//
//  Created by Christian Howard on 12/14/19.
//  Copyright © 2019 Christian Howard. All rights reserved.
//

#ifndef classicalLSH_hpp
#define classicalLSH_hpp

#include <vector>
#include <set>
#include "lsh_defn.hpp"

namespace classical {

    // data structure to approximately solve r-near problem
    // for some input dataset that is a subset of [U]^d for some
    // universe size U and some dimension d. This data structure
    // uses the classical LSH approach via i.i.d. randomized hash functions
    // and making enough corresponding hash tables to make the error probability
    // extremely small. Note that data structure makes use of threading to
    // improve the construction and querying of the data structure.
    class near_lsh {
    public:
        
        // ctor/dtor
        near_lsh() = default;
        ~near_lsh() = default;
        
        // set the parameters for the data
        void set_parameters(uint32_t r,
                            uint32_t d,
                            uint32_t U,
                            uint32_t c = 2);
        
        // set the dataset
        void set_dataset( std::vector<lsh::point_t>& ds );
        
        // build the hash table
        void build();
        
        // query methods
        int k_near( const lsh::point_t& q,
                   int k,
                   std::set<lsh::tuple_t>& nearest_ind ) const;
        
    private:
        
        // internal state
        uint32_t                        c, r, d, U;
        std::vector<lsh::point_t>       *ds_ref;
        std::vector<util::bitvec>       bit_dataset;
        std::vector<lsh::hashtable_t>   hash_tables;
        
        // useful helper functions
        void orig2bitvec(const lsh::point_t& p, util::bitvec& bv) const;
        
    };

}

#endif /* classicalLSH_hpp */
