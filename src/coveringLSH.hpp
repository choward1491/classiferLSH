//
//  coveringLSH.hpp
//  classifierLSH
//
//  Created by Christian Howard on 12/14/19.
//  Copyright © 2019 Christian Howard. All rights reserved.
//

#ifndef coveringLSH_hpp
#define coveringLSH_hpp

#include <vector>
#include <set>
#include "lsh_defn.hpp"

namespace covering {

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
        void set_v_vector(util::bitvec& bv, uint32_t vidx) const;
        
    };

}

#endif /* coveringLSH_hpp */
