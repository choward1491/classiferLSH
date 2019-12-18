//
//  lsh_defn.hpp
//  classifierLSH
//
//  Created by Christian Howard on 12/15/19.
//  Copyright © 2019 Christian Howard. All rights reserved.
//

#ifndef lsh_defn_hpp
#define lsh_defn_hpp

#include "bitvector.hpp"
#include <unordered_map>
#include "cjh_hashtable.hpp"

namespace lsh {

    // definition for a point type
    using point_t = std::vector<uint32_t>;

    // definition for a (distance,ID) output tuple
    struct tuple_t {
        uint32_t dist_val;
        size_t id;
        bool operator<(const tuple_t& v) const {
            if( dist_val < v.dist_val ){ return true; }
            else if( dist_val > v.dist_val ){ return false; }
            else{ // dist_val == v.dist_val
                return (id < v.id);
            }
        }
        bool operator==(const tuple_t& v) const {
            return (dist_val == v.dist_val) && (id == v.id);
        }
    };

    // definition for a hash function
    struct hash_func {
        hash_func() = default;
        hash_func(const util::bitvec& bv_):bv(bv_){}
        size_t operator()( const util::bitvec& v) const {
            auto and_bv = (v & bv);
            auto num_ints = and_bv.num_compressed_ints();
            size_t hash = 1654033;
            for(size_t i = 0; i < num_ints; ++i){
              hash ^= and_bv.get_int(i) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            }
            return hash;
        }
        
        util::bitvec bv;
    };

    // definition for a hash table type for the LSH purposes
    using hashtable_t = std::unordered_map<
                                    util::bitvec,
                                    size_t,
                                    hash_func
                                    >;

    /*using hashtable_t = cjh::hashtable<
                                    util::bitvec,
                                    size_t,
                                    hash_func
                                    >;*/

}

#endif /* lsh_defn_h */
