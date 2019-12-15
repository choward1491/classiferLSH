//
//  bitvector.hpp
//  classifierLSH
//
//  Created by Christian Howard on 12/13/19.
//  Copyright © 2019 Christian Howard. All rights reserved.
//

#ifndef bitvector_hpp
#define bitvector_hpp

#include <array>
#include <vector>

namespace util {
    
    class bitvec {
    public:
        
        // ctor
        bitvec();
        
        // useful operations
        void set(int dim, uint32_t v);
        void set_int(int idx, uint32_t v);
        uint32_t get(int dim) const;
        uint32_t num_dims() const;
        uint32_t univ_size() const;
        size_t num_bits() const;
        size_t num_compressed_ints() const;
        void set_parameters(uint32_t dims, uint32_t Usize);
        
        // useful operators and operations
        bitvec operator&( const bitvec& bv ) const;
        bitvec operator^( const bitvec& bv ) const;
        bool operator==(const bitvec& bv) const;
        uint32_t dot( const bitvec& bv ) const;
        uint32_t norm() const;
        
    private:
        
        uint32_t d, U, NUM_INTS, map[256];
        std::vector<uint32_t> bits;
        
        // helper methods
        void build_map();
        uint32_t num_nonzero_bits8(uint8_t n) const;
        uint32_t num_nonzero_bits32(uint32_t n) const;
        
    };


    uint32_t dist(const bitvec& bv1, const bitvec& bv2);
    
}

#endif /* bitvector_h */
