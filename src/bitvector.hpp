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
    
// class that represents the bitvector transformation
// for some vector v \in [U]^d for some universe size U
// and some dimension d. This representation allows for
// computing the necessary operations using a bitvector
// in Hamming space, like bitwise AND and XOR, as well as
// norms wrt Hamming measure. Can also set bits based on
// values in vector in original universe space
    class bitvec {
    public:
        
        // ctor
        bitvec();
        ~bitvec() = default;
        
        // useful operations
        void        set(int dim, uint32_t v);
        void        set_int(int idx, uint32_t v);
        uint32_t    get(int dim) const;
        uint32_t    get_int(int idx) const;
        uint32_t    num_dims() const;
        uint32_t    univ_size() const; // universe size
        size_t      num_bits() const;
        size_t      num_compressed_ints() const;
        void        set_parameters(uint32_t dims, uint32_t Usize);
        
        // useful operators and operations
        bitvec      operator&( const bitvec& bv ) const;
        bitvec      operator^( const bitvec& bv ) const;
        bool        operator==(const bitvec& bv) const;
        uint32_t    dot( const bitvec& bv ) const;
        uint32_t    norm() const;
        
    private:
        
        uint32_t d, U, NUM_INTS, map[256];
        std::vector<uint32_t> bits;
        
        // helper methods
        void        build_map();
        uint32_t    num_nonzero_bits8(uint8_t n) const;
        uint32_t    num_nonzero_bits32(uint32_t n) const;
        
    };

    // Hamming distance between two vectors
    uint32_t dist(const bitvec& bv1, const bitvec& bv2);
    
}

#endif /* bitvector_h */
