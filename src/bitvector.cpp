//
//  bitvector.cpp
//  classifierLSH
//
//  Created by Christian Howard on 12/13/19.
//  Copyright © 2019 Christian Howard. All rights reserved.
//

#include "bitvector.hpp"

namespace util {
    
    // ctor
    bitvec::bitvec(){ d = 1; U = 32; NUM_INTS = (U*d)/32 + 1; build_map(); }
    
    // useful operations
    void bitvec::set(int dim, uint32_t v) {
        int i1 = (U*dim)/32, d1 = (U*dim) % 32;
        int is = i1, ds = d1;
        uint32_t max = -1;
        for(uint32_t u = 0; u < U; ++u){
            if( u < v ){
                bits[is] = bits[is] | (1 << ds);
            }else{
                bits[is] = bits[is] & ((1 << ds) ^ max);
            }
            /*if( is >= NUM_INTS ){
                printf("exceeding!\n");
            }*/
            ++ds;
            if( ds == 32 ){ ds = ds % 32; ++is; }
        }
    }
    void bitvec::set_int(int idx, uint32_t v) {
        bits[idx] = v;
    }
    uint32_t bitvec::get_int(int idx) const {
        return bits[idx];
    }

    void bitvec::set_parameters(uint32_t dims, uint32_t univ_size) {
        d = dims;
        U = univ_size;
        NUM_INTS = (U*d)/32 + 1;
        bits.resize(NUM_INTS,0);
    }

    uint32_t bitvec::get(int dim) const {
        uint32_t out = 0;
        int i1 = (U*dim)/32, d1 = (U*dim) % 32;
        int is = i1, ds = d1;
        for(uint32_t u = 0; u < U; ++u){
            if( bits[is] & (1 << ds) ){ ++out; }
            ++ds;
            if( ds == 32 ){ ds = ds % 32; ++is; }
        }
        
        return out;
    }

    uint32_t bitvec::num_dims() const {
        return d;
    }
    uint32_t bitvec::univ_size() const {
        return U;
    }

    size_t bitvec::num_bits() const {
        return d * U;
    }
    size_t bitvec::num_compressed_ints() const {
        return NUM_INTS;
    }
    
    bitvec bitvec::operator&( const bitvec& bv ) const {
        assert(bv.num_compressed_ints() == num_compressed_ints());
        bitvec output; output.set_parameters(d, U);
        for(int i = 0; i < NUM_INTS; ++i){
            output.bits[i] = bits[i] & bv.bits[i];
        }
        uint32_t mask = -1;
        mask = (mask >> (32*NUM_INTS - d*U));
        output.bits[NUM_INTS-1] &= mask;
        return output;
    }
    
    bitvec bitvec::operator^( const bitvec& bv ) const {
        assert(bv.num_compressed_ints() == num_compressed_ints());
        bitvec output; output.set_parameters(d, U);
        for(int i = 0; i < NUM_INTS; ++i){
            output.bits[i] = bits[i] ^ bv.bits[i];
        }
        uint32_t mask = -1;
        mask = (mask >> (32*NUM_INTS - d*U));
        output.bits[NUM_INTS-1] &= mask;
        return output;
    }

    bool bitvec::operator==(const bitvec& bv) const {
        assert(bv.num_compressed_ints() == num_compressed_ints());
        bool output = true;
        for(int i = 0; i < NUM_INTS; ++i){
            output &= (bits[i] == bv.bits[i]);
        }
        return output;
    }
    
    uint32_t bitvec::dot( const bitvec& bv ) const {
        uint32_t dp = 0;
        for(int i = 0; i < NUM_INTS; ++i){
            dp += num_nonzero_bits32(bv.bits[i] & bits[i]);
        }
        return dp;
    }
    
    uint32_t bitvec::norm() const {
        return dot(*this);
    }
    
    void bitvec::build_map() {
        for(uint32_t i = 0; i < 256; ++i){
            map[i] = 0;
            for(int j = 0; j < 8; ++j){
                map[i] += ((i & (1<<j)) > 0);
            }// compute num nonzero bits for number i
        }// loop over number i
    }
    
    uint32_t bitvec::num_nonzero_bits8(uint8_t n) const {
        return map[n];
    }
    uint32_t bitvec::num_nonzero_bits32(uint32_t n) const {
        uint32_t nbits = 0;
        uint8_t* ptr = (uint8_t*)&n;
        for(int i = 0; i < 4; ++i){
            nbits += num_nonzero_bits8(*ptr);
            ptr = ptr + 1;
        }
        return nbits;
    }

    uint32_t dist(const bitvec& bv1, const bitvec& bv2){
        return (bv1 ^ bv2).norm();
    }
    
}
