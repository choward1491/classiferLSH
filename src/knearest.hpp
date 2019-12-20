//
//  knearest.hpp
//  classifierLSH
//
//  Created by Christian Howard on 12/19/19.
//  Copyright © 2019 Christian Howard. All rights reserved.
//

#ifndef knearest_hpp
#define knearest_hpp

#include <vector>
#include <algorithm>

namespace approx {

    template<typename r_near_DS>
    class knearest {
    public:
        
        // ctor/dtor
        knearest() = default;
        ~knearest() = default;
        
        void set_parameters(const std::vector<uint32_t>& radii_,
                            uint32_t d_,
                            uint32_t U_,
                            uint32_t c_ = 2)
        {
            radii   = radii_;
            std::sort(radii.begin(), radii.end());
            d       = d_;
            U       = U_;
            c       = c_;
        }
        
        void set_dataset( std::vector<lsh::point_t>& ds ) {
            ds_ref = &ds;
        }
        
        // build the data structure
        void build(){
            
            // resize the rnear Data Structure instances
            rnear_ds_instances.resize(radii.size());
            
            // loop over the radii and construct the associated DS
            for(size_t idx = 0; idx < radii.size(); ++idx){
                auto& rnear = rnear_ds_instances[idx];
                rnear.set_parameters(radii[idx], d, U, c);
                rnear.set_dataset(*ds_ref);
                rnear.build();
            }
        }
        
        // k nearest neighbor
        int k_nearest(const lsh::point_t& q,
                      int k,
                      std::set<lsh::tuple_t>& nearest_ind ) const
        {
            // clear the nearest neighbor size
            nearest_ind.clear();
            
            // loop over data structures for fixed radii until
            // we obtain k points near our given point point
            for(size_t idx = 0; idx < rnear_ds_instances.size(); ++idx){
                auto num = rnear_ds_instances[idx].k_near(q, k, nearest_ind);
                if( num == k ){ break; }
            }
            
            // return number of nearest neighbor points that were found
            return static_cast<int>(nearest_ind.size());
        }
        
        
    private:
        
        // internal state
        uint32_t                        c, d, U;
        std::vector<lsh::point_t>       *ds_ref;
        std::vector<uint32_t>           radii;
        std::vector<r_near_DS>          rnear_ds_instances;
        
    };


}


#endif /* knearest_h */
