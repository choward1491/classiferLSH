//
//  cjh_hashtable.hpp
//  classifierLSH
//
//  Created by Christian Howard on 12/17/19.
//  Copyright © 2019 Christian Howard. All rights reserved.
//

#ifndef cjh_hashtable_hpp
#define cjh_hashtable_hpp

#include <vector>
#include <list>

namespace cjh {

// class to represent a simple hash table because it seems
// the std::unordered_map is doing something weird with
// buckets and collisions. Let's see how this does
template<   typename key_type,
            typename data_type,
            typename hashfunc >
class hashtable {
public:
    
    // define useful type
    struct node {
        node() = default;
        node(size_t hval, const key_type& k):hash(hval){
            key = k;
        }
        size_t      hash;
        key_type    key;
        data_type   data;
    };
    
    // ctor/dtor
    hashtable() = default;
    hashtable(size_t n, const hashfunc& hf_):table(n) {
        hf = hf_;
    }
    ~hashtable() = default;
    
    // getter/setters
    void resize(size_t num_buckets) {
        // do later
    }
    void set_hashfunc(const hashfunc& hf_) {
        hf = hf_;
    }
    hashfunc& get_hashfunc() { return hf; }
    const hashfunc& get_hashfunc() const { return hf; }
    size_t size() const { return num_data; }
    size_t num_buckets() const { return table.size(); }
    
    // useful methods
    data_type& operator[](const key_type& k) {
        auto& list = get_bucket_list(k);
        auto it = list.cbegin(), end = list.cend();
        
        // loop over bucket elements
        for(; it != end; ++it){
            if( it->key == k ){
                return it->data;
            }
        }// end for
        
        // if did not find the element, add it
        auto& newnode = list.emplace_front(hf(k), k);
        
        // return the new data
        return newnode.data;
    }
    const data_type& operator[](const key_type& k, bool* success = nullptr) const {
        auto& list = get_bucket_list(k);
        auto it = list.cbegin(), end = list.cend();
        
        // loop over bucket elements
        for(; it != end; ++it){
            if( it->key == k ){
                return it->data;
            }
        }// end for
        
        if( success ){
            *success = false;
        }
        return dummy;
        
    }
    const std::list<struct node>& get_bucket_list(const key_type& k) const {
        size_t hash = hf(k);
        return table[hash % table.size()];
    }
    
private:
    
    size_t num_data;
    hash_func hf;
    data_type dummy;
    std::vector<std::list<struct node>> table;
    
    // useful private methods
    std::list<struct node>& get_bucket_list(const key_type& k) {
        size_t hash = hf(k);
        return table[hash % table.size()];
    }
    
};


}

#endif /* cjh_hashtable_h */
