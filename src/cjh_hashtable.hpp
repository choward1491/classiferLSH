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
            typename hash_func >
class hashtable {
public:
    
    // define useful type
    using list_t = std::list<struct node>;
    struct node {
        node() = default;
        node(size_t hval, const key_type& k):hash(hval){
            first = k;
        }
        size_t      hash;
        key_type    first;
        data_type   second;
    };
    
    // ctor/dtor
    hashtable() = default;
    hashtable(size_t n, const hash_func& hf_):table(n) {
        hf = hf_;
    }
    ~hashtable() = default;
    
    // getter/setters
    void resize(size_t num_buckets) {
        // do later
    }
    void set_hashfunc(const hash_func& hf_) {
        hf = hf_;
    }
    hash_func& get_hashfunc() { return hf; }
    hash_func& hash_function() { return hf; }
    const hash_func& get_hashfunc() const { return hf; }
    const hash_func& hash_function() const { return hf; }
    size_t bucket(const key_type& k) const { return hf(k) % table.size(); }
    size_t size() const { return num_data; }
    size_t num_buckets() const { return table.size(); }
    
    auto cbegin(size_t bucket_id) const {
        return table[bucket_id].cbegin();
    }
    auto cend(size_t bucket_id) const {
        return table[bucket_id].cend();
    }
    
    // useful methods
    data_type& operator[](const key_type& k) {
        auto& list = get_bucket_list(k);
        auto it = list.begin(), end = list.end();
        
        // loop over bucket elements
        for(; it != end; ++it){
            if( it->first == k ){
                return it->second;
            }
        }// end for
        
        // if did not find the element, add it
        list.emplace_front(hf(k), k);
        num_data++;
        
        // return the new data
        return list.front().second;
    }
    const data_type& operator[](const key_type& k) const {
        auto& list = get_bucket_list(k);
        auto it = list.cbegin(), end = list.cend();
        
        // loop over bucket elements
        for(; it != end; ++it){
            if( it->key == k ){
                return it->data;
            }
        }// end for
        
        return dummy;
        
    }
    const std::list<struct node>& get_bucket_list(const key_type& k) const {
        last_hash = hf(k);
        last_bucket = last_hash % table.size();
        return table[last_bucket];
    }
    
private:
    
    size_t last_hash, last_bucket;
    size_t num_data;
    hash_func hf;
    data_type dummy;
    std::vector<std::list<struct node>> table;
    
    // useful private methods
    std::list<struct node>& get_bucket_list(const key_type& k) {
        last_hash = hf(k);
        last_bucket = last_hash % table.size();
        return table[last_bucket];
    }
    
};


}

#endif /* cjh_hashtable_h */
