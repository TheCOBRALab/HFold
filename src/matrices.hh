#ifndef MATRICES_H
#define MATRICES_H

#include "base_types.hh"
#include <vector>
#include <iostream>
#include <limits>
#include <cassert>
#define INF 10000000 /* (INT_MAX/10) */

typedef std::vector<size_t> index_offset_t;

class TriangleMatrix {
public:
    TriangleMatrix() {}

    void init(cand_pos_t n, std::vector<cand_pos_t> &index, energy_t return_val=INF){
        n_ = n;
        index_ = &index;
        return_val_ = return_val;

        size_t tl=total_length(n);
        m_.clear();
        m_.resize(tl,INF+1);
    }

    cand_pos_t ij(cand_pos_t i, cand_pos_t j) const {
        return (*index_)[i]+j-i;
    }

    //! unchecked get
    energy_t get_uc(cand_pos_t i, cand_pos_t j) const {
        assert (i<=j);
        return m_[ij(i,j)];
    }

    energy_t get (cand_pos_t i, cand_pos_t j) const {
        if (i>j) return return_val_;
        return get_uc(i,j);
    }


    energy_t& operator [] (cand_pos_t ij) {
        return m_[ij];
    }

    energy_t& set (cand_pos_t i, cand_pos_t j) {
        return m_[ij(i,j)];
    }

    void print() {
        for (cand_pos_t i=0; i<n_; i++) {
            std::cout << i << ": ";
            for (cand_pos_t j=i; j<n_; j++) {
                std::cout << get(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }

    static cand_pos_t total_length(cand_pos_t n) {
        return (n *(n+1))/2;
    }

    static void new_index(std::vector<cand_pos_t> &index, cand_pos_t n) {
        index.resize(n);
        index[1] = 0;
        for (cand_pos_t i=2; i < n; i++) {
            index[i] = index[i-1]+n-i+1;
        }
    }

private:
    cand_pos_t n_;
    energy_t return_val_;
    std::vector<cand_pos_t> *index_;
    std::vector<energy_t> m_;
};
#endif