#ifndef ENERGY_MATRIX_H
#define ENERGY_MATRIX_H

#include "base_types.hh"
#include "sparse_tree.hh"
#include "SHAPE.hh"
#include "matrices.hh"
#include <string>
#include <vector>

#include "ViennaRNA/loops.hh"
#include "ViennaRNA/pair_mat.hh"
#include "ViennaRNA/params/io.hh"

class s_energy_matrix {
  public:

    s_energy_matrix(std::string seq, cand_pos_t length, SHAPEData *ShapeData, short *S, short *S1, vrna_param_t *params);
    // The constructor

    ~s_energy_matrix();
    // The destructor

    vrna_param_t *params_;

    short *S_;
    short *S1_;
    SHAPEData *ShapeData;
    // VM_sub should be NULL if you don't want suboptimals

    // void compute_energy (int i, int j);
    // compute the V(i,j) value

    free_energy_node *get_node(cand_pos_t i, cand_pos_t j) {
        cand_pos_t ij = index[i] + j - i;
        return &nodes[ij];
    }
    // return the node at (i,j)

    // May 15, 2007. Added "if (i>=j) return INF;"  below. It was miscalculating the backtracked structure.
    energy_t get_energy(cand_pos_t i, cand_pos_t j) {
        if (i >= j) return INF;
        cand_pos_t ij = index[i] + j - i;
        return nodes[ij].energy;
    }

    // return the value at V(i,j)

    char get_type(cand_pos_t i, cand_pos_t j) {
        cand_pos_t ij = index[i] + j - i;
        return nodes[ij].type;
    }
    // return the type at V(i,j)
    // Mateo 13 Sept 2023
    void compute_hotspot_energy(cand_pos_t i, cand_pos_t j, bool is_stack);
    energy_t HairpinE(const std::string &seq, const short *S, const short *S1, const vrna_param_t *params, cand_pos_t i, cand_pos_t j);
    energy_t compute_stack(cand_pos_t i, cand_pos_t j, const vrna_param_t *params);
    // better to have protected variable rather than private, it's necessary for Hfold
  protected:
    // private:

    std::string seq_;
    cand_pos_t n; // sequence length
    std::vector<cand_pos_t> index;
    // int *index;                // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
    std::vector<free_energy_node> nodes; // the free energy and type (i.e. base pair closing a hairpin loops, stacked pair etc), for each i and j
};

#endif
