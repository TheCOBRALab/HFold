/***************************************************************************
                          s_energy_matrix.cpp  -  description
                             -------------------
    begin                : Fri Apr 12 2002
    copyright            : (C) 2002 by Mirela Andronescu
    email                : andrones@cs.ubc.ca
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

// This is the V matrix

#include "constants.hh"
#include "h_globals.hh"
#include "h_struct.hh"

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

#include "s_energy_matrix.hh"

s_energy_matrix::s_energy_matrix(std::string seq, cand_pos_t length,  SHAPEData *ShapeData, short *S, short *S1, vrna_param_t *params)
// The constructor
{
    params_ = params;
    make_pair_matrix();
    S_ = S;
    S1_ = S1;

    n = length;
    seq_ = seq;
    this->ShapeData = ShapeData;

    // an vector with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
    index.resize(n + 1);
    cand_pos_t total_length = ((n + 1) * (n + 2)) / 2;
    TriangleMatrix::new_index(index,n+1);
    // this array holds V(i,j), and what (i,j) encloses: hairpin loop, stack pair, internal loop or multi-loop
    nodes.resize(total_length);
}

s_energy_matrix::~s_energy_matrix()
// The destructor
{}

/**
 * @brief This code returns the hairpin energy for a given base pair.
 * @param i The left index in the base pair
 * @param j The right index in the base pair
 */
energy_t s_energy_matrix::HairpinE(const std::string &seq, const short *S, const short *S1, const vrna_param_t *params, cand_pos_t i, cand_pos_t j) {

    const int ptype_closing = pair[S[i]][S[j]];

    if (ptype_closing == 0) return INF;

    return E_Hairpin(j - i - 1, ptype_closing, S1[i + 1], S1[j - 1], &seq.c_str()[i - 1], const_cast<vrna_param_t *>(params));
}

energy_t s_energy_matrix::compute_stack(cand_pos_t i, cand_pos_t j, const vrna_param_t *params) {

    const int ptype_closing = pair[S_[i]][S_[j]];
    cand_pos_t k = i + 1;
    cand_pos_t l = j - 1;
    return E_IntLoop(k - i - 1, j - l - 1, ptype_closing, rtype[pair[S_[k]][S_[l]]], S1_[i + 1], S1_[j - 1], S1_[k - 1], S1_[l + 1],
                     const_cast<vrna_param_t *>(params))
           + get_energy(k, l) + ShapeData->get_calculated(i) + ShapeData->get_calculated(j);
}

// Mateo 13 Sept 2023
void s_energy_matrix::compute_hotspot_energy(cand_pos_t i, cand_pos_t j, bool is_stack) {
    // printf("in compute_hotspot_energy i:%d j:%d\n",i,j);
    energy_t energy = 0;
    if (is_stack) {
        energy = compute_stack(i, j, params_);
        energy+=ShapeData->get_calculated(i);
        energy+=ShapeData->get_calculated(j);
        // printf("stack: %d\n",energy);
    } else {
        energy = 0; // HairpinE(seq_,S_,S1_,params_,i,j);
    }

    cand_pos_t ij = index[i] + j - i;
    nodes[ij].energy = energy;
    return;
}
