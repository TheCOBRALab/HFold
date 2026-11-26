#ifndef HFOLD_VIENNA_HH
#define HFOLD_VIENNA_HH

// HFold's interface to ViennaRNA, which comes from ThermoRNA's ViennaRNA component.
// It keeps the names HFold's code has always used and maps them onto ThermoRNA.

#include <ViennaRNA/eval/exterior.hpp>
#include <ViennaRNA/eval/hairpin.hpp>
#include <ViennaRNA/eval/internal.hpp>
#include <ViennaRNA/eval/multibranch.hpp>
#include <ViennaRNA/model.hpp>
#include <ViennaRNA/params/basic.hpp>
#include <ViennaRNA/params/constants.hpp>
#include <ViennaRNA/params/io.hpp>
#include <ViennaRNA/sequences/alphabet.hpp>
#include <ViennaRNA/utils/utils.hpp>

#include <cstring>

using thermorna::viennarna::vrna_message_warning;
using thermorna::viennarna::vrna_param_s;
using thermorna::viennarna::vrna_param_t;
using thermorna::viennarna::vrna_params;
using thermorna::viennarna::vrna_params_load;
using thermorna::viennarna::vrna_params_load_DNA_Mathews2004;

// Pair types, as ViennaRNA's old pair_mat.h provided them: pair[a][b] is the type of the pair
// of base codes a and b (0 if they cannot pair), and rtype[t] is the type of the reversed pair.
inline int pair[MAXALPHA + 1][MAXALPHA + 1];
inline unsigned int rtype[8];

inline void make_pair_matrix() {
    thermorna::viennarna::vrna_md_t md;
    thermorna::viennarna::vrna_md_set_default(&md);
    std::memcpy(pair, md.pair, sizeof(pair));
    std::memcpy(rtype, md.rtype, sizeof(rtype));
}

// how == 0: base codes (S). how == 1: base codes with aliases, used for mismatches (S1).
// S[0] holds the length for how == 0, as in ViennaRNA.
inline short *encode_sequence(const char *sequence, short how) {
    thermorna::viennarna::vrna_md_t md;
    thermorna::viennarna::vrna_md_set_default(&md);
    return how ? thermorna::viennarna::vrna_seq_encode(sequence, &md)
               : thermorna::viennarna::vrna_seq_encode_simple(sequence, &md);
}

// Loop energies under their old ViennaRNA names
inline int E_Hairpin(int size, int type, int si1, int sj1, const char *string, vrna_param_t *P) {
    return thermorna::viennarna::vrna_E_hairpin(size, type, si1, sj1, string, P);
}

inline int E_IntLoop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1, vrna_param_t *P) {
    return thermorna::viennarna::vrna_E_internal(n1, n2, type, type_2, si1, sj1, sp1, sq1, P);
}

inline int E_ExtLoop(int type, int si1, int sj1, vrna_param_t *P) {
    return thermorna::viennarna::vrna_E_exterior_stem(type, si1, sj1, P);
}

inline int E_MLstem(int type, int si1, int sj1, vrna_param_t *P) {
    return thermorna::viennarna::vrna_E_multibranch_stem(type, si1, sj1, P);
}

#endif
