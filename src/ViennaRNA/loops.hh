#ifndef LOOPS
#define LOOPS

#include "ViennaRNA/params/basic.hh"
#include "ViennaRNA/params/default.hh"
#include <cmath>
#include <cstring>
#include <algorithm>

/**                                External                                       */

inline int E_ExtLoop(int type, int si1, int sj1, vrna_param_t *P){
  int energy = 0;

  if (si1 >= 0 && sj1 >= 0)
    energy += P->mismatchExt[type][si1][sj1];
  else if (si1 >= 0)
    energy += P->dangle5[type][si1];
  else if (sj1 >= 0)
    energy += P->dangle3[type][sj1];

  if (type > 2)
    energy += P->TerminalAU;

  return energy;
}

inline double exp_E_ExtLoop(unsigned int type, int n5d, int n3d, vrna_exp_param_t *p){
  double energy = 1.0;

  if (n5d >= 0 && n3d >= 0)
    energy = p->expmismatchExt[type][n5d][n3d];
  else if (n5d >= 0)
    energy = p->expdangle5[type][n5d];
  else if (n3d >= 0)
    energy = p->expdangle3[type][n3d];

  if (type > 2)
    energy *= p->expTermAU;

  return (double) energy;
}


/**                                Hairpin                                       */






/**
 *  @brief Compute the Energy of a hairpin-loop
 *
 *  To evaluate the free energy of a hairpin-loop, several parameters have to be known.
 *  A general hairpin-loop has this structure:<BR>
 *  <PRE>
 *        a3 a4
 *      a2     a5
 *      a1     a6
 *        X - Y
 *        |   |
 *        5'  3'
 *  </PRE>
 *  where X-Y marks the closing pair [e.g. a <B>(G,C)</B> pair]. The length of this loop is 6 as there are
 *  six unpaired nucleotides (a1-a6) enclosed by (X,Y). The 5' mismatching nucleotide is
 *  a1 while the 3' mismatch is a6. The nucleotide sequence of this loop is &quot;a1.a2.a3.a4.a5.a6&quot; <BR>
 *  @note The parameter sequence should contain the sequence of the loop in capital letters of the nucleic acid
 *  alphabet if the loop size is below 7. This is useful for unusually stable tri-, tetra- and hexa-loops
 *  which are treated differently (based on experimental data) if they are tabulated.
 *  @see scale_parameters()
 *  @see vrna_param_t
 *  @warning Not (really) thread safe! A threadsafe implementation will replace this function in a future release!\n
 *  Energy evaluation may change due to updates in global variable "tetra_loop"
 *
 *  @param  size  The size of the loop (number of unpaired nucleotides)
 *  @param  type  The pair type of the base pair closing the hairpin
 *  @param  si1   The 5'-mismatching nucleotide
 *  @param  sj1   The 3'-mismatching nucleotide
 *  @param  string  The sequence of the loop (May be @p NULL, otherwise mst be at least @f$size + 2@f$ long)
 *  @param  P     The datastructure containing scaled energy parameters
 *  @return The Free energy of the Hairpin-loop in dcal/mol
 */
inline int E_Hairpin(int size, int type, int si1, int sj1,const char *string,vrna_param_t *P){
  int energy;

  if (size <= 30)
    energy = P->hairpin[size];
  else
    energy = P->hairpin[30] + (int)(P->lxc * log((size) / 30.));

  if (size < 3)
    return energy;            /* should only be the case when folding alignments */

  if ((string) && (P->model_details.special_hp)) {
    if (size == 4) {
      /* check for tetraloop bonus */
      char tl[7] = {
        0
      }, *ts;
      memcpy(tl, string, sizeof(char) * 6);
      tl[6] = '\0';
      if ((ts = strstr(P->Tetraloops, tl)))
        return P->Tetraloop_E[(ts - P->Tetraloops) / 7];
    } else if (size == 6) {
      char tl[9] = {
        0
      }, *ts;
      memcpy(tl, string, sizeof(char) * 8);
      tl[8] = '\0';
      if ((ts = strstr(P->Hexaloops, tl)))
        return energy = P->Hexaloop_E[(ts - P->Hexaloops) / 9];
    } else if (size == 3) {
      char tl[6] = {
        0
      }, *ts;
      memcpy(tl, string, sizeof(char) * 5);
      tl[5] = '\0';
      if ((ts = strstr(P->Triloops, tl)))
        return P->Triloop_E[(ts - P->Triloops) / 6];

      return energy + (type > 2 ? P->TerminalAU : 0);
    }
  }

  energy += P->mismatchH[type][si1][sj1];

  return energy;
}

/**
 *  @brief Compute Boltzmann weight @f$e^{-\Delta G/kT} @f$ of a hairpin loop
 *
 *  multiply by scale[u+2]
 *  @see get_scaled_pf_parameters()
 *  @see vrna_exp_param_t
 *  @see E_Hairpin()
 *  @warning Not (really) thread safe! A threadsafe implementation will replace this function in a future release!\n
 *  Energy evaluation may change due to updates in global variable "tetra_loop"
 *
 *  @param  u       The size of the loop (number of unpaired nucleotides)
 *  @param  type    The pair type of the base pair closing the hairpin
 *  @param  si1     The 5'-mismatching nucleotide
 *  @param  sj1     The 3'-mismatching nucleotide
 *  @param  string  The sequence of the loop (May be @p NULL, otherwise mst be at least @f$size + 2@f$ long)
 *  @param  P       The datastructure containing scaled Boltzmann weights of the energy parameters
 *  @return The Boltzmann weight of the Hairpin-loop
 */
inline double exp_E_Hairpin(int u, int type, short si1, short sj1, const char *string, vrna_exp_param_t *P){
  double q, kT;

  kT = P->kT;   /* kT in cal/mol  */

  if (u <= 30)
    q = P->exphairpin[u];
  else
    q = P->exphairpin[30] * exp(-(P->lxc * log(u / 30.)) * 10. / kT);

  if (u < 3)
    return (double)q;         /* should only be the case when folding alignments */

  if ((string) && (P->model_details.special_hp)) {
    if (u == 4) {
      char tl[7] = {
        0
      }, *ts;
      memcpy(tl, string, sizeof(char) * 6);
      tl[6] = '\0';
      if ((ts = strstr(P->Tetraloops, tl))) {
        if (type != 7)
          return (double)(P->exptetra[(ts - P->Tetraloops) / 7]);
        else
          q *= P->exptetra[(ts - P->Tetraloops) / 7];
      }
    } else if (u == 6) {
      char tl[9] = {
        0
      }, *ts;
      memcpy(tl, string, sizeof(char) * 8);
      tl[8] = '\0';
      if ((ts = strstr(P->Hexaloops, tl)))
        return (double)(P->exphex[(ts - P->Hexaloops) / 9]);
    } else if (u == 3) {
      char tl[6] = {
        0
      }, *ts;
      memcpy(tl, string, sizeof(char) * 5);
      tl[5] = '\0';
      if ((ts = strstr(P->Triloops, tl)))
        return (double)(P->exptri[(ts - P->Triloops) / 6]);

      if (type > 2)
        return (double)(q * P->expTermAU);
      else
        return (double)q;
    }
  }

  q *= P->expmismatchH[type][si1][sj1];

  return (double)q;
}


/**                                Internal                                       */


inline int E_IntLoop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1, vrna_param_t *P){
  /* compute energy of degree 2 loop (stack bulge or interior) */

   if (n1 < 0 || n2 < 0) {
        std::fprintf(
            stderr,
            "E_IntLoop: invalid lengths: n1=%d, n2=%d\n",
            n1,
            n2
        );
        std::abort();
    }

  int nl, ns, u, energy;

  energy = INF;

  if (n1 > n2) {
    nl  = n1;
    ns  = n2;
  } else {
    nl  = n2;
    ns  = n1;
  }

  if (nl == 0)
    return P->stack[type][type_2];  /* stack */

  if (ns == 0) {
    /* bulge */
    energy = (nl <= MAXLOOP) ? P->bulge[nl] :
             (P->bulge[30] + (int)(P->lxc * log(nl / 30.)));
    if (nl == 1) {
      energy += P->stack[type][type_2];
    } else {
      if (type > 2)
        energy += P->TerminalAU;

      if (type_2 > 2)
        energy += P->TerminalAU;
    }

    return energy;
  } else {
    /* interior loop */
    if (ns == 1) {
      if (nl == 1)                    /* 1x1 loop */
        return P->int11[type][type_2][si1][sj1];

      if (nl == 2) {
        /* 2x1 loop */
        if (n1 == 1)
          energy = P->int21[type][type_2][si1][sq1][sj1];
        else
          energy = P->int21[type_2][type][sq1][si1][sp1];

        return energy;
      } else {
        /* 1xn loop */
        energy =
          (nl + 1 <=
           MAXLOOP) ? (P->internal_loop[nl + 1]) : (P->internal_loop[30] +
                                                    (int)(P->lxc * log((nl + 1) / 30.)));
        energy  += std::min(MAX_NINIO, (nl - ns) * P->ninio[2]);
        energy  += P->mismatch1nI[type][si1][sj1] + P->mismatch1nI[type_2][sq1][sp1];
        return energy;
      }
    } else if (ns == 2) {
      if (nl == 2) {
        /* 2x2 loop */
        return P->int22[type][type_2][si1][sp1][sq1][sj1];
      } else if (nl == 3) {
        /* 2x3 loop */
        energy  = P->internal_loop[5] + P->ninio[2];
        energy  += P->mismatch23I[type][si1][sj1] + P->mismatch23I[type_2][sq1][sp1];
        return energy;
      }
    }

    {
      /* generic interior loop (no else here!)*/
      u       = nl + ns;
      energy  =
        (u <=
         MAXLOOP) ? (P->internal_loop[u]) : (P->internal_loop[30] + (int)(P->lxc * log((u) / 30.)));

      energy += std::min(MAX_NINIO, (nl - ns) * P->ninio[2]);

      energy += P->mismatchI[type][si1][sj1] + P->mismatchI[type_2][sq1][sp1];
    }
  }

  return energy;
}


inline double exp_E_IntLoop(int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1,vrna_exp_param_t *P){
  int     ul, us, no_close = 0;
  double  z           = 0.;
  int     noGUclosure = P->model_details.noGUclosure;

  if ((noGUclosure) && ((type2 == 3) || (type2 == 4) || (type == 3) || (type == 4)))
    no_close = 1;

  if (u1 > u2) {
    ul  = u1;
    us  = u2;
  } else {
    ul  = u2;
    us  = u1;
  }

  if (ul == 0) {
    /* stack */
    z = P->expstack[type][type2];
  } else if (!no_close) {
    if (us == 0) {
      /* bulge */
      z = P->expbulge[ul];
      if (ul == 1) {
        z *= P->expstack[type][type2];
      } else {
        if (type > 2)
          z *= P->expTermAU;

        if (type2 > 2)
          z *= P->expTermAU;
      }

      return (double)z;
    } else if (us == 1) {
      if (ul == 1)                     /* 1x1 loop */
        return (double)(P->expint11[type][type2][si1][sj1]);

      if (ul == 2) {
        /* 2x1 loop */
        if (u1 == 1)
          return (double)(P->expint21[type][type2][si1][sq1][sj1]);
        else
          return (double)(P->expint21[type2][type][sq1][si1][sp1]);
      } else {
        /* 1xn loop */
        z = P->expinternal[ul + us] * P->expmismatch1nI[type][si1][sj1] *
            P->expmismatch1nI[type2][sq1][sp1];
        return (double)(z * P->expninio[2][ul - us]);
      }
    } else if (us == 2) {
      if (ul == 2) {
        /* 2x2 loop */
        return (double)(P->expint22[type][type2][si1][sp1][sq1][sj1]);
      } else if (ul == 3) {
        /* 2x3 loop */
        z = P->expinternal[5] * P->expmismatch23I[type][si1][sj1] *
            P->expmismatch23I[type2][sq1][sp1];
        return (double)(z * P->expninio[2][1]);
      }
    }

    /* generic interior loop (no else here!)*/
    z = P->expinternal[ul + us] * P->expmismatchI[type][si1][sj1] *
        P->expmismatchI[type2][sq1][sp1];
    return (double)(z * P->expninio[2][ul - us]);
  }

  return (double)z;
}

/**                                Multiloop                                       */

/**
 *  @def E_MLstem(A,B,C,D)
 *  <H2>Compute the Energy contribution of a Multiloop stem</H2>
 *  This definition is a wrapper for the E_Stem() funtion.
 *  It is substituted by an E_Stem() funtion call with argument
 *  extLoop=0, so the energy contribution returned reflects a
 *  stem introduced in a multiloop.<BR>
 *  As for the parameters B (si1) and C (sj1) of the substituted
 *  E_Stem() function, you can inhibit to take 5'-, 3'-dangles
 *  or mismatch contributions to be taken into account by passing
 *  -1 to these parameters.
 *
 *  @see    E_Stem()
 *  @param  A The pair type of the stem-closing pair
 *  @param  B The 5'-mismatching nucleotide
 *  @param  C The 3'-mismatching nucleotide
 *  @param  D The datastructure containing scaled energy parameters
 *  @return   The energy contribution of the introduced multiloop stem
 */
inline int E_MLstem(int type,int si1,int sj1,vrna_param_t *P){
  int energy = 0;

  if (si1 >= 0 && sj1 >= 0)
    energy += P->mismatchM[type][si1][sj1];
  else if (si1 >= 0)
    energy += P->dangle5[type][si1];
  else if (sj1 >= 0)
    energy += P->dangle3[type][sj1];

  if (type > 2)
    energy += P->TerminalAU;

  energy += P->MLintern[type];

  return energy;
}

/**
 *  @def exp_E_MLstem(A,B,C,D)
 *  This is the partition function variant of @ref E_MLstem()
 *  @see E_MLstem()
 *  @return The Boltzmann weighted energy contribution of the introduced multiloop stem
 */
inline double exp_E_MLstem(int type, int si1, int sj1, vrna_exp_param_t *P){
  double energy = 1.0;

  if (si1 >= 0 && sj1 >= 0)
    energy = P->expmismatchM[type][si1][sj1];
  else if (si1 >= 0)
    energy = P->expdangle5[type][si1];
  else if (sj1 >= 0)
    energy = P->expdangle3[type][sj1];

  if (type > 2)
    energy *= P->expTermAU;

  energy *= P->expMLintern[type];
  return (double)energy;
}

#endif