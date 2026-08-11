#include "ViennaRNA/params/default.hh"
#include "ViennaRNA/params/basic.hh"
#include "ViennaRNA/params/io.hh"
#include "ViennaRNA/utils.hh"
#include "ViennaRNA/model.hh"
#include <cstring>
#include <cmath>
#include <algorithm>


/**
 *** \file ViennaRNA/params/basic.c
 *** <P>
 *** This file provides functions that return temperature scaled energy parameters and
 *** Boltzmann weights packed in datastructures
 *** </P>
 ***/

/*------------------------------------------------------------------------*/
#define SCALE 10

/**
 *** dangling ends should never be destabilizing, i.e. expdangle>=1<BR>
 *** specific heat needs smooth function (2nd derivative)<BR>
 *** we use a*(sin(x+b)+1)^2, with a=2/(3*sqrt(3)), b=Pi/6-sqrt(3)/2,
 *** in the interval b<x<sqrt(3)/2
 */
#define CLIP_NEGATIVE(X) ((X) < 0 ? 0 : (X))
#define SMOOTH(X) ((!pf_smooth)               ?   CLIP_NEGATIVE(X) : \
                   ((X) / SCALE < -1.2283697) ?                 0  : \
                   ((X) / SCALE > 0.8660254) ?                (X) : \
                   SCALE *0.38490018   \
                   * (sin((X) / SCALE - 0.34242663) + 1) \
                   * (sin((X) / SCALE - 0.34242663) + 1) \
                   )

/* #define SMOOTH(X) ((X)<0 ? 0 : (X)) */

/*
 * If the global use_mfelike_energies flag is set, truncate doubles to int
 * values and cast back to double. This makes the energy parameters of the
 * partition (folding get_scaled_exp_params()) compatible with the mfe folding
 * parameters (get_scaled_exp_params()), e.g. for explicit partition function
 * computations.
 */
#define TRUNC_MAYBE(X) ((!pf_smooth) ? (double)((int)(X)) : (X))


/* Rescale Free energy contribution according to deviation of temperature from measurement conditions */
#define RESCALE_dG(dG, dH, dT)   ((dH) - ((dH) - (dG)) * dT)

/*
 * Rescale Free energy contribution according to deviation of temperature from measurement conditions
 * and convert it to Boltzmann Factor for specific kT
 */
#define RESCALE_BF(dG, dH, dT, kT)          ( \
    exp( \
      -TRUNC_MAYBE((double)RESCALE_dG((dG), (dH), (dT))) \
      * 10. \
      / kT \
      ) \
    )

#define RESCALE_BF_SMOOTH(dG, dH, dT, kT)   ( \
    exp(  \
      SMOOTH( \
        -TRUNC_MAYBE((double)RESCALE_dG((dG), (dH), (dT))) \
        ) \
      * 10. \
      / kT \
      ) \
    )

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */
static int id = -1;

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

vrna_param_t* get_scaled_params(vrna_md_t *md);


vrna_exp_param_t* get_scaled_exp_params(vrna_md_t *md, double pfs);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
vrna_param_t* vrna_params(vrna_md_t *md){
  if (md) {
    return get_scaled_params(md);
  } else {
    vrna_md_t md;
    vrna_md_set_default(&md);
    return get_scaled_params(&md);
  }
}

vrna_exp_param_t* vrna_exp_params(vrna_md_t *md){
  if (md) {
    return get_scaled_exp_params(md, -1.);
  } else {
    vrna_md_t md;
    vrna_md_set_default(&md);
    return get_scaled_exp_params(&md, -1.);
  }
}

/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
vrna_param_t* get_scaled_params(vrna_md_t *md){
  double        tempf;
  vrna_param_t  *params;

  params = (vrna_param_t *)vrna_alloc(sizeof(vrna_param_t));

  memset(params->param_file, '\0', 256);
  if (last_parameter_file() != NULL)
    strncpy(params->param_file, last_parameter_file(), 255);

  params->model_details = *md;  /* copy over the model details */
  params->temperature   = md->temperature;
  tempf                 = ((params->temperature + K0) / Tmeasure);

  params->ninio[2]              = RESCALE_dG(ninio37, niniodH, tempf);
  params->lxc                   = lxc37 * tempf;
  params->TripleC               = RESCALE_dG(TripleC37, TripleCdH, tempf);
  params->MultipleCA            = RESCALE_dG(MultipleCA37, MultipleCAdH, tempf);
  params->MultipleCB            = RESCALE_dG(MultipleCB37, MultipleCBdH, tempf);
  params->TerminalAU            = RESCALE_dG(TerminalAU37, TerminalAUdH, tempf);
  params->DuplexInit            = RESCALE_dG(DuplexInit37, DuplexInitdH, tempf);
  params->MLbase                = RESCALE_dG(ML_BASE37, ML_BASEdH, tempf);
  // Pseudoknot -- added by Mateo 2/14/23
  params->PS_penalty            = -138;
  params->PSM_penalty           = 1007; 
  params->PSP_penalty           = 1500; 
  params->PB_penalty            = 246; 
  params->PUP_penalty           = 6; 
  params->PPS_penalty           = 96; 
  params->e_stP_penalty         = .89; 
  params->e_intP_penalty        = .74;
  params->ap_penalty            = 341;
  params->bp_penalty            = 56;
  params->cp_penalty            = 12;
  params->a_penalty            = 339; // Fix later
  params->b_penalty            = 3;
  params->c_penalty            = 2;
  //
  params->MLclosing             = RESCALE_dG(ML_closing37, ML_closingdH, tempf);
  params->gquadLayerMismatch    = RESCALE_dG(GQuadLayerMismatch37, GQuadLayerMismatchH, tempf);
  params->gquadLayerMismatchMax = GQuadLayerMismatchMax;

  for (unsigned int i = VRNA_GQUAD_MIN_STACK_SIZE; i <= VRNA_GQUAD_MAX_STACK_SIZE; i++)
    for (unsigned int j = 3 * VRNA_GQUAD_MIN_LINKER_LENGTH; j <= 3 * VRNA_GQUAD_MAX_LINKER_LENGTH; j++) {
      double  GQuadAlpha_T  = RESCALE_dG(GQuadAlpha37, GQuadAlphadH, tempf);
      double  GQuadBeta_T   = RESCALE_dG(GQuadBeta37, GQuadBetadH, tempf);
      params->gquad[i][j] = (int)GQuadAlpha_T * (i - 1) + (int)(((double)GQuadBeta_T) * log(j - 2));
    }

  for (unsigned int i = 0; i < 31; i++)
    params->hairpin[i] = RESCALE_dG(hairpin37[i], hairpindH[i], tempf);

  for (unsigned int i = 0; i <= MAX_INTLOOP; i++) {
    params->bulge[i]          = RESCALE_dG(bulge37[i], bulgedH[i], tempf);
    params->internal_loop[i]  = RESCALE_dG(internal_loop37[i], internal_loopdH[i], tempf);
  }

  for (unsigned int i = MAX_INTLOOP + 1; i <= MAXLOOP; i++) {
    params->bulge[i] = params->bulge[30] +
                       (int)(params->lxc * log((double)(i) / 30.));
    params->internal_loop[i] = params->internal_loop[30] +
                               (int)(params->lxc * log((double)(i) / 30.));
  }

  for (unsigned int i = 0; (i * 7) < strlen(Tetraloops); i++)
    params->Tetraloop_E[i] = RESCALE_dG(Tetraloop37[i], TetraloopdH[i], tempf);

  for (unsigned int i = 0; (i * 5) < strlen(Triloops); i++)
    params->Triloop_E[i] = RESCALE_dG(Triloop37[i], TriloopdH[i], tempf);

  for (unsigned int i = 0; (i * 9) < strlen(Hexaloops); i++)
    params->Hexaloop_E[i] = RESCALE_dG(Hexaloop37[i], HexaloopdH[i], tempf);

  for (unsigned int i = 0; i <= NBPAIRS; i++)
    params->MLintern[i] = RESCALE_dG(ML_intern37, ML_interndH, tempf);

  /* stacks    G(T) = H - [H - G(T0)]*T/T0 */
  for (unsigned int i = 0; i <= NBPAIRS; i++)
    for (unsigned int j = 0; j <= NBPAIRS; j++)
      params->stack[i][j] = RESCALE_dG(stack37[i][j],
                                       stackdH[i][j],
                                       tempf);

  /* mismatches */
  for (unsigned int i = 0; i <= NBPAIRS; i++)
    for (unsigned int j = 0; j < 5; j++)
      for (unsigned int k = 0; k < 5; k++) {
        int mm;
        params->mismatchI[i][j][k] = RESCALE_dG(mismatchI37[i][j][k],
                                                mismatchIdH[i][j][k],
                                                tempf);
        params->mismatchH[i][j][k] = RESCALE_dG(mismatchH37[i][j][k],
                                                mismatchHdH[i][j][k],
                                                tempf);
        params->mismatch1nI[i][j][k] = RESCALE_dG(mismatch1nI37[i][j][k],
                                                  mismatch1nIdH[i][j][k],
                                                  tempf);
        params->mismatch23I[i][j][k] = RESCALE_dG(mismatch23I37[i][j][k],
                                                  mismatch23IdH[i][j][k],
                                                  tempf);
        if (md->dangles) {
          mm = RESCALE_dG(mismatchM37[i][j][k],
                          mismatchMdH[i][j][k],
                          tempf);
          params->mismatchM[i][j][k]  = (mm > 0) ? 0 : mm;
          mm                          = RESCALE_dG(mismatchExt37[i][j][k],
                                                   mismatchExtdH[i][j][k],
                                                   tempf);
          params->mismatchExt[i][j][k] = (mm > 0) ? 0 : mm;
        } else {
          params->mismatchM[i][j][k] = params->mismatchExt[i][j][k] = 0;
        }
      }

  /* dangles */
  for (unsigned int i = 0; i <= NBPAIRS; i++)
    for (unsigned int j = 0; j < 5; j++) {
      int dd;
      dd = RESCALE_dG(dangle5_37[i][j],
                      dangle5_dH[i][j],
                      tempf);
      params->dangle5[i][j] = (dd > 0) ? 0 : dd;  /* must be <= 0 */
      dd                    = RESCALE_dG(dangle3_37[i][j],
                                         dangle3_dH[i][j],
                                         tempf);
      params->dangle3[i][j] = (dd > 0) ? 0 : dd;  /* must be <= 0 */
    }

  /* interior 1x1 loops */
  for (unsigned int i = 0; i <= NBPAIRS; i++)
    for (unsigned int j = 0; j <= NBPAIRS; j++)
      for (unsigned int k = 0; k < 5; k++)
        for (unsigned int l = 0; l < 5; l++)
          params->int11[i][j][k][l] = RESCALE_dG(int11_37[i][j][k][l],
                                                 int11_dH[i][j][k][l],
                                                 tempf);

  /* interior 2x1 loops */
  for (unsigned int i = 0; i <= NBPAIRS; i++)
    for (unsigned int j = 0; j <= NBPAIRS; j++)
      for (unsigned int k = 0; k < 5; k++)
        for (unsigned int l = 0; l < 5; l++) {
          int m;
          for (m = 0; m < 5; m++)
            params->int21[i][j][k][l][m] = RESCALE_dG(int21_37[i][j][k][l][m],
                                                      int21_dH[i][j][k][l][m],
                                                      tempf);
        }

  /* interior 2x2 loops */
  for (unsigned int i = 0; i <= NBPAIRS; i++)
    for (unsigned int j = 0; j <= NBPAIRS; j++)
      for (unsigned int k = 0; k < 5; k++)
        for (unsigned int l = 0; l < 5; l++) {
          for (unsigned int m = 0; m < 5; m++)
            for (unsigned int n = 0; n < 5; n++)
              params->int22[i][j][k][l][m][n] = RESCALE_dG(int22_37[i][j][k][l][m][n],
                                                           int22_dH[i][j][k][l][m][n],
                                                           tempf);
        }

  strncpy(params->Tetraloops, Tetraloops, 281);
  strncpy(params->Triloops, Triloops, 241);
  strncpy(params->Hexaloops, Hexaloops, 361);

  params->id = ++id;
  return params;
}


vrna_exp_param_t* get_scaled_exp_params(vrna_md_t *md,double pfs){
  int               pf_smooth;
  double            kT, TT;
  double            GT;
  vrna_exp_param_t  *pf;

  pf = (vrna_exp_param_t *)vrna_alloc(sizeof(vrna_exp_param_t));

  memset(pf->param_file, '\0', 256);
  if (last_parameter_file() != NULL)
    strncpy(pf->param_file, last_parameter_file(), 255);

  pf->model_details = *md;
  pf->temperature   = md->temperature;
  pf->alpha         = md->betaScale;
  pf->kT            = kT = md->betaScale * (md->temperature + K0) * GASCONST; /* kT in cal/mol  */
  pf->pf_scale      = pfs;
  pf_smooth         = md->pf_smooth;
  TT                = (md->temperature + K0) / (Tmeasure);

  pf->lxc                   = lxc37 * TT;
  pf->expDuplexInit         = RESCALE_BF(DuplexInit37, DuplexInitdH, TT, kT);
  pf->expTermAU             = RESCALE_BF(TerminalAU37, TerminalAUdH, TT, kT);
  pf->expMLbase             = RESCALE_BF(ML_BASE37, ML_BASEdH, TT, kT);
  pf->expMLclosing          = RESCALE_BF(ML_closing37, ML_closingdH, TT, kT);
  pf->expgquadLayerMismatch = RESCALE_BF(GQuadLayerMismatch37, GQuadLayerMismatchH, TT, kT);
  pf->gquadLayerMismatchMax = GQuadLayerMismatchMax;

  for (unsigned int i = VRNA_GQUAD_MIN_STACK_SIZE; i <= VRNA_GQUAD_MAX_STACK_SIZE; i++)
    for (unsigned int j = 3 * VRNA_GQUAD_MIN_LINKER_LENGTH; j <= 3 * VRNA_GQUAD_MAX_LINKER_LENGTH; j++) {
      double  GQuadAlpha_T  = RESCALE_dG(GQuadAlpha37, GQuadAlphadH, TT);
      double  GQuadBeta_T   = RESCALE_dG(GQuadBeta37, GQuadBetadH, TT);
      GT = ((double)GQuadAlpha_T) * ((double)(i - 1)) + ((double)GQuadBeta_T) *
           log(((double)j) - 2.);
      pf->expgquad[i][j] = exp(-TRUNC_MAYBE(GT) * 10. / kT);
    }

  /* loop energies: hairpins, bulges, interior, mulit-loops */
  for (unsigned int i = 0; i < 31; i++)
    pf->exphairpin[i] = RESCALE_BF(hairpin37[i], hairpindH[i], TT, kT);

  for (unsigned int i = 0; i <= MAX_INTLOOP; i++) {
    pf->expbulge[i]     = RESCALE_BF(bulge37[i], bulgedH[i], TT, kT);
    pf->expinternal[i]  = RESCALE_BF(internal_loop37[i], internal_loopdH[i], TT, kT);
  }

  /* special case of size 2 interior loops (single mismatch) */
  if (james_rule)
    pf->expinternal[2] = exp(-80 * 10. / kT);

  GT = RESCALE_dG(bulge37[30],
                  bulgedH[30],
                  TT);
  for (unsigned int i = 31; i <= MAXLOOP; i++)
    pf->expbulge[i] = exp(-TRUNC_MAYBE(GT + (pf->lxc * log(i / 30.))) * 10. / kT);

  GT = RESCALE_dG(internal_loop37[30],
                  internal_loopdH[30],
                  TT);
  for (unsigned int i = 31; i <= MAXLOOP; i++)
    pf->expinternal[i] = exp(-TRUNC_MAYBE(GT + (pf->lxc * log(i / 30.))) * 10. / kT);

  GT = RESCALE_dG(ninio37, niniodH, TT);
  for (unsigned int j = 0; j <= MAXLOOP; j++)
    pf->expninio[2][j] = exp(-std::min(MAX_NINIO, (int)(j * TRUNC_MAYBE(GT))) * 10. / kT);

  for (unsigned int i = 0; (i * 7) < strlen(Tetraloops); i++)
    pf->exptetra[i] = RESCALE_BF(Tetraloop37[i], TetraloopdH[i], TT, kT);

  for (unsigned int i = 0; (i * 5) < strlen(Triloops); i++)
    pf->exptri[i] = RESCALE_BF(Triloop37[i], TriloopdH[i], TT, kT);

  for (unsigned int i = 0; (i * 9) < strlen(Hexaloops); i++)
    pf->exphex[i] = RESCALE_BF(Hexaloop37[i], HexaloopdH[i], TT, kT);

  for (unsigned int i = 0; i <= NBPAIRS; i++)
    pf->expMLintern[i] = RESCALE_BF(ML_intern37, ML_interndH, TT, kT);

  /* if dangles==0 just set their energy to 0,
   * don't let dangle energies become > 0 (at large temps),
   * but make sure go smoothly to 0                        */
  for (unsigned int i = 0; i <= NBPAIRS; i++)
    for (unsigned int j = 0; j <= 4; j++) {
      if (md->dangles) {
        pf->expdangle5[i][j]  = RESCALE_BF_SMOOTH(dangle5_37[i][j], dangle5_dH[i][j], TT, kT);
        pf->expdangle3[i][j]  = RESCALE_BF_SMOOTH(dangle3_37[i][j], dangle3_dH[i][j], TT, kT);
      } else {
        pf->expdangle3[i][j] = pf->expdangle5[i][j] = 1;
      }
    }

  /* stacking energies */
  for (unsigned int i = 0; i <= NBPAIRS; i++)
    for (unsigned int j = 0; j <= NBPAIRS; j++)
      pf->expstack[i][j] = RESCALE_BF(stack37[i][j], stackdH[i][j], TT, kT);

  /* mismatch energies */
  for (unsigned int i = 0; i <= NBPAIRS; i++)
    for (unsigned int j = 0; j < 5; j++)
      for (unsigned int k = 0; k < 5; k++) {
        pf->expmismatchI[i][j][k] = RESCALE_BF(mismatchI37[i][j][k],
                                               mismatchIdH[i][j][k],
                                               TT,
                                               kT);
        pf->expmismatch1nI[i][j][k] = RESCALE_BF(mismatch1nI37[i][j][k],
                                                 mismatch1nIdH[i][j][k],
                                                 TT,
                                                 kT);
        pf->expmismatchH[i][j][k] = RESCALE_BF(mismatchH37[i][j][k],
                                               mismatchHdH[i][j][k],
                                               TT,
                                               kT);
        pf->expmismatch23I[i][j][k] = RESCALE_BF(mismatch23I37[i][j][k],
                                                 mismatch23IdH[i][j][k],
                                                 TT,
                                                 kT);

        if (md->dangles) {
          pf->expmismatchM[i][j][k] = RESCALE_BF_SMOOTH(mismatchM37[i][j][k],
                                                        mismatchMdH[i][j][k],
                                                        TT,
                                                        kT);
          pf->expmismatchExt[i][j][k] = RESCALE_BF_SMOOTH(mismatchExt37[i][j][k],
                                                          mismatchExtdH[i][j][k],
                                                          TT,
                                                          kT);
        } else {
          pf->expmismatchM[i][j][k] = pf->expmismatchExt[i][j][k] = 1.;
        }
      }

  /* interior lops of length 2 */
  for (unsigned int i = 0; i <= NBPAIRS; i++)
    for (unsigned int j = 0; j <= NBPAIRS; j++)
      for (unsigned int k = 0; k < 5; k++)
        for (unsigned int l = 0; l < 5; l++) {
          pf->expint11[i][j][k][l] = RESCALE_BF(int11_37[i][j][k][l],
                                                int11_dH[i][j][k][l],
                                                TT,
                                                kT);
        }

  /* interior 2x1 loops */
  for (unsigned int i = 0; i <= NBPAIRS; i++)
    for (unsigned int j = 0; j <= NBPAIRS; j++)
      for (unsigned int k = 0; k < 5; k++)
        for (unsigned int l = 0; l < 5; l++) {
          int m;
          for (m = 0; m < 5; m++) {
            pf->expint21[i][j][k][l][m] = RESCALE_BF(int21_37[i][j][k][l][m],
                                                     int21_dH[i][j][k][l][m],
                                                     TT,
                                                     kT);
          }
        }

  /* interior 2x2 loops */
  for (unsigned int i = 0; i <= NBPAIRS; i++)
    for (unsigned int j = 0; j <= NBPAIRS; j++)
      for (unsigned int k = 0; k < 5; k++)
        for (unsigned int l = 0; l < 5; l++) {
          int m, n;
          for (m = 0; m < 5; m++)
            for (n = 0; n < 5; n++) {
              pf->expint22[i][j][k][l][m][n] = RESCALE_BF(int22_37[i][j][k][l][m][n],
                                                          int22_dH[i][j][k][l][m][n],
                                                          TT,
                                                          kT);
            }
        }

  memcpy(pf->Tetraloops, Tetraloops, 281);
  pf->Tetraloops[280] = '\0';

  memcpy(pf->Triloops, Triloops, 241);
  pf->Triloops[240] = '\0';

  memcpy(pf->Hexaloops, Hexaloops, 361);
  pf->Hexaloops[360] = '\0';

  return pf;
}