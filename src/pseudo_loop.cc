#include "pseudo_loop.hh"
#include "h_globals.hh"
#include <algorithm>
#include <iostream>
#include <string>

pseudo_loop::pseudo_loop(std::string seq, std::string res, sparse_tree &tree, SHAPEData &ShapeData, bool pk_free, bool pk_only, int dangle) : params_(vrna_params(NULL)) {
    this->seq_ = seq;
    this->res = res;
    this->tree = &tree;
    this->ShapeData = &ShapeData;
    params_->model_details.dangles = dangle;
    make_pair_matrix();
    S_ = encode_sequence(seq.c_str(), 0);
    S1_ = encode_sequence(seq.c_str(), 1);
    this->pk_free = pk_free;
    this->pk_only = pk_only;
    allocate_space();
}

void pseudo_loop::allocate_space() {
    n = seq_.length();

    structure = std::string(n+1, '.');
    index.resize(n + 1);
    TriangleMatrix::new_index(index,n+1);
    cand_pos_t total_length = ((n + 1) * (n + 2)) / 2;
    W.resize(n+1,0);
    WM.init(n+1,index);
    WMv.init(n+1,index);
    WMp.init(n+1,index);
    V.resize(total_length);

    WI.init(n+1,index,0);
    WIP.init(n+1,index);
    VP.init(n+1,index);
    VPL.init(n+1,index);
    VPR.init(n+1,index);
    WMB.init(n+1,index);
    WMBP.init(n+1,index);
    WMBW.init(n+1,index);
    BE.init(n+1,index,0);
}

pseudo_loop::~pseudo_loop() {
    free(S_);
    free(S1_);
    free(params_);
}

double pseudo_loop::hfold() {

    for (int i = n; i >= 1; --i) {
        for (int j = i; j <= n; ++j) // for (i=0; i<=j; i++)
        {
            const bool evaluate = tree->weakly_closed(i, j);
            const pair_type ptype_closing = pair[S_[i]][S_[j]];
            const bool restricted = tree->tree[i].pair == -1 || tree->tree[j].pair == -1;
            const bool paired = (tree->tree[i].pair == j && tree->tree[j].pair == i);

            const bool pkonly = (!pk_only || paired);

            if (ptype_closing > 0 && evaluate && !restricted && pkonly) compute_energy_restricted(i, j);

            if (!pk_free) compute_energies_PK(i, j);

            compute_WMv_WMp(i, j);
            compute_energy_WM_restricted(i, j);
        }
    }
    for (cand_pos_t j = TURN + 1; j <= n; j++) {
        energy_t m1 = INF;
        energy_t m2 = INF;
        energy_t m3 = INF;
        if (tree->tree[j].pair < 0) m1 = W[j - 1];

        for (cand_pos_t k = 1; k <= j - TURN - 1; ++k) {
            if (tree->weakly_closed(1, k - 1)) {
                energy_t acc = (k > 1) ? W[k - 1] : 0;
                m2 = std::min(m2, acc + E_ext_Stem(get_energy(k, j), get_energy(k + 1, j), get_energy(k, j - 1), get_energy(k + 1, j - 1), k, j));
                if (k == 1 || tree->weakly_closed(k, j)) m3 = std::min(m3, acc + WMB.get(k, j) + PS_penalty);
            }
        }
        W[j] = std::min({m1, m2, m3});
    }
    double energy = W[n] / 100.0;
    backtrack();
    this->structure = structure.substr(1, n);
    return energy;
}

void pseudo_loop::compute_WMv_WMp(cand_pos_t i, cand_pos_t j) {
    if (j - i + 1 < 4) return;
    energy_t m1 = E_MLStem(get_energy(i,j), get_energy(i+1,j), get_energy(i,j-1), get_energy(i+1,j-1), i, j);
    energy_t m2 = WMB.get(i,j) + PSM_penalty + b_penalty;
    if (tree->tree[j].pair <= -1) {
        m1 = std::min(m1,WMv.get(i,j-1) + params_->MLbase);
        m2 = std::min(m2,WMp.get(i,j-1) + params_->MLbase);
    }
    WMv.set(i,j) = m1;
    WMp.set(i,j) = m2;
}

void pseudo_loop::compute_energy_WM_restricted(cand_pos_t i, cand_pos_t j)
// compute de MFE of a partial multi-loop closed at (i,j), the restricted case
{
    if (j - i + 1 < 4) return;
    energy_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;

    for (cand_pos_t k = j - TURN - 1; k >= i; --k) {
        energy_t wm_kj = E_MLStem(get_energy(k,j), get_energy(k+1,j), get_energy(k,j-1), get_energy(k+1,j-1), k, j);
        energy_t wmb_kj = WMB.get(k,j) + PSM_penalty + b_penalty;
        bool can_pair = tree->up[k - 1] >= (k - i);
        if (can_pair) m1 = std::min(m1, static_cast<energy_t>((k - i) * params_->MLbase) + wm_kj);
        if (can_pair) m2 = std::min(m2, static_cast<energy_t>((k - i) * params_->MLbase) + wmb_kj);
        m3 = std::min(m3, WM.get(i,k-1) + wm_kj);
        m4 = std::min(m4, WM.get(i,k-1) + wmb_kj);
    }
    if (tree->tree[j].pair <= -1) m5 = WM.get(i,j-1) + params_->MLbase;
    WM.set(i,j) = std::min({m1, m2, m3, m4, m5});
}

energy_t pseudo_loop::compute_energy_VM_restricted(cand_pos_t i, cand_pos_t j)
// compute the MFE of a multi-loop closed at (i,j), the restricted case
{
    energy_t min = INF;
    for (cand_pos_t k = i + 1; k <= j - TURN - 1; ++k) {
        energy_t WM2ij = WM.get(i+1,k-1) + WMv.get(k,j-1);
        WM2ij = std::min(WM2ij, WM.get(i + 1, k - 1) + WMp.get(k, j - 1));
        if (tree->up[k - 1] >= (k - (i + 1))) WM2ij = std::min(WM2ij, static_cast<energy_t>((k - i - 1) * params_->MLbase) + WMp.get(k, j - 1));

        energy_t WM2ip1j = WM.get(i + 2, k - 1) + WMv.get(k, j - 1);
        WM2ip1j = std::min(WM2ip1j, WM.get(i + 2, k - 1) + WMp.get(k - 1, j - 1));
        if (tree->up[k - 1] >= (k - (i + 1)))
            WM2ip1j = std::min(WM2ip1j, static_cast<energy_t>((k - (i + 1) - 1) * params_->MLbase) + WMp.get(k, j - 1));

        energy_t WM2ijm1 = WM.get(i + 1, k - 1) + WMv.get(k, j - 2);
        WM2ijm1 = std::min(WM2ijm1, WM.get(i + 1, k - 1) + WMp.get(k, j - 2));
        if (tree->up[k - 1] >= (k - (i + 2)))
            WM2ijm1 = std::min(WM2ijm1, static_cast<energy_t>((k - i - 1) * params_->MLbase) + WMp.get(k, j - 2));

        energy_t WM2ip1jm1 = WM.get(i + 2, k - 1) + WMv.get(k, j - 2);
        WM2ip1jm1 = std::min(WM2ip1jm1, WM.get(i + 2, k - 1) + WMp.get(k, j - 2));
        if (tree->up[k - 2] >= (k - (i + 2)))
            WM2ip1jm1 = std::min(WM2ip1jm1, static_cast<energy_t>((k - (i + 1) - 1) * params_->MLbase) + WMp.get(k, j - 2));
        min = std::min(min, E_MbLoop(WM2ij, WM2ip1j, WM2ijm1, WM2ip1jm1, i, j));
    }
    return min;
}

/**
 * @brief restricted version
 */
energy_t pseudo_loop::compute_internal_restricted(cand_pos_t i, cand_pos_t j) {
    energy_t v_iloop = INF;
    cand_pos_t max_k = std::min(j - TURN - 2, i + MAXLOOP + 1);
    const int ptype_closing = pair[S_[i]][S_[j]];
    for (cand_pos_t k = i + 1; k <= max_k; ++k) {
        cand_pos_t min_l = std::max(k + TURN + 1 + MAXLOOP + 2, k + j - i) - MAXLOOP - 2;
        if ((tree->up[k - 1] >= (k - i - 1))) {
            for (cand_pos_t l = j - 1; l >= min_l; --l) {
                if (tree->up[j - 1] >= (j - l - 1)) {
                    energy_t v_iloop_kl = E_IntLoop(k - i - 1, j - l - 1, ptype_closing, rtype[pair[S_[k]][S_[l]]], S1_[i + 1], S1_[j - 1],
                                                    S1_[k - 1], S1_[l + 1], const_cast<vrna_param_t *>(params_))
                                          + get_energy(k, l);
                    if(i+1==k && j-1==l) v_iloop_kl += ShapeData->get_calculated(i) + ShapeData->get_calculated(j);
                    v_iloop = std::min(v_iloop, v_iloop_kl);
                }
            }
        }
    }
    return v_iloop;
}

void pseudo_loop::compute_energy_restricted(cand_pos_t i, cand_pos_t j)
// compute the V(i,j) value, if the structure must be restricted
{
    energy_t min, min_en[3];
    cand_pos_t k, min_rank;
    char type;

    min_rank = -1;
    min = INF / 2;
    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;

    const bool unpaired = (tree->tree[i].pair < -1 && tree->tree[j].pair < -1);
    const bool paired = (tree->tree[i].pair == j && tree->tree[j].pair == i);
    if (paired || unpaired) // if i and j can pair
    {
        bool canH = !(tree->up[j - 1] < (j - i - 1));
        if (canH) min_en[0] = HairpinE(seq_, i, j);

        min_en[1] = compute_internal_restricted(i, j);
        min_en[2] = compute_energy_VM_restricted(i, j);
    }
    for (k = 0; k < 3; k++) {
        if (min_en[k] < min) {
            min = min_en[k];
            min_rank = k;
        }
    }

    switch (min_rank) {
    case 0:
        type = HAIRP;
        break;
    case 1:
        type = INTER;
        break;
    case 2:
        type = MULTI;
        break;
    default:
        type = NONE;
    }

    if (min < INF / 2) {
        int ij = index[i] + j - i;
        V[ij].energy = min;
        V[ij].type = type;
    }
}

void pseudo_loop::compute_energies_PK(cand_pos_t i, cand_pos_t j) {
    cand_pos_t ij = index[i] + j - i;
    const pair_type ptype_closing = pair[S_[i]][S_[j]];
    bool weakly_closed_ij = tree->weakly_closed(i, j);
    // base cases:
    // a) i == j => VP[ij] = INF
    // b) [i,j] is a weakly_closed region => VP[ij] = INF
    // c) i or j is paired in original structure => VP[ij] = INF
    if (!(i == j || j - i < 4 || weakly_closed_ij)){
        if (ptype_closing > 0 && tree->tree[i].pair < -1 && tree->tree[j].pair < -1) compute_VP(i, j);
        if (tree->tree[j].pair < -1) compute_VPL(i, j);
        if (tree->tree[j].pair < j) compute_VPR(i, j);
    }

    if (!((j - i - 1) <= TURN || (tree->tree[i].pair >= -1 && tree->tree[i].pair > j) || (tree->tree[j].pair >= -1 && tree->tree[j].pair < i)
          || (tree->tree[i].pair >= -1 && tree->tree[i].pair < i) || (tree->tree[j].pair >= -1 && j < tree->tree[j].pair))) {
        compute_WMBW(i, j);
        compute_WMBP(i, j);
        compute_WMB(i, j);
    }

    if (!weakly_closed_ij) {
        WI[ij] = INF; // WIP is already INF
    } else {
        compute_WI(i, j);
        compute_WIP(i, j);
    }

    cand_pos_t ip = tree->tree[i].pair; // i's pair ip should be right side so ip = )
    cand_pos_t jp = tree->tree[j].pair; // j's pair jp should be left side so jp = (

    compute_BE(i, ip, jp, j);
}
// Added +1 to fres/tree indices as they are 1 ahead at the moment
void pseudo_loop::compute_WI(cand_pos_t i, cand_pos_t j) {
    energy_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;
    // branch 4, one base
    if (i == j) {
        WI.set(i,j) = PUP_penalty;
        return;
    }
    for (cand_pos_t k = i + 1; k < j - TURN - 1; ++k) {
        energy_t wi_1 = WI.get(i, k - 1);
        energy_t v_energy = wi_1 + get_energy(k, j);
        energy_t wmb_energy = wi_1 + WMB.get(k,j);
        m1 = std::min(m1, v_energy);
        m2 = std::min(m2, wmb_energy);
    }
    m1 += PPS_penalty;
    m2 += PSP_penalty + PPS_penalty;
    if (tree->tree[j].pair < 0) m3 = WI.get(i, j - 1) + PUP_penalty;
    m4 = get_energy(i, j) + PPS_penalty;
    m5 = WMB.get(i, j) + PSP_penalty + PPS_penalty;

    WI.set(i,j) = std::min({m1, m2, m3, m4, m5});
}

void pseudo_loop::compute_WIP(cand_pos_t i, cand_pos_t j) {
    energy_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF, m6 = INF, m7 = INF;
    // branch 1:
    for (cand_pos_t k = i + 1; k < j - TURN - 1; ++k) {
        bool can_pair = tree->up[k - 1] >= (k - i);
        energy_t wi_1 = WIP.get(i,k-1);
        energy_t v_energy = get_energy(k, j);
        energy_t wmb_energy = WMB.get(k,j);
        m1 = std::min(m1, wi_1 + v_energy);
        m2 = std::min(m2, wi_1 + wmb_energy);
        if (can_pair) m3 = std::min(m3, static_cast<energy_t>((k - i) * cp_penalty) + v_energy);
        if (can_pair) m4 = std::min(m4, static_cast<energy_t>((k - i) * cp_penalty) + wmb_energy);
    }
    m1 += bp_penalty;
    m2 += PSM_penalty + bp_penalty;
    m3 += bp_penalty;
    m4 += PSM_penalty + bp_penalty;
    // branch 2:
    if (tree->tree[j].pair < 0) m5 = WIP.get(i, j - 1) + cp_penalty;
    m6 = get_energy(i,j) + bp_penalty;
    m7 = WMB.get(i,j) + PSM_penalty + bp_penalty;

    WIP.set(i,j) = std::min({m1, m2, m3, m4, m5, m6, m7});
}

void pseudo_loop::compute_VPL(cand_pos_t i, cand_pos_t j) {
    energy_t m1 = INF;

    cand_pos_t min_Bp_j = std::min((cand_pos_tu)tree->b(i, j), (cand_pos_tu)tree->Bp(i, j));
    for (cand_pos_t k = i + 1; k < min_Bp_j; ++k) {
        bool can_pair = tree->up[k - 1] >= (k - i);
        if (can_pair) m1 = std::min(m1, static_cast<energy_t>((k - i) * cp_penalty) + VP.get(k,j));
    }
    VPL.set(i,j) = m1;
}

void pseudo_loop::compute_VPR(cand_pos_t i, cand_pos_t j) {
    energy_t m1 = INF, m2 = INF;

    cand_pos_t max_i_bp = std::max(tree->B(i, j), tree->bp(i, j));

    for (cand_pos_t k = max_i_bp + 1; k < j; ++k) {
        energy_t VP_energy = VP.get(i,k);
        bool can_pair = tree->up[j - 1] >= (j - k);

        m1 = std::min(m1, VP_energy + WIP.get(k + 1, j));
        if (can_pair) m2 = std::min(m2, VP_energy + static_cast<energy_t>((j - k) * cp_penalty));
    }
    VPR.set(i,j) = std::min(m1, m2);
}

void pseudo_loop::compute_VP(cand_pos_t i, cand_pos_t j) {
    energy_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF, m6 = INF, m7 = INF, m8 = INF, m9 = INF; // different branches

    // Borders -- added one to i and j to make it fit current bounds but also subtracted 1 from answer as the tree bounds are shifted as well
    cand_pos_t Bp_ij = tree->Bp(i, j);
    cand_pos_t B_ij = tree->B(i, j);
    cand_pos_t b_ij = tree->b(i, j);
    cand_pos_t bp_ij = tree->bp(i, j);

    //  1) inArc(i) and NOT_inArc(j)
    //  WI(i+1)(B'(i,j)-1)+WI(B(i,j)+1)(j-1)
    if ((tree->tree[i].parent->index) > 0 && (tree->tree[j].parent->index) < (tree->tree[i].parent->index) && Bp_ij >= 0 && B_ij >= 0 && bp_ij < 0) {
        energy_t WI_ipus1_BPminus = WI.get(i + 1, Bp_ij - 1);
        energy_t WI_Bplus_jminus = WI.get(B_ij + 1, j - 1);
        m1 = WI_ipus1_BPminus + WI_Bplus_jminus;
    }

    // 2) NOT_inArc(i) and inArc(j)
    // WI(i+1)(b(i,j)-1)+WI(b'(i,j)+1)(j-1)
    if ((tree->tree[i].parent->index) < (tree->tree[j].parent->index) && (tree->tree[j].parent->index) > 0 && b_ij >= 0 && bp_ij >= 0 && Bp_ij < 0) {
        energy_t WI_i_plus_b_minus = WI.get(i + 1, b_ij - 1);
        energy_t WI_bp_plus_j_minus = WI.get(bp_ij + 1, j - 1);
        m2 = WI_i_plus_b_minus + WI_bp_plus_j_minus;
    }

    // 3) inArc(i) and inArc(j)
    // WI(i+1)(B'(i,j)-1)+WI(B(i,j)+1)(b(i,j)-1)+WI(b'(i,j)+1)(j-1)
    if ((tree->tree[i].parent->index) > 0 && (tree->tree[j].parent->index) > 0 && Bp_ij >= 0 && B_ij >= 0 && b_ij >= 0 && bp_ij >= 0) {
        energy_t WI_i_plus_Bp_minus = WI.get(i + 1, Bp_ij - 1);
        energy_t WI_B_plus_b_minus = WI.get(B_ij + 1, b_ij - 1);
        energy_t WI_bp_plus_j_minus = WI.get(bp_ij + 1, j - 1);
        m3 = WI_i_plus_Bp_minus + WI_B_plus_b_minus + WI_bp_plus_j_minus;
    }

    // 4) NOT_paired(i+1) and NOT_paired(j-1) and they can pair together
    pair_type ptype_closingip1jm1 = pair[S_[i + 1]][S_[j - 1]];
    if ((tree->tree[i + 1].pair) < -1 && (tree->tree[j - 1].pair) < -1 && ptype_closingip1jm1 > 0) {
        m4 = get_e_stP(i, j) + VP.get(i + 1, j - 1) + ShapeData->get_calculated(i) + ShapeData->get_calculated(j);
    }

    // 5) NOT_paired(r) and NOT_paired(rp)
    //  VP(i,j) = e_intP(i,ip,jp,j) + VP(ip,jp)
    // Hosna, April 6th, 2007
    // whenever we use get_borders we have to check for the correct values
    cand_pos_t min_borders = std::min((cand_pos_tu)Bp_ij, (cand_pos_tu)b_ij);
    cand_pos_t edge_i = std::min(i + MAXLOOP + 1, j - TURN - 1);
    min_borders = std::min({min_borders, edge_i});
    for (cand_pos_t k = i + 1; k < min_borders; ++k) {
        if (tree->tree[k].pair < -1 && (tree->up[(k)-1] >= ((k) - (i)-1))) {
            cand_pos_t max_borders = std::max(bp_ij, B_ij) + 1;
            cand_pos_t edge_j = k + j - i - MAXLOOP - 2;
            max_borders = std::max({max_borders, edge_j});
            for (cand_pos_t l = j - 1; l > max_borders; --l) {

                pair_type ptype_closingkj = pair[S_[k]][S_[l]];
                if (tree->tree[l].pair < -1 && ptype_closingkj > 0 && (tree->up[(j)-1] >= ((j) - (l)-1))) {
                    // i and ip and j and jp should be in the same arc -- If it's unpaired between them, they have to be
                    energy_t tmp = get_e_intP(i, k, l, j) + VP.get(k,l);
                    m5 = std::min(m5, tmp);
                }
            }
        }
    }

    cand_pos_t min_Bp_j = std::min((cand_pos_tu)tree->b(i, j), (cand_pos_tu)tree->Bp(i, j));
    cand_pos_t max_i_bp = std::max(tree->B(i, j), tree->bp(i, j));

    for (cand_pos_t k = i + 1; k < min_Bp_j; ++k) {
        m6 = std::min(m6, WIP.get(i + 1, k - 1) + VP.get(k, j - 1));
    }
    m6 += ap_penalty + 2 * bp_penalty;

    for (cand_pos_t k = max_i_bp + 1; k < j; ++k) {
        m7 = std::min(m7, VP.get(i + 1, k) + WIP.get(k + 1, j - 1));
    }
    m7 += ap_penalty + 2 * bp_penalty;

    for (cand_pos_t k = i + 1; k < min_Bp_j; ++k) {
        m8 = std::min(m8, WIP.get(i + 1, k - 1) + VPR.get(k, j - 1));
    }
    m8 += ap_penalty + 2 * bp_penalty;

    for (cand_pos_t k = max_i_bp + 1; k < j; ++k) {
        m9 = std::min(m9, VPL.get(i + 1, k) + WIP.get(k + 1, j - 1));
    }
    m9 += ap_penalty + 2 * bp_penalty;

    // finding the min energy
    energy_t vp_h = std::min({m1, m2, m3});
    energy_t vp_iloop = std::min({m4, m5});
    energy_t vp_split = std::min({m6, m7, m8, m9});
    energy_t min = std::min({vp_h, vp_iloop, vp_split});

    VP.set(i,j) = min;
}

void pseudo_loop::compute_WMBW(cand_pos_t i, cand_pos_t j) {
    energy_t m1 = INF;
    if (tree->tree[j].pair < j) {
        for (cand_pos_t l = i + 1; l < j; l++) {
            if (tree->tree[l].pair < 0 && tree->tree[l].parent->index > -1 && tree->tree[j].parent->index > -1
                && tree->tree[j].parent->index == tree->tree[l].parent->index) {
                energy_t tmp = WMBP.get(i, l) + WI.get(l + 1, j);
                m1 = std::min(m1, tmp);
            }
        }
    }
    WMBW.set(i,j) = m1;
}

void pseudo_loop::compute_WMBP(cand_pos_t i, cand_pos_t j) {
    energy_t m1 = INF, m2 = INF, m4 = INF;
    // 1)
    if (tree->tree[j].pair < 0) {
        energy_t tmp = INF;
        cand_pos_t b_ij = tree->b(i, j);
        for (cand_pos_t l = i + 1; l < j; l++) {
            // Hosna: April 19th, 2007
            // the chosen l should be less than border_b(i,j) -- should be greater than border_b(i,l)
            // Mateo Jan 2025 Added exterior cases to consider when looking at band borders. Solved case of [.(.].[.).]
            int ext_case = compute_exterior_cases(l, j);
            if ((b_ij > 0 && l < b_ij) || (b_ij < 0 && ext_case == 0)) {
                cand_pos_t bp_il = tree->bp(i, l);
                cand_pos_t Bp_lj = tree->Bp(l, j);
                if (bp_il >= 0 && l > bp_il && Bp_lj > 0 && l < Bp_lj) { // bp(i,l) < l < Bp(l,j)
                    cand_pos_t B_lj = tree->B(l, j);
                    // Hosna: July 5th, 2007:
                    // as long as we have i <= arc(l)< j we are fine
                    if (i <= tree->tree[l].parent->index && tree->tree[l].parent->index < j && l + TURN <= j) {
                        energy_t sum = get_BE(tree->tree[B_lj].pair, B_lj, tree->tree[Bp_lj].pair, Bp_lj) + WMBP.get(i, l - 1) + VP.get(l, j);
                        tmp = std::min(tmp, sum);
                    }
                }
            }
            m1 = 2 * PB_penalty + tmp;
        }
    }
    // 2) WMB(i,j) = min_{i<l<j}{WMB(i,l)+WI(l+1,j)} if bp(j)<j
    // Hosna: Feb 5, 2007
    if (tree->tree[j].pair < 0) {
        energy_t tmp = INF;
        cand_pos_t b_ij = tree->b(i, j);
        for (cand_pos_t l = i + 1; l < j; l++) {
            // Hosna, April 6th, 2007
            // whenever we use get_borders we have to check for the correct values
            cand_pos_t bp_il = tree->bp(i, l);
            cand_pos_t Bp_lj = tree->Bp(l, j);
            // Hosna: April 19th, 2007
            // the chosen l should be less than border_b(i,j) -- should be greater than border_b(i,l)
            // Mateo Jan 2025 Added exterior cases to consider when looking at band borders. Solved case of [.(.].[.).]
            // Mateo May 2025 I'm adding to this -- exterior cases are for [.(.]..[.).] where b_ij = -2 when it should be inf and allow
            // the for loop to happen. If b_ij>0 though, the exterior cases shouldn't play a role.
            int ext_case = compute_exterior_cases(l, j);
            if ((b_ij > 0 && l < b_ij) || (b_ij < 0 && ext_case == 0)) {
                if (bp_il >= 0 && l > bp_il && Bp_lj > 0 && l < Bp_lj) { // bp(i,l) < l < Bp(l,j)
                    cand_pos_t B_lj = tree->B(l, j);
                    if (i <= tree->tree[l].parent->index && tree->tree[l].parent->index < j && l + TURN <= j) {
                        energy_t sum = get_BE(tree->tree[B_lj].pair, B_lj, tree->tree[Bp_lj].pair, Bp_lj) + WMBW.get(i, l - 1) + VP.get(l, j);
                        tmp = std::min(tmp, sum);
                    }
                }
            }
            m2 = 2 * PB_penalty + tmp;
        }
    }
    // 3) WMB(i,j) = VP(i,j) + P_b
    energy_t m3 = VP.get(i, j) + PB_penalty;
    // if not paired(j) and paired(i) then
    // WMBP(i,j) = 2*Pb + min_{i<l<bp(i)}(BE(i,bp(i),b'(i,l),bp(b'(i,l)))+WI(b'+1,l-1)+VP(l,j))
    if (tree->tree[j].pair < 0 && tree->tree[i].pair >= 0) {
        energy_t tmp = INF;
        // Hosna: June 29, 2007
        // if j is inside i's arc then the l should be
        // less than j not bp(i)
        // check with Anne
        for (cand_pos_t l = i + 1; l < j; l++) {
            // Hosna, April 9th, 2007
            // checking the borders as they may be negative
            cand_pos_t bp_il = tree->bp(i, l);
            if (bp_il >= 0 && bp_il < n && l + TURN <= j) {
                energy_t BE_energy = get_BE(i, tree->tree[i].pair, bp_il, tree->tree[bp_il].pair);
                energy_t WI_energy = WI.get(bp_il + 1, l - 1);
                energy_t VP_energy = VP.get(l, j);
                energy_t sum = BE_energy + WI_energy + VP_energy;
                tmp = std::min(tmp, sum);
            }
        }
        m4 = 2 * PB_penalty + tmp;
    }
    // get the min for WMB
    WMBP.set(i,j) = std::min({m1, m2, m3, m4});
}

void pseudo_loop::compute_WMB(cand_pos_t i, cand_pos_t j) {
    energy_t m2 = INF, mWMBP = INF;
    if (tree->tree[j].pair >= 0 && j > tree->tree[j].pair && tree->tree[j].pair > i) {
        cand_pos_t bp_j = tree->tree[j].pair;
        for (cand_pos_t l = (bp_j + 1); (l < j); l++) {
            // Hosna: April 24, 2007
            // correct case 2 such that a multi-pseudoknotted
            // loop would not be treated as case 2
            cand_pos_t Bp_lj = tree->Bp(l, j);

            if (Bp_lj >= 0 && Bp_lj < n) {
                energy_t sum = get_BE(bp_j, j, tree->tree[Bp_lj].pair, Bp_lj) + WMBP.get(i, l) + WI.get(l + 1, Bp_lj - 1);
                m2 = std::min(m2, sum);
            }
        }
        m2 += PB_penalty;
    }
    // check the WMBP value
    mWMBP = WMBP.get(i, j);
    // get the min for WMB
    WMB.set(i,j) = std::min(m2, mWMBP);
}

void pseudo_loop::compute_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp) {

    if (!(i >= 1 && i <= ip && ip < jp && jp <= j && j <= n && tree->tree[i].pair > 0 && tree->tree[j].pair > 0 && tree->tree[ip].pair > 0
          && tree->tree[jp].pair > 0 && tree->tree[i].pair == j && tree->tree[j].pair == i && tree->tree[ip].pair == jp
          && tree->tree[jp].pair == ip)) { // impossible cases
        return;
    }
    // base case: i.j and ip.jp must be in G
    if (tree->tree[i].pair != j || tree->tree[ip].pair != jp) {
        BE.set(i,ip) = INF;
        return;
    }

    // base case:
    if (i == ip && j == jp && i < j) {
        BE.set(i,ip) = 0;
        return;
    }

    energy_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;
    // 1) bp(i+1) == j-1
    if (tree->tree[i + 1].pair == j - 1) {
        m1 = get_e_stP(i, j) + get_BE(i + 1, j - 1, ip, jp) + ShapeData->get_calculated(i) + ShapeData->get_calculated(j);
    }

    // cases 2-5 are all need an l s.t. i<l<=ip and jp<=bp(l)<j
    for (cand_pos_t l = i + 1; l <= ip; l++) {

        // Hosna: March 14th, 2007
        if (tree->tree[l].pair >= -1 && jp <= tree->tree[l].pair && tree->tree[l].pair < j) {
            cand_pos_t lp = tree->tree[l].pair;

            bool empty_region_il = (tree->up[(l)-1] >= l - i - 1);       // empty between i+1 and lp-1
            bool empty_region_lpj = (tree->up[(j)-1] >= j - lp - 1);     // empty between l+1 and ip-1
            bool weakly_closed_il = tree->weakly_closed(i + 1, l - 1);   // weakly closed between i+1 and lp-1
            bool weakly_closed_lpj = tree->weakly_closed(lp + 1, j - 1); // weakly closed between l+1 and ip-1

            if (empty_region_il && empty_region_lpj) { //&& !(ip == (i+1) && jp==(j-1)) && !(l == (i+1) && lp == (j-1))){
                energy_t tmp = get_e_intP(i, l, lp, j) + get_BE(l, lp, ip, jp);
                m2 = std::min(m2, tmp);
            }
            // 3)
            if (weakly_closed_il && weakly_closed_lpj) {
                energy_t tmp = WIP.get(i + 1, l - 1) + get_BE(l, lp, ip, jp) + WIP.get(lp + 1, j - 1) + ap_penalty + 2 * bp_penalty;
                m3 = std::min(m3, tmp);
            }
            // 4)
            if (weakly_closed_il && empty_region_lpj) {
                energy_t tmp = WIP.get(i + 1, l - 1) + get_BE(l, lp, ip, jp) + cp_penalty * (j - lp - 1) + ap_penalty + 2 * bp_penalty;
                m4 = std::min(m4, tmp);
            }
            // 5)
            if (empty_region_il && weakly_closed_lpj) {
                energy_t tmp = ap_penalty + 2 * bp_penalty + (cp_penalty * (l - i - 1)) + get_BE(l, lp, ip, jp) + WIP.get(lp + 1, j - 1);
                m5 = std::min(m5, tmp);
            }
        }
    }
    // finding the min and putting it in BE[iip]
    BE.set(i,ip) = std::min({m1, m2, m3, m4, m5});
}
///////////////////////////// Util /////////////////////////////////
/**
 * @brief This code returns the hairpin energy for a given base pair.
 * @param i The left index in the base pair
 * @param j The right index in the base pair
 */
energy_t pseudo_loop::HairpinE(const std::string &seq, cand_pos_t i, cand_pos_t j) {

    const int ptype_closing = pair[S_[i]][S_[j]];

    if (ptype_closing == 0) return INF;
    return E_Hairpin(j - i - 1, ptype_closing, S1_[i + 1], S1_[j - 1], &seq.c_str()[i - 1], const_cast<vrna_param_t *>(params_));
}
/**
 * In cases where the band border is not found, if specific cases are met, the value is Inf(i.e n) not -1.
 * When applied to WMBP, if all cases are 0, then we can proceed with WMBP
 * Mateo Jan 2025: Added to Fix WMBP problem
 */
int pseudo_loop::compute_exterior_cases(cand_pos_t l, cand_pos_t j) {
    // Case 1 -> l is not covered
    bool case1 = tree->tree[l].parent->index <= 0;
    // Case 2 -> l is paired
    bool case2 = tree->tree[l].pair > 0;
    // Case 3 -> l is part of a closed subregion
    // bool case3 = 0;
    // Case 4 -> l.bp(l) i.e. l.j does not cross anything -- could I compare parents instead?
    bool case4 = j < tree->Bp(l, j);
    // By bitshifting each one, we have a more granular idea of what cases fail and is faster than branching
    return (case1 << 2) | (case2 << 1) | case4;
}
energy_t pseudo_loop::compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) {
    const pair_type ptype_closing = pair[S_[i]][S_[j]];
    return E_IntLoop(k - i - 1, j - l - 1, ptype_closing, rtype[pair[S_[k]][S_[l]]], S1_[i + 1], S1_[j - 1], S1_[k - 1], S1_[l + 1],
                     const_cast<vrna_param_t *>(params_));
}
energy_t pseudo_loop::get_e_stP(cand_pos_t i, cand_pos_t j) {
    if (i + 1 == j - 1) { // TODO: do I need something like that or stack is taking care of this?
        return INF;
    }
    energy_t ss = compute_int(i, j, i + 1, j - 1) + ShapeData->get_calculated(i) + ShapeData->get_calculated(j);
    return lrint(e_stP_penalty * ss);
}

energy_t pseudo_loop::get_e_intP(cand_pos_t i, cand_pos_t ip, cand_pos_t jp, cand_pos_t j) {
    energy_t e_int = compute_int(i, j, ip, jp);
    energy_t energy = lrint(e_intP_penalty * e_int);
    return energy;
}

/**
 * @brief Gives the W(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
*/
energy_t pseudo_loop::E_ext_Stem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1, const cand_pos_t i,const cand_pos_t j){

	energy_t e = INF;

    auto consider = [&](energy_t v, bool valid, pair_type tt, base_type s5, base_type s3) {
        if (!valid || v == INF) return;
        e = std::min(e, v + E_ExtLoop(tt, s5, s3, params_));
    };
	base_type si1  = i > 1 ? S_[i-1] : -1;
    base_type sj1  = j < n ? S_[j+1] : -1;
    base_type si = S_[i];
    base_type sj = S_[j];

	bool dangle2 = params_->model_details.dangles == 2;
    bool dangle1 = params_->model_details.dangles == 1;

	consider(vij, ((tree->tree[i].pair < -1 && tree->tree[j].pair < -1) || (tree->tree[i].pair == j && tree->tree[j].pair == i)), pair[S_[i]][S_[j]], dangle2 ? si1 : -1, dangle2 ? sj1 : -1);
	if (dangle1) {
        consider(vi1j,j-i-1>TURN && (((tree->tree[i + 1].pair < -1 && tree->tree[j].pair < -1) || (tree->tree[i + 1].pair == j)) && tree->tree[i].pair < 0), pair[S_[i+1]][S_[j]], si, -1);
        consider(vij1,j-1-i>TURN && (((tree->tree[i].pair < -1 && tree->tree[j - 1].pair < -1) || (tree->tree[i].pair == j - 1)) && tree->tree[j].pair < 0), pair[S_[i]][S_[j-1]], -1, sj);
        consider(vi1j1,j-1-i-1>TURN && (((tree->tree[i + 1].pair < -1 && tree->tree[j - 1].pair < -1) || (tree->tree[i + 1].pair == j - 1)) &&tree-> tree[i].pair < 0 && tree->tree[j].pair < 0), pair[S_[i+1]][S_[j-1]], si, sj);
    }
	return e;
}

/**
 * @brief Gives the WM(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
*/
energy_t pseudo_loop::E_MLStem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,cand_pos_t i, cand_pos_t j){

	energy_t e = INF;

    auto consider = [&](energy_t v, bool valid, pair_type type, base_type s5, base_type s3, int ml_count) {
        if (!valid || v == INF) return;
        e = std::min(e, v + E_MLstem(type, s5, s3, params_) + ml_count * params_->MLbase);
    };

	base_type si1  = i > 1 ? S_[i-1] : -1;
    base_type sj1  = j < n ? S_[j+1] : -1;
    base_type si = S_[i];
    base_type sj = S_[j];

	bool dangle2 = params_->model_details.dangles == 2;
    bool dangle1 = params_->model_details.dangles == 1;

	consider(vij, (tree->tree[i].pair < -1 && tree->tree[j].pair < -1) || (tree->tree[i].pair == j), pair[S_[i]][S_[j]], dangle2 ? si1 : -1, dangle2 ? sj1 : -1, 0);
	if (dangle1) {
		consider(vi1j,j-i-1>TURN && (((tree->tree[i + 1].pair < -1 && tree->tree[j].pair < -1) || (tree->tree[i + 1].pair == j)) && tree->tree[i].pair < 0), pair[S_[i+1]][S_[j]], si, -1, 1);
        consider(vij1,j-1-i>TURN && (((tree->tree[i].pair < -1 && tree->tree[j - 1].pair < -1) || (tree->tree[i].pair == j - 1)) && tree->tree[j].pair < 0), pair[S_[i]][S_[j-1]], -1, sj, 1);
        consider(vi1j1,j-1-i-1>TURN && (((tree->tree[i + 1].pair < -1 && tree->tree[j - 1].pair < -1) || (tree->tree[i + 1].pair == j - 1)) && tree->tree[i].pair < 0 && tree->tree[j].pair < 0), pair[S_[i+1]][S_[j-1]], si, sj, 2);
	}
    return e;
}

/**
* @brief Computes the multiloop V contribution. This gives back essentially VM(i,j).
* 
*/
energy_t pseudo_loop::E_MbLoop(const energy_t WM2ij, const energy_t WM2ip1j, const energy_t WM2ijm1, const energy_t WM2ip1jm1, cand_pos_t i, cand_pos_t j){
	energy_t e = INF;

    bool pairable = (tree->tree[i].pair < -1 && tree->tree[j].pair < -1) || (tree->tree[i].pair == j);
    pair_type tt = pair[S_[j]][S_[i]];
    base_type si1 = S_[i+1];
    base_type sj1 = S_[j-1];

	auto consider = [&](energy_t v, bool check, base_type s5, base_type s3, int ml_count) {
        if (check && v == INF) return;
        e = std::min(e, v + E_MLstem(tt, s5, s3, params_) + params_->MLclosing + ml_count * params_->MLbase);
    };

	bool dangle2 = params_->model_details.dangles == 2;
    bool dangle1 = params_->model_details.dangles == 1;

	consider(WM2ij,pairable, dangle2 ? sj1 : -1, dangle2 ? si1 : -1, 0);
	if(dangle1){
		// ML pair 5 — closing (i,j) with mb part [i+2, j-1]
		consider(WM2ip1j,pairable && tree->tree[i+1].pair < 0, -1, si1, 1);
        // ML pair 3 — closing (i,j) with mb part [i+1, j-2]
        consider(WM2ijm1,pairable && tree->tree[j-1].pair < 0, sj1, -1, 1);
        // ML pair 53 — closing (i,j) with mb part [i+2, j-2]
        consider(WM2ip1jm1,pairable && tree->tree[i+1].pair < 0 && tree->tree[j-1].pair < 0, sj1, si1, 2);
	}
	return e;
}






// Mateo 13 Sept 2023
// return number of bases in between the two inclusive index
int distance(int left, int right) { return (right - left - 1); }

// Mateo 13 Sept 2023
// given a initial hotspot which is a hairpin loop, keep trying to add a arc to form a larger stack
void expand_hotspot(s_energy_matrix *V, Hotspot &hotspot, int n) {
    // printf("\nexpanding hotspot: i: %d j: %d\n",hotspot->get_left_inner_index(),hotspot->get_right_inner_index());
    // calculation for the hairpin that is already in there
    V->compute_hotspot_energy(hotspot.get_left_outer_index(), hotspot.get_right_outer_index(), 0);

    // try to expand by adding a arc right beside the current out most arc
    while (hotspot.get_left_outer_index() - 1 >= 1 && hotspot.get_right_outer_index() + 1 <= n) {
        base_type sim1 = V->S_[hotspot.get_left_outer_index() - 1];
        base_type sjp1 = V->S_[hotspot.get_right_outer_index() + 1];
        pair_type ptype_closing = pair[sim1][sjp1];
        if (ptype_closing > 0) {
            hotspot.move_left_outer_index();
            hotspot.move_right_outer_index();
            hotspot.increment_size();
            V->compute_hotspot_energy(hotspot.get_left_outer_index(), hotspot.get_right_outer_index(), 1);
        } else {
            break;
        }
    }
    base_type i = hotspot.get_left_outer_index();
    base_type j = hotspot.get_right_outer_index();
    pair_type tt = pair[V->S_[i]][V->S_[j]];
    base_type si1 = i > 1 ? V->S_[i - 1] : -1;
    base_type sj1 = j < n ? V->S_[j + 1] : -1;
    energy_t dangle_penalty = E_ExtLoop(tt, si1, sj1, V->params_);

    double energy = V->get_energy(hotspot.get_left_outer_index(), hotspot.get_right_outer_index());

    // printf("here and %d\n",energy);
    // printf("energy: %lf, AU_total: %d, dangle_top_total: %d, dangle_bot_total: %d\n",energy,non_gc_penalty,dangle_top_penalty,dangle_bot_penalty);

    energy = (energy + dangle_penalty) / 100;

    hotspot.set_energy(energy);
    // printf("done: %d %d %d
    // %d\n",hotspot->get_left_outer_index(),hotspot->get_left_inner_index(),hotspot->get_right_inner_index(),hotspot->get_right_outer_index());
    return;
}

// Mateo 13 Sept 2023
// look for every possible hairpin loop, and try to add a arc to form a larger stack with at least min_stack_size bases
void get_hotspots(std::string seq, std::vector<Hotspot> &hotspot_list, SHAPEData &ShapeData, int max_hotspot, vrna_param_s *params) {

    int n = seq.length();
    s_energy_matrix *V;
    make_pair_matrix();
    short *S_ = encode_sequence(seq.c_str(), 0);
    short *S1_ = encode_sequence(seq.c_str(), 1);
    V = new s_energy_matrix(seq, n, &ShapeData, S_, S1_, params);
    int min_bp_distance = 3;
    int min_stack_size = 3; // the hotspot must be a stack of size >= 3
    // Hotspot current_hotspot;
    // start at min_stack_size-1 and go outward to try to add more arcs to form bigger stack because we cannot expand more than min_stack_size from
    // there anyway
    for (int i = min_stack_size; i <= n; i++) {
        for (int j = i; j <= n; j++) {
            int ptype_closing = pair[V->S_[i]][V->S_[j]];
            if (ptype_closing > 0 && distance(i, j) >= min_bp_distance) {
                // current_hotspot = new Hotspot(i,j,nb_nucleotides);

                Hotspot current_hotspot(i, j, n);

                expand_hotspot(V, current_hotspot, n);

                if (current_hotspot.get_size() < min_stack_size || current_hotspot.is_invalid_energy()) {

                } else {

                    current_hotspot.set_structure();
                    hotspot_list.push_back(current_hotspot);
                }
            }
        }
    }

    // make sure we only keep top 20 hotspot with lowest energy
    std::sort(hotspot_list.begin(), hotspot_list.end(), compare_hotspot_ptr);
    while ((int)hotspot_list.size() > max_hotspot) {
        hotspot_list.pop_back();
    }

    // if no hotspot found, add all _ as restricted
    if ((int)hotspot_list.size() == 0) {
        Hotspot hotspot(1, n, n + 1);
        hotspot.set_default_structure();
        hotspot_list.push_back(hotspot);
    }
    delete V;
    free(S_);
    free(S1_);

    return;
}

bool compare_hotspot_ptr(Hotspot &a, Hotspot &b) {
    if (a.get_energy() != b.get_energy()) {
        return a.get_energy() < b.get_energy();
    }
    return a.get_structure() < b.get_structure();
}