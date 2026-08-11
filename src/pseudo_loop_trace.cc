#include "pseudo_loop.hh"
#include "h_globals.hh"
#include <algorithm>
#include <iostream>
#include <string>



////////////////////////// Traceback ////////////////////////////////
void pseudo_loop::backtrack(){
   Trace_W(1,n,W[n]);
   return;
}

void pseudo_loop::Trace_W(cand_pos_t i, cand_pos_t j, energy_t e){
	if (debug) printf("W at %d and %d with %d\n", i, j, e);
	if (j<=i) return;

	energy_t acc = INF;

	// this case is for j unpaired, so I have to check that.
	energy_t tmp = W[j-1];
	if (e==tmp){
		Trace_W(i,j-1,W[j-1]);
		return;
	}
	for (cand_pos_t i=1; i<=j-1; i++){
		acc = (i>1) ? W[i-1] : 0;
		base_type si1 = i>1 ? S_[i-1] : -1;
		base_type sj1 = j<n ? S_[j+1] : -1;
		tmp = acc + get_energy(i,j) + ((params_->model_details.dangles == 2) ? E_ExtLoop(pair[S_[i]][S_[j]],si1,sj1,params_) : E_ExtLoop(pair[S_[i]][S_[j]],-1,-1,params_));
		if(e==tmp){
			Trace_W(1,i-1,W[i-1]);
			Trace_V(i,j,get_energy(i,j));
			return;
		}
		if(params_->model_details.dangles ==1){
			tmp = acc + get_energy(i+1,j) + E_ExtLoop(pair[S_[i+1]][S_[j]],S_[i],-1,params_);
			if(e==tmp){
				Trace_W(1,i-1,W[i-1]);
				Trace_V(i+1,j,get_energy(i+1,j));
				return;
			}
			tmp = acc + get_energy(i,j-1) + E_ExtLoop(pair[S_[i]][S_[j-1]],-1,S_[j],params_);
			if(e==tmp){
				Trace_W(1,i-1,W[i-1]);
				Trace_V(i,j-1,get_energy(i,j-1));
				return;
			}
			tmp = acc + get_energy(i+1,j-1) + E_ExtLoop(pair[S_[i+1]][S_[j-1]],S_[i],S_[j],params_);
			if(e==tmp){
				Trace_W(1,i-1,W[i-1]);
				Trace_V(i+1,j-1,get_energy(i+1,j-1));
				return;
			}
		}
	}
	for (cand_pos_t i=1; i<=j-1; i++){
		acc = (i-1>0) ? W[i-1]: 0;
		tmp = acc + WMB.get(i,j)+ PS_penalty;
		if(e==tmp){
			Trace_W(1,i-1,W[i-1]);
			Trace_WMB(i,j,WMB.get(i,j));
			return;
		}
	}
	UNREACHABLE();
}
void pseudo_loop::Trace_V(cand_pos_t i, cand_pos_t j, energy_t e){
	if (debug) printf("V at %d and %d as type: %c with %d\n", i, j,get_type(i,j), e);
	structure[i] = '(';
	structure[j] = ')';
	char type = get_type(i,j);

	switch(type){
		case HAIRP:{
			return;
		}
		case INTER:{
			cand_pos_t max_k = std::min(j-TURN-2,i+MAXLOOP+1);
			for (cand_pos_t k = i+1; k <= max_k; ++k){
				cand_pos_t min_l=std::max(k+TURN+1 + MAXLOOP+2, k+j-i) - MAXLOOP-2;
				for (cand_pos_t l = j-1; l >= min_l; --l)
				{
					energy_t tmp = compute_int(i,j,k,l) + get_energy(k,l);
                    if(i+1==k && j-1==l) tmp += ShapeData->get_calculated(i) + ShapeData->get_calculated(j);
					if (e == tmp)
					{
						Trace_V(k,l,get_energy(k,l));
						return;
					}
				}
				
			}
		}
		break;
		case MULTI: {
			energy_t tmp = INF;
			for (cand_pos_t k = i+1; k <= j-1; ++k){
                if(i==58 && j==334) printf("k is %d and e is %d and WM is %d and WMv is %d and total is %d\n",k,e,WM.get(i+1,k-1),WMv.get(k,j-1),WM.get(i+1,k-1) + std::min(WMv.get(k,j-1),WMp.get(k,j-1)) + params_->MLclosing + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_));
				tmp = WM.get(i+1,k-1) + std::min(WMv.get(k,j-1),WMp.get(k,j-1)) + params_->MLclosing;
				if(params_->model_details.dangles == 2){
					tmp += E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_);
				} else {
					tmp += E_MLstem(pair[S_[j]][S_[i]],-1,-1,params_);
				}
				if (e==tmp){
					tmp -= params_->MLclosing;
					tmp -= (params_->model_details.dangles == 2 ? E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_) : E_MLstem(pair[S_[j]][S_[i]],-1,-1,params_));
					Trace_WM(i+1,k-1,WM.get(i+1,k-1));
					if(tmp == WM.get(i+1,k-1) + WMv.get(k,j-1)){
						Trace_WMv(k,j-1,WMv.get(k, j-1));
					} else {
						Trace_WMp(k,j-1,WMp.get(k, j-1));
					}
					return;
				}
				tmp = static_cast<energy_t>((k-i-1)*params_->MLbase + WMp.get(k,j-1)) + params_->MLclosing;
                tmp += (params_->model_details.dangles == 2 ? E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_) : E_MLstem(pair[S_[j]][S_[i]],-1,-1,params_));
				if (e==tmp){
					Trace_WMp(k,j-1,WMp.get(k,j-1));
					return;
				}
				if(params_->model_details.dangles ==1){
					tmp = WM.get(i+2,k-1) + std::min(WMv.get(k,j-1),WMp.get(k,j-1)) + E_MLstem(pair[S_[j]][S_[i]],-1,S_[i+1],params_) + params_->MLclosing + params_->MLbase;
					if (e==tmp){
						tmp -= (params_->MLclosing + E_MLstem(pair[S_[j]][S_[i]],-1,S_[i+1],params_) + params_->MLbase);
						Trace_WM(i+2,k-1,WM.get(i+2,k-1));
						if(tmp == WM.get(i+2,k-1) + WMv.get(k,j-1)){
							Trace_WMv(k,j-1,WMv.get(k, j-1));
						} else {
							Trace_WMp(k,j-1,WMv.get(k, j-1));
						}
						return;
					}
					tmp = WM.get(i+1,k-1) + std::min(WMv.get(k,j-2),WMp.get(k,j-2)) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],-1,params_) + params_->MLclosing + params_->MLbase;
					if (e==tmp){
						tmp -= (params_->MLclosing + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],-1,params_) + params_->MLbase);
						Trace_WM(i+1,k-1,WM.get(i+1,k-1));
						if(tmp == WM.get(i+1,k-1) + WMv.get(k,j-2)){
							Trace_WMv(k,j-2,WMv.get(k, j-2));
						} else {
							Trace_WMp(k,j-2,WMv.get(k, j-2));
						}
						return;
					}
					tmp = WM.get(i+2,k-1) + std::min(WMv.get(k,j-2),WMp.get(k,j-2)) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_) + params_->MLclosing + 2*params_->MLbase;
					if (e==tmp){
						tmp -= (params_->MLclosing + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_) + 2*params_->MLbase);
						Trace_WM(i+2,k-1,WM.get(i+2,k-1));
						if(tmp == WM.get(i+2,k-1) + WMv.get(k,j-2)){
							Trace_WMv(k,j-2,WMv.get(k, j-2));
						} else {
							Trace_WMp(k,j-2,WMv.get(k, j-2));
						}
						return;
					}
					if((k-(i+1)-1) >=0) tmp = static_cast<energy_t>((k-(i+1)-1)*params_->MLbase) + WMp.get(k,j-1) + E_MLstem(pair[S_[j]][S_[i]],-1,S_[i+1],params_) + params_->MLclosing + params_->MLbase;
					if (e==tmp){
						Trace_WMp(k,j-1,WMp.get(k, j-1));
						return;
					}
					tmp = static_cast<energy_t>((k-i-1)*params_->MLbase) + WMp.get(k,j-2) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],-1,params_) + params_->MLclosing + params_->MLbase;
					if (e==tmp){
						Trace_WMp(k,j-2,WMp.get(k, j-2));
						return;
					}
					if((k-(i+1)-1) >=0) tmp = static_cast<energy_t>((k-(i+1)-1)*params_->MLbase) + WMp.get(k,j-2) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_) + params_->MLclosing + 2*params_->MLbase;
					if (e==tmp){
						Trace_WMp(k,j-2,WMp.get(k, j-2));
						return;
					}
				}				
			}
		}
		break;
	}
	UNREACHABLE();
}
void pseudo_loop::Trace_WM(cand_pos_t i, cand_pos_t j, energy_t e){
	if (debug) printf("WM at %d and %d with %d\n", i, j, e);
	energy_t tmp = INF;

	tmp = WM.get(i,j-1)+params_->MLbase;
	if(e==tmp){
		Trace_WM(i,j-1,WM.get(i,j-1));
		return;
	}
	for (cand_pos_t k=i; k <= j-TURN-1; k++){	
		tmp = static_cast<energy_t>((k-i)*params_->MLbase) + WMv.get(k,j);
		if(e==tmp){
			Trace_WMv(k,j,WMv.get(k,j));
			return;
		}
		tmp = static_cast<energy_t>((k-i)*params_->MLbase) + WMp.get(k,j);
		if(e==tmp){
			Trace_WMp(k,j,WMp.get(k,j));
			return;
		}
		tmp = WM.get(i,k-1) + WMv.get(k,j);
		if(e==tmp){
			Trace_WM(i,k-1,WM.get(i, k-1));
			Trace_WMv(k,j,WMv.get(k,j));
			return;
		}
		tmp = WM.get(i,k-1) + WMp.get(k,j);
		if(e==tmp){
			Trace_WM(i,k-1,WM.get(i, k-1));
			Trace_WMp(k,j,WMp.get(k,j));
			return;
		}
	}
	UNREACHABLE();
}
void pseudo_loop::Trace_WMv(cand_pos_t i, cand_pos_t j, energy_t e){
	if (debug) printf("WMv at %d and %d with %d\n", i, j, e);
	cand_pos_t si = S_[i];
	cand_pos_t sj = S_[j];
	cand_pos_t si1 = (i>1) ? S_[i-1] : -1;
	cand_pos_t sj1 = (j<n) ? S_[j+1] : -1;
	pair_type tt = pair[S_[i]][S_[j]];
	energy_t tmp = get_energy(i,j) + ((params_->model_details.dangles == 2) ? E_MLstem(tt,si1,sj1,params_) : E_MLstem(tt,-1,-1,params_));
	if(e==tmp){
		Trace_V(i,j,get_energy(i,j));
		return;
	}

	if(params_->model_details.dangles == 1){
		tt = pair[S_[i+1]][S_[j]];
		energy_t tmp = get_energy(i+1,j) + E_MLstem(tt,si,-1,params_) + params_->MLbase;
		if(e==tmp){
			Trace_V(i+1,j,get_energy(i+1,j));
			return;
		}
		tt = pair[S_[i]][S_[j-1]];
		tmp = get_energy(i,j-1) + E_MLstem(tt,-1,sj,params_) + params_->MLbase;
		if(e==tmp){
			Trace_V(i,j-1,get_energy(i,j-1));
			return;
		}
		tt = pair[S_[i+1]][S_[j-1]];
		tmp = get_energy(i+1,j-1) + E_MLstem(tt,si,sj,params_) + 2*params_->MLbase;
		if(e==tmp){
			Trace_V(i+1,j-1,get_energy(i+1,j-1));
			return;
		}
	}

	tmp = WMv.get(i,j-1) + params_->MLbase;
	if(e==tmp){
		Trace_WMv(i,j-1,WMv.get(i,j-1));
		return;
	}
	UNREACHABLE();
}
void pseudo_loop::Trace_WMp(cand_pos_t i, cand_pos_t j, energy_t e){
	if (debug) printf("WMp at %d and %d with %d\n", i, j, e);
	energy_t tmp = WMB.get(i,j) + PSM_penalty + b_penalty;
	if(e==tmp){
		Trace_WMB(i,j,WMB.get(i,j));
		return;
	}
	tmp = WMp.get(i,j-1) + params_->MLbase;
	if(e==tmp){
		Trace_WMp(i,j-1,WMp.get(i,j-1));
		return;
	}
	UNREACHABLE();
}

void pseudo_loop::Trace_WI(cand_pos_t i, cand_pos_t j, energy_t e){
    if (debug) printf("WI at %d and %d with %d\n", i, j, e);
    if (i >= j) return;
    if (tree->tree[j].pair < 0 && e == WI.get(i,j-1) + PUP_penalty) {
        Trace_WI(i,j-1,WI.get(i,j-1));
        return;
    }
    if(e==get_energy(i, j) + PPS_penalty){
        Trace_V(i,j,get_energy(i, j));
        return;
    };
    if (e == WMB.get(i,j) + PSP_penalty + PPS_penalty) {
        Trace_WMB(i,j,WMB.get(i,j));
        return;
    }
    for (cand_pos_t k = i+1; k < j; ++k) {
        if (e == WI.get(i,k-1) + get_energy(k,j) + PPS_penalty) {
            Trace_WI(i,k-1,WI.get(i,k-1));
            Trace_V(k,j,get_energy(k,j));
            return;
        }
    }
    for (cand_pos_t k = i+1; k < j; ++k) {
        if (e == WI.get(i,k-1) + WMB.get(k,j) + PSP_penalty + PPS_penalty) {
            Trace_WI(i,k-1,WI.get(i,k-1));
            Trace_WMB(k,j,WMB.get(k,j));
            return;
        }
    }
    UNREACHABLE();
}

void pseudo_loop::Trace_WIP(cand_pos_t i, cand_pos_t j, energy_t e){
    if (debug) printf("WIP at %d and %d with %d\n", i, j, e);
    if (tree->tree[j].pair < 0 && e == WIP.get(i,j-1) + cp_penalty) {
        Trace_WIP(i,j-1,WIP.get(i,j-1));
        return;
    }
    if (e == get_energy(i,j) + bp_penalty) {
        Trace_V(i,j,get_energy(i,j));
        return;
    }
    if (e == WMB.get(i,j) + PSM_penalty + bp_penalty) {
        Trace_WMB(i,j,WMB.get(i,j));
        return;
    }
    for (cand_pos_t k = i + 1; k < j - TURN - 1; ++k) {
        if (e == WIP.get(i,k-1) + get_energy(k,j) + bp_penalty) {
            Trace_WIP(i,k-1,WIP.get(i,k-1));
            Trace_V(k,j,get_energy(k,j));
            return;
        }
    }
    for (cand_pos_t k = i + 1; k < j - TURN - 1; ++k) {
        if (e == WIP.get(i,k-1) + WMB.get(k,j) + PSM_penalty + bp_penalty) {
            Trace_WIP(i,k-1,WIP.get(i,k-1));
            Trace_WMB(k,j,WMB.get(k,j));
            return;
        }
    }
    for (cand_pos_t k = i + 1; k < j - TURN - 1; ++k) {
        bool can_pair = tree->up[k - 1] >= (k - i);
        if (can_pair && e == static_cast<energy_t>((k-i)*cp_penalty) + get_energy(k,j) + bp_penalty){
            Trace_V(k,j,get_energy(k,j));
            return;
        }
    }
    for (cand_pos_t k = i + 1; k < j - TURN - 1; ++k) {
        bool can_pair = tree->up[k - 1] >= (k - i);
        if (can_pair && e == static_cast<energy_t>((k-i)*cp_penalty) + WMB.get(k,j) + PSM_penalty + bp_penalty){
            Trace_WMB(k,j,WMB.get(k,j));
            return;
        }
    }
    UNREACHABLE();
}
void pseudo_loop::Trace_WMB(cand_pos_t i, cand_pos_t j, energy_t e){
    if (debug) printf("WMB at %d and %d with %d\n", i, j, e);
    if(e == WMBP.get(i,j)){
        Trace_WMBP(i,j,WMBP.get(i,j));
        return;
    }
    if (tree->tree[j].pair >= 0 && j > tree->tree[j].pair && tree->tree[j].pair > i) {
        cand_pos_t bp_j = tree->tree[j].pair;
        for (cand_pos_t l = bp_j + 1; l < j; l++) {
            cand_pos_t Bp_lj = tree->Bp(l, j);
            if (Bp_lj >= 0 && Bp_lj < n) {
                if (e == get_BE(bp_j,j,tree->tree[Bp_lj].pair,Bp_lj) + WMBP.get(i,l) + WI.get(l+1,Bp_lj-1) + PB_penalty) {
                    Trace_BE(bp_j,j,tree->tree[Bp_lj].pair, Bp_lj,BE.get(bp_j,tree->tree[Bp_lj].pair));
                    Trace_WMBP(i,l,WMBP.get(i,l));
                    Trace_WI(l+1,Bp_lj-1,WI.get(l+1,Bp_lj-1));
                    return;   
                }
            }
        }
    }
    if(i==391 && j==444) printf("i is %d and j is %d and e is %d and WMBP is %d\n",i,j,e,WMBP.get(i,j));
    UNREACHABLE();
}
void pseudo_loop::Trace_WMBP(cand_pos_t i, cand_pos_t j, energy_t e){
    if (debug) printf("WMBP at %d and %d with %d\n", i, j, e);
    if (e == VP.get(i,j) + PB_penalty) {
        Trace_VP(i,j,VP.get(i,j));
        return;
    }
    if (tree->tree[j].pair < 0) {
        cand_pos_t b_ij = tree->b(i, j);
        for (cand_pos_t l = i + 1; l < j; l++) {
            cand_pos_t bp_il = tree->bp(i, l);
            cand_pos_t Bp_lj = tree->Bp(l, j);
            // Mateo Jan 2025 Added exterior cases to consider when looking at band borders. Solved case of [.(.].[.).]
            int ext_case = compute_exterior_cases(l, j);
            if ((b_ij > 0 && l < b_ij) || (b_ij < 0 && ext_case == 0)) {
                if (bp_il >= 0 && l > bp_il && Bp_lj > 0 && l < Bp_lj) { // bp(i,l) < l < Bp(l,j)
                    cand_pos_t B_lj = tree->B(l, j);
                    if (i <= tree->tree[l].parent->index && tree->tree[l].parent->index < j && l + TURN <= j) {
                        if(e == get_BE(tree->tree[B_lj].pair, B_lj, tree->tree[Bp_lj].pair, Bp_lj) + WMBP.get(i,l-1) + VP.get(l,j) + 2*PB_penalty){
                            Trace_BE(tree->tree[B_lj].pair, B_lj, tree->tree[Bp_lj].pair, Bp_lj,get_BE(tree->tree[B_lj].pair, B_lj, tree->tree[Bp_lj].pair, Bp_lj));
                            Trace_WMBP(i,l-1,WMBP.get(i,l-1));
                            Trace_VP(l,j,VP.get(l,j));
                            return;
                        }
                    }
                }
            }
        }
    }
    if (tree->tree[j].pair < 0) {
        cand_pos_t b_ij = tree->b(i, j);
        for (cand_pos_t l = i + 1; l < j; l++) {
            cand_pos_t bp_il = tree->bp(i, l);
            cand_pos_t Bp_lj = tree->Bp(l, j);
            // Mateo Jan 2025 Added exterior cases to consider when looking at band borders. Solved case of [.(.].[.).]
            int ext_case = compute_exterior_cases(l, j);
            if ((b_ij > 0 && l < b_ij) || (b_ij < 0 && ext_case == 0)) {
                if (bp_il >= 0 && l > bp_il && Bp_lj > 0 && l < Bp_lj) { // bp(i,l) < l < Bp(l,j)
                    cand_pos_t B_lj = tree->B(l, j);
                    if (i <= tree->tree[l].parent->index && tree->tree[l].parent->index < j && l + TURN <= j) {
                        if(e == get_BE(tree->tree[B_lj].pair, B_lj, tree->tree[Bp_lj].pair, Bp_lj) + WMBW.get(i,l-1) + VP.get(l,j) + 2*PB_penalty){
                            Trace_BE(tree->tree[B_lj].pair, B_lj, tree->tree[Bp_lj].pair, Bp_lj,get_BE(tree->tree[B_lj].pair, B_lj, tree->tree[Bp_lj].pair, Bp_lj));
                            Trace_WMBW(i,l-1,WMBW.get(i,l-1));
                            Trace_VP(l,j,VP.get(l,j));
                            return;
                        }
                    }
                }
            }
        }
    }
    if (tree->tree[j].pair < 0 && tree->tree[i].pair >= 0) {
            for (cand_pos_t l = i + 1; l < j; l++) {
                cand_pos_t bp_il = tree->bp(i, l);
                if (bp_il >= 0 && bp_il < n && l + TURN <= j) {
                    if (e == get_BE(i,tree->tree[i].pair, bp_il, tree->tree[bp_il].pair) + WI.get(bp_il+1,l-1) + VP.get(l,j) + 2*PB_penalty) {
                        Trace_BE(i,tree->tree[i].pair, bp_il, tree->tree[bp_il].pair,get_BE(i,tree->tree[i].pair, bp_il, tree->tree[bp_il].pair));
                        Trace_WI(bp_il+1,l-1,WI.get(bp_il+1,l-1));
                        Trace_VP(l,j,VP.get(l,j));
                        return;
                    }
                }
            }
        }

    UNREACHABLE();
}
void pseudo_loop::Trace_WMBW(cand_pos_t i, cand_pos_t j, energy_t e){
    if (debug) printf("WMBW at %d and %d with %d\n", i, j, e);
    if (tree->tree[j].pair < j) {
        for (cand_pos_t l = i + 1; l < j; l++) {
            if (tree->tree[l].pair < 0 && tree->tree[l].parent->index > -1 && tree->tree[j].parent->index > -1
                && tree->tree[j].parent->index == tree->tree[l].parent->index) {
                if (e == WMBP.get(i,l) + WI.get(l+1,j)) {
                    Trace_WMBP(i,l,WMBP.get(i,l));
                    Trace_WI(l+1,j,WI.get(l+1,j));
                    return;
                }
            }
        }
    }

    UNREACHABLE();
}
void pseudo_loop::Trace_VP(cand_pos_t i, cand_pos_t j, energy_t e){
    if (debug) printf("VP at %d and %d with %d\n", i, j, e);
    structure[i] = '[';
    structure[j] = ']';
    cand_pos_t Bp_ij = tree->Bp(i, j), B_ij = tree->B(i, j), b_ij = tree->b(i, j), bp_ij = tree->bp(i, j);

    if ((tree->tree[i].parent->index) > 0 && (tree->tree[j].parent->index) < (tree->tree[i].parent->index) && Bp_ij >= 0 && B_ij >= 0 && bp_ij < 0) {
        if (e==WI.get(i+1,Bp_ij-1)+WI.get(B_ij+1,j-1)) {
            Trace_WI(i+1,Bp_ij-1,WI.get(i+1,Bp_ij-1));
            Trace_WI(B_ij+1,j-1,WI.get(B_ij+1,j-1));
            return;
        }
    }
    if ((tree->tree[i].parent->index) < (tree->tree[j].parent->index) && (tree->tree[j].parent->index) > 0 && b_ij >= 0 && bp_ij >= 0 && Bp_ij < 0) {
        if (e == WI.get(i+1,b_ij-1)+WI.get(bp_ij+1,j-1)) {
            Trace_WI(i+1,b_ij-1,WI.get(i+1,b_ij-1));
            Trace_WI(bp_ij+1,j-1,WI.get(bp_ij+1,j-1));
            return;
        }
    }
    if ((tree->tree[i].parent->index) > 0 && (tree->tree[j].parent->index) > 0 && Bp_ij >= 0 && B_ij >= 0 && b_ij >= 0 && bp_ij >= 0) {
        if (e == WI.get(i+1,Bp_ij-1)+WI.get(B_ij+1,b_ij-1)+WI.get(bp_ij+1,j-1)) {
            Trace_WI(i+1,Bp_ij-1,WI.get(i+1,Bp_ij-1));
            Trace_WI(B_ij+1,b_ij-1,WI.get(B_ij+1,b_ij-1));
            Trace_WI(bp_ij+1,j-1,WI.get(bp_ij+1,j-1));
            return;
        }
    }
    pair_type ptype_closingip1jm1 = pair[S_[i+1]][S_[j-1]];
    if (tree->tree[i+1].pair < 0 && tree->tree[j-1].pair < 0 && ptype_closingip1jm1 > 0) {
        if (e == (energy_t) (get_e_stP(i,j) + VP.get(i+1,j-1) + ShapeData->get_calculated(i) + ShapeData->get_calculated(j))) {
            Trace_VP(i+1,j-1,VP.get(i+1,j-1));
            return;
        }
    }
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
                if (tree->tree[l].pair < -1 && ptype_closingkj > 0 && (tree->up[(j)-1] >= ((j)-(l)-1))) {
                    if (e == get_e_intP(i,k,l,j) + VP.get(k,l)) {
                        Trace_VP(k,l,VP.get(k,l));
                        return;
                    }
                }
            }
        }
    }
    cand_pos_t min_Bp_j = std::min((cand_pos_tu)tree->b(i, j),(cand_pos_tu)tree->Bp(i, j));
    cand_pos_t max_i_bp = std::max(tree->B(i, j),tree->bp(i, j));

    for (cand_pos_t k = i + 1; k < min_Bp_j; ++k) {
        if (e == WIP.get(i+1,k-1) + VP.get(k,j-1) + ap_penalty + 2 * bp_penalty) {
            Trace_WIP(i+1,k-1,WIP.get(i+1,k-1));
            Trace_VP(k,j-1,VP.get(k,j-1));
            return;
        }
    }
    for (cand_pos_t k = max_i_bp + 1; k < j; ++k) {
        if (e == VP.get(i+1,k) + WIP.get(k+1,j-1) + ap_penalty + 2 * bp_penalty) {
            Trace_VP(i+1,k,VP.get(i+1,k));
            Trace_WIP(k+1,j-1,WIP.get(k+1,j-1));
            return;
        }
    }
    for (cand_pos_t k = i+1; k < min_Bp_j; ++k) {
        if (e == WIP.get(i+1,k-1) + VPR.get(k,j-1) + ap_penalty + 2 * bp_penalty) {
            Trace_WIP(i+1,k-1,WIP.get(i+1,k-1));
            Trace_VPR(k,j-1,VPR.get(k,j-1));
            return;
        }
    }
    for (cand_pos_t k = max_i_bp + 1; k < j; ++k) {
        if (e == VPL.get(i+1,k) + WIP.get(k+1,j-1) + ap_penalty + 2 * bp_penalty) {
            Trace_VPL(i+1,k,VPL.get(i+1,k));
            Trace_WIP(k+1,j-1,WIP.get(k+1,j-1));
            return;
        }
    }
    UNREACHABLE();
}
void pseudo_loop::Trace_VPL(cand_pos_t i, cand_pos_t j, energy_t e){
    if (debug) printf("VPL at %d and %d with %d\n", i, j, e);
    cand_pos_t min_Bp_j = std::min((cand_pos_tu)tree->b(i, j), (cand_pos_tu)tree->Bp(i, j)); //unsigned is a trick to make it take the min but exclude negatives
    for (cand_pos_t k = i + 1; k < min_Bp_j; ++k) {
        bool can_pair = tree->up[k - 1] >= (k - i);
        if(can_pair && e==static_cast<energy_t>((k-i)*cp_penalty) + VP.get(k,j)){
            Trace_VP(k,j,VP.get(k,j));
            return;
        }
    }
    UNREACHABLE();
}
void pseudo_loop::Trace_VPR(cand_pos_t i, cand_pos_t j, energy_t e){
    if (debug) printf("VPR at %d and %d with %d\n", i, j, e);
    cand_pos_t max_i_bp = std::max(tree->B(i, j), tree->bp(i, j));

    for (cand_pos_t k = max_i_bp + 1; k < j; ++k) {
        if (e == VP.get(i,k) + WIP.get(k+1,j)) {
            Trace_VP(i,k,VP.get(i,k));
            Trace_WIP(k+1,j,WIP.get(k+1,j));
            return;
        }
    }
    for (cand_pos_t k = max_i_bp + 1; k < j; ++k) {
        bool can_pair = tree->up[j - 1] >= (j - k);
        if (can_pair && e==VP.get(i,k) + static_cast<energy_t>((j-k)*cp_penalty)){
            Trace_VP(i,k,VP.get(i,k));
            return;
        }
    }
    UNREACHABLE();
}
void pseudo_loop::Trace_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp, energy_t e){
    if (debug) printf("BE at %d and %d with %d\n", i, j, e);
    this->structure[i] = '(';
    this->structure[j] = ')';
    if(i==ip && j==jp) return; // This should ensure that the symbols are at every position
    if (tree->tree[i+1].pair == j-1 && e == (energy_t) (get_e_stP(i,j) + get_BE(i+1,j-1,ip,jp) + ShapeData->get_calculated(i) + ShapeData->get_calculated(j))) {
        Trace_BE(i+1,j-1,ip,jp,get_BE(i+1,j-1,ip,jp));
        return;
    }
    for (cand_pos_t l = i + 1; l <= ip; l++) {
        if (tree->tree[l].pair >= 0 && jp <= tree->tree[l].pair && tree->tree[l].pair < j) {
            cand_pos_t lp = tree->tree[l].pair;
            if(lp<0) continue; // No need to look at indices where there is no pair in BE

            bool empty_region_il = (tree->up[(l)-1] >= l - i - 1);       // empty between i+1 and lp-1
            bool empty_region_lpj = (tree->up[(j)-1] >= j - lp - 1);     // empty between l+1 and ip-1
            bool weakly_closed_il = tree->weakly_closed(i + 1, l - 1);   // weakly closed between i+1 and lp-1
            bool weakly_closed_lpj = tree->weakly_closed(lp + 1, j - 1); // weakly closed between l+1 and ip-1

            if (empty_region_il && empty_region_lpj) {
                if (e==get_e_intP(i,l,lp,j) + get_BE(l,lp,ip,jp)) {
                    Trace_BE(l,lp,ip,jp,get_BE(l,lp,ip,jp));
                    return;
                }
            }
            if (weakly_closed_il && weakly_closed_lpj) {
                if (e == WIP.get(i+1,l-1) + get_BE(l,lp,ip,jp) + WIP.get(lp+1,j-1) + ap_penalty + 2*bp_penalty) {
                    Trace_WIP(i+1,l-1,WIP.get(i+1,l-1));
                    Trace_BE(l,lp,ip,jp,get_BE(l,lp,ip,jp));
                    Trace_WIP(lp+1,j-1,WIP.get(lp+1,j-1));
                    return;
                }
            }
            if (weakly_closed_il && empty_region_lpj) {
                if (e == WIP.get(i+1,l-1) + get_BE(l,lp,ip,jp) + cp_penalty*(j-lp-1) + ap_penalty + 2*bp_penalty) {
                    Trace_WIP(i+1,l-1,WIP.get(i+1,l-1));
                    Trace_BE(l,lp,ip,jp,get_BE(l,lp,ip,jp));
                    return;
                }
            }
            if (empty_region_il && weakly_closed_lpj) {
                if (e == cp_penalty*(l-i-1) + get_BE(l,lp,ip,jp) + WIP.get(lp+1,j-1) + ap_penalty + 2*bp_penalty) {
                    Trace_BE(l,lp,ip,jp,get_BE(l,lp,ip,jp));
                    Trace_WIP(lp+1,j-1,WIP.get(lp+1,j-1));
                    return;
                }
            }
        }
    }
    UNREACHABLE();
}