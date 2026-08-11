#ifndef PSEUDO_LOOP_H_
#define PSEUDO_LOOP_H_
#include "base_types.hh"
#include "h_struct.hh"
#include "SHAPE.hh"
#include "matrices.hh"
#include "sparse_tree.hh"
#include "Hotspot.hh"
#include "s_energy_matrix.hh"

#include "ViennaRNA/loops.hh"
#include "ViennaRNA/pair_mat.hh"
#include "ViennaRNA/params/io.hh"
#include <string>

#define debug 0

#ifdef NDEBUG
	#define UNREACHABLE() __builtin_unreachable()
#else
	#define UNREACHABLE() \
		do { \
			std::cerr << "Reached unreachable at line " << __LINE__ << " in File: " << __FILE__ << std::endl; \
			abort(); \
		} while(0)
#endif

void get_hotspots(std::string seq, std::vector<Hotspot> &hotspot_list, SHAPEData &ShapeData, int max_hotspot, vrna_param_t *params);
int distance(int left, int right);
void expand_hotspot(s_energy_matrix *V, Hotspot &hotspot, int n);
// Mateo 2024
// comparison function for hotspot so we can use it when sorting
bool compare_hotspot_ptr(Hotspot &a, Hotspot &b);
class pseudo_loop {

  public:
    // constructor
    pseudo_loop(std::string seq, std::string res,sparse_tree &tree, SHAPEData &ShapeData, bool pk_free, bool pk_only, int dangle);

    // destructor
    ~pseudo_loop();

    energy_t get_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp) {
        // Hosna, March 16, 2012,
        // i and j should be at least 3 bases apart
        if (j - i >= TURN && i >= 1 && i <= ip && ip < jp && jp <= j && j <= n && tree->tree[i].pair >= 0 && tree->tree[j].pair >= 0
            && tree->tree[ip].pair >= 0 && tree->tree[jp].pair >= 0 && tree->tree[i].pair == j && tree->tree[j].pair == i && tree->tree[ip].pair == jp
            && tree->tree[jp].pair == ip) {
            if (i == ip && j == jp && i < j) {
                return 0;
            }
            cand_pos_t iip = index[i] + ip - i;

            return BE[iip];
        } else {
            return INF;
        }
    }
    std::string structure;
    double hfold();
  private:
    cand_pos_t n;
    std::string res;
    std::string seq_;
    sparse_tree* tree;

    vrna_param_t *params_;
    SHAPEData *ShapeData;
    bool pk_free = false;
    bool pk_only = false;

    std::vector<energy_t> W;
    std::vector<free_energy_node> V;
    TriangleMatrix WM;
    TriangleMatrix WMv;
    TriangleMatrix WMp;

    TriangleMatrix WI;   // the loop inside a pseudoknot (in general it looks like a W but is inside a pseudoknot)
    TriangleMatrix VP;   // the loop corresponding to the pseudoknotted region of WMB
    TriangleMatrix VPL;  // the loop corresponding to the pseudoknotted region of WMB
    TriangleMatrix VPR;  // the loop corresponding to the pseudoknotted region of WMB
    TriangleMatrix WMB; // the main loop for pseudoloops and bands
    TriangleMatrix WMBP; // the main loop to calculate WMB
    TriangleMatrix WMBW;
    TriangleMatrix WIP;     // the loop corresponding to WI'
    TriangleMatrix BE;      // the loop corresponding to BE
    std::vector<cand_pos_t> index; // the array to keep the index of two dimensional arrays like WI and weakly_closed

    short *S_;
    short *S1_;

    void compute_energy_restricted(cand_pos_t i, cand_pos_t j);
    energy_t compute_internal_restricted(cand_pos_t i, cand_pos_t j);
    void compute_energy_WM_restricted(cand_pos_t i, cand_pos_t j);
    energy_t compute_energy_VM_restricted(cand_pos_t i, cand_pos_t j);
    void compute_WMv_WMp(cand_pos_t i, cand_pos_t j);

    void compute_energies_PK(cand_pos_t i, cand_pos_t j);
    void compute_WMB(cand_pos_t i, cand_pos_t j);
    void compute_WI(cand_pos_t i, cand_pos_t j);
    void compute_VP(cand_pos_t i, cand_pos_t j);
    void compute_VPL(cand_pos_t i, cand_pos_t j);
    void compute_VPR(cand_pos_t i, cand_pos_t j);
    void compute_WMBP(cand_pos_t i, cand_pos_t j);
    void compute_WMBW(cand_pos_t i, cand_pos_t j);
    void compute_WIP(cand_pos_t i, cand_pos_t j);
    void compute_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp);

    // Traceback //
	void backtrack();
	void Trace_W(cand_pos_t i, cand_pos_t j, energy_t e);
	void Trace_V(cand_pos_t i, cand_pos_t j, energy_t e);
	void Trace_WM(cand_pos_t i, cand_pos_t j, energy_t e);
	void Trace_WMv(cand_pos_t i, cand_pos_t j, energy_t e);
	void Trace_WMp(cand_pos_t i, cand_pos_t j, energy_t e);

    void Trace_WI(cand_pos_t i, cand_pos_t j, energy_t e);
    void Trace_WIP(cand_pos_t i, cand_pos_t j, energy_t e);
    void Trace_WMB(cand_pos_t i, cand_pos_t j, energy_t e);
    void Trace_WMBP(cand_pos_t i, cand_pos_t j, energy_t e);
    void Trace_WMBW(cand_pos_t i, cand_pos_t j, energy_t e);
    void Trace_VP(cand_pos_t i, cand_pos_t j, energy_t e);
    void Trace_VPL(cand_pos_t i, cand_pos_t j, energy_t e);
    void Trace_VPR(cand_pos_t i, cand_pos_t j, energy_t e);
    void Trace_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp, energy_t e);
    



    // Util
    void allocate_space();
    energy_t get_energy (cand_pos_t i, cand_pos_t j) { if (i>=j) return INF; cand_pos_t ij = index[i]+j-i; return V[ij].energy;}
    char get_type (cand_pos_t i, cand_pos_t j) {cand_pos_t ij = index[i]+j-i; return V[ij].type;}
    int compute_exterior_cases(cand_pos_t l, cand_pos_t j);
    energy_t HairpinE(const std::string &seq, cand_pos_t i, cand_pos_t j);
    energy_t compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l);
    energy_t get_e_stP(cand_pos_t i, cand_pos_t j);
    energy_t get_e_intP(cand_pos_t i, cand_pos_t ip, cand_pos_t jp, cand_pos_t j);
    energy_t E_ext_Stem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1, const cand_pos_t i,const cand_pos_t j);
    energy_t E_MLStem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,cand_pos_t i, cand_pos_t j);
    energy_t E_MbLoop(const energy_t WM2ij, const energy_t WM2ip1j, const energy_t WM2ijm1, const energy_t WM2ip1jm1, cand_pos_t i, cand_pos_t j);
};
#endif /*PSEUDO_LOOP_H_*/
