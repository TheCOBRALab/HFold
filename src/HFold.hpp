#ifndef HFOLD_H
#define HFOLD_H
#include "sparse_tree.hpp"
#include "cmdline.hpp"
#include "Result.hpp"
#include "Hotspot.hpp"
#include <string>

struct RNAEntry {
    std::string name;
    std::string sequence;
    std::string structure;

    // Constructor that takes all three fields
    RNAEntry(std::string n, std::string s, std::string st)
        : name(std::move(n)), sequence(std::move(s)), structure(std::move(st)) {}

    // Default constructor (needed for vector resizing or default initialization)
    RNAEntry() = default;
};

std::string getSequence(args_info a);
void trim(std::string& s);
std::vector<RNAEntry> get_all_file_entries(const std::string& file);
std::vector<RNAEntry> get_all_inputs(const std::string& fileI, const std::string& seq, const std::string& restricted);
void validateStructure(std::string sequence, std::string structure);
void validateSequence(std::string sequence);
void seqtoRNA(std::string &seq);
std::string hfold(std::string seq, std::string res, double &energy, sparse_tree &tree, bool pk_free, bool pk_only, int dangles);
void handle_output_file(std::string& fileO);
void preprocess_sequence(std::string& seq, std::string& restricted, bool noConv);
void load_energy_parameters(const std::string& paramFile, const std::string& seq);
std::vector<Hotspot> build_hotspots(const std::string& seq, const std::string& restricted, int suboptCount);
std::vector<Result> fold_hotspots(
    const std::string& seq, 
    const std::vector<Hotspot>& hotspots,
    bool pk_free, bool pk_only, int dangles,
	bool input_structure_given
);
void output_results(
    const std::string& seq,
    const std::vector<Result>& results,
    const std::string& fileO,
    int suboptCount,
    const std::string& name = "",
    const std::size_t input_count = 1
);
#endif // HFOLD_H