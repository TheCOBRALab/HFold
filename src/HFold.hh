#ifndef HFOLD_H
#define HFOLD_H
#include "sparse_tree.hh"
#include <string>

bool exists (const std::string path);
void get_input(std::string file, std::string &sequence, std::string &structure);
void validateStructure(std::string sequence, std::string structure);
void validateSequence(std::string sequence);
void seqtoRNA(std::string &seq);
std::string hfold(std::string seq, std::string structure, double &energy, sparse_tree &tree, bool pk_free, bool pk_only, int dangles);

#endif // HFOLD_H