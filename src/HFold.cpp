// Iterative HFold files
#include "hotspot.hpp"
#include "Result.hpp"
#include "cmdline.hpp"
#include "W_final.hpp"
#include "h_globals.hpp"
// a simple driver for the HFold
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>
#include <string>
#include <getopt.h>

bool exists (const std::string path) {
  struct stat buffer;   
  return (stat (path.c_str(), &buffer) == 0); 
}

void get_input(std::string file, std::string &sequence, std::string &structure){
	if(!exists(file)){
		std::cout << "Input file does not exist" << std::endl;
		exit(EXIT_FAILURE);
	}
	std::ifstream in(file.c_str());
	std::string str;
	int i = 0;
	while(getline(in,str)){
		if(str[0] == '>') continue;
		if(i==0) sequence = str;
		if(i==1) structure = str;
		++i;
	}
	in.close();
}

//check length and if any characters other than ._()
void validateStructure(std::string sequence, std::string structure){
	if(structure.length() != sequence.length()){
		std::cout << " The length of the sequence and corresponding structure must have the same length" << std::endl;
		exit(EXIT_FAILURE);
	}

	//check if any characters are not ._()
	for(char c : structure) {
		if (!(c == '.' || c == '_' || c == '(' || c == ')' || c == 'x')){
			std::cout << "Structure must only contain ._()x: " << c << std::endl;
			exit(EXIT_FAILURE);
		}
	}
}

//check if sequence is valid with regular expression
//check length and if any characters other than GCAUT
void validateSequence(std::string sequence){

	if(sequence.length() == 0){
		std::cout << "sequence is missing" << std::endl;
		exit(EXIT_FAILURE);
	}
  // return false if any characters other than GCAUT -- future implement check based on type
  for(char c : sequence) {
    if (!(c == 'G' || c == 'C' || c == 'A' || c == 'U' || c == 'T')) {
		std::cout << "Sequence contains character " << c << " that is not G,C,A,U, or T." << std::endl;
		exit(EXIT_FAILURE);
    }
  }
}

std::string hfold(std::string seq,std::string res, double &energy, sparse_tree &tree, bool pk_free, bool pk_only, int dangle){
	W_final min_fold(seq,res, pk_free, pk_only, dangle);
	energy = min_fold.hfold(tree);
    std::string structure = min_fold.structure;
    return structure;
}

void seqtoRNA(std::string &sequence){
	/**
	 * @brief Converts a DNA sequence to RNA by replacing T with U.
	 * 
	 * This function modifies the input string in-place and logs the 
	 * conversion if any T's are found.
	 * 
	 * @param sequence The DNA sequence to be converted (e.g. "ATCGT").
	 *                 The result (RNA) will be stored in the same string.
	 */
	bool isRNA = true;
	std::string original_sequence = sequence;
    for (char &c : sequence) {
      	if (c == 'T') {
			c = 'U';
			isRNA = false;
		}
    }
	if(!isRNA){
		std::cout << "Input sequence contains T's, converting to U's" << std::endl;
		std::cout << "Original sequence: " << original_sequence << std::endl;
		std::cout << "Converted sequence: " << sequence << std::endl;
	}
}
