// Iterative HFold files
#include "HFold.hpp"
#include "cmdline.hpp"
#include "W_final.hpp"
#include "h_globals.hpp"
#include "cmdline.hpp"
#include "Result.hpp"
#include "Hotspot.hpp"
// a simple driver for the HFold
#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>

// filesystem::exists not supported in older macOS which is needed for the Conda Build
bool file_exists(const std::string& name) { 
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

std::string getSequence(args_info a) {
	if (a.inputs_num > 0) {
		return a.inputs[0];
	} else if (!a.input_file_given) {
		std::cout << "Sequence: ";
		std::string seq;
		std::getline(std::cin, seq);
		return seq;
	}
    return "";
}

void trim(std::string& s) {
    /**
     * @brief Trims leading and trailing whitespace from a string.
     * 
     * This function modifies the input string in-place to remove any leading
     * and trailing whitespace characters, including spaces, tabs, newlines,
     * carriage returns, form feeds, and vertical tabs.
     * 
     * @param s The string to be trimmed.
     * @return void
     */
    s.erase(0, s.find_first_not_of(" \t\n\r\f\v"));
    s.erase(s.find_last_not_of(" \t\n\r\f\v") + 1);
}

std::vector<RNAEntry> get_all_file_entries(const std::string& file){
    if(!file_exists(file)){
        std::cerr << "Error: Input file not found: " << file << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // state machine to parse the file
    #define UNINITIALIZED -1
    #define NAME 0
    #define SEQUENCE 1
    #define STRUCTURE 2

    std::ifstream in(file.c_str());
    std::string line;
    RNAEntry current;
    std::vector<RNAEntry> entries;
    int state = UNINITIALIZED;
    int line_number = 0;

    while(getline(in, line)){
        ++line_number;
        trim(line);
        if (line.empty()) continue;

        // Check if the line is the name of the entry
        if ((state == NAME || state == UNINITIALIZED) && (line[0] != '>')) {
            std::cerr << "Error: Expected '>' at the beginning of the line: " << line << ". Line number: " << line_number <<std::endl;
            exit(EXIT_FAILURE);
        }

        if (line[0] == '>'){

            if (!current.name.empty() && !current.sequence.empty() && current.structure.empty()) {
                current.structure = "";
            }

            if (state != UNINITIALIZED) { // valid entry, save it
                if (current.sequence.empty() && current.structure.empty()) {
                    std::cerr << "Warning: Sequence and structure are empty for entry: " << current.name << ". Line number: " << line_number << ". Skipping..."<<  std::endl;
                }
                entries.push_back(current);
                current = RNAEntry();
            }

            current.name = line.substr(1);
            state = SEQUENCE;

        } else if (state == SEQUENCE){
            if (!validateSequence(line, false)) {
                std::cerr << "Error: Sequence is invalid for entry: " << current.name  << ". Line number: " << line_number << std::endl;
                exit(EXIT_FAILURE);
            }
            current.sequence = line;
            state = STRUCTURE;
            
        } else if (state == STRUCTURE) {
            if (!validateStructure(current.sequence, line, false)) {
                std::cerr << "Error: Structure is invalid for entry: " << current.name << ". Line number: " << line_number  << std::endl;
                exit(EXIT_FAILURE);
            }
            current.structure = line;
            state = NAME;
        } else {
            // Should never reach here
            std::cerr << "Error: Unexpected state at line " << line_number << ": " << line << ". Line number: " << line_number  << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    // Handle the last entry
    if (!current.name.empty() && !current.sequence.empty()) {
        if (current.structure.empty()) {
            current.structure = "";
        }
        entries.push_back(current);
    }

    return entries;
}

std::vector<RNAEntry> get_all_inputs(const std::string& fileI, const std::string& seq, const std::string& restricted) {
    std::vector<RNAEntry> entries;
    if (!seq.empty()) {
        entries.emplace_back("Console Sequence", seq, restricted);
    }
    if (!fileI.empty()){
        std::vector<RNAEntry> file_entries = get_all_file_entries(fileI);
        entries.insert(entries.end(), file_entries.begin(), file_entries.end());
    }
    if (entries.empty()) throw std::runtime_error("Sequence is missing");
    return entries;
}

//check length and if any characters other than ._()
bool validateStructure(std::string& sequence, std::string& structure, bool exit_on_invalid){
	if(structure.length() != sequence.length()){
        if (exit_on_invalid) {
            std::cerr << "Error: The length of the sequence and structure must be the same." << std::endl;
            exit(EXIT_FAILURE);
        }
        return false;
	}

	//check if any characters are not ._()
	for(char c : structure) {
		if (!(c == '.' || c == '_' || c == '(' || c == ')' || c == 'x')){
            if (exit_on_invalid) {
                std::cerr << "Error: Structure contains invalid character " << c << ". Allowed: ., _, (, ), x." << std::endl;
                exit(EXIT_FAILURE);
            }
            return false;
		}
	}
    return true;
}

//check if sequence is valid with regular expression
//check length and if any characters other than GCAUT
bool validateSequence(std::string& sequence, bool exit_on_invalid){ 
	if(sequence.length() == 0){
        if (exit_on_invalid) {
            std::cerr << "Error: Sequence is missing." << std::endl;
            exit(EXIT_FAILURE);
        }
        return false;
	}

  // return false if any characters other than GCAUT -- future implement check based on type
  for(char c : sequence) {
    if (!(c == 'G' || c == 'C' || c == 'A' || c == 'U' || c == 'T')) {
        if (exit_on_invalid) {
            std::cerr  << "Error: Sequence contains invalid character " << c << ". Allowed: G, C, A, U, T." << std::endl;
            exit(EXIT_FAILURE);
        }
        return false;
    }
  }
    return true;
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

void preprocess_sequence(std::string& seq, std::string& restricted, bool noConv) {
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper); // convert sequence to uppercase
    if (!noConv) seqtoRNA(seq);
    validateSequence(seq);
    if (!restricted.empty()) validateStructure(seq, restricted);
}

void load_energy_parameters(const std::string& paramFile, const std::string& seq) {
    if (paramFile.empty()) return; // No parameter file provided
    
    if (file_exists(paramFile)) {
        vrna_params_load(paramFile.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
    } else if (seq.find('T') != std::string::npos) { // if T is present, load DNA parameters
        vrna_params_load_DNA_Mathews2004();
    } else {
        std::cerr << "Error: Input file not found: " << paramFile << std::endl;
    }
}

std::vector<Hotspot> build_hotspots(const std::string& seq, const std::string& restricted, int suboptCount) {
    std::vector<Hotspot> hotspots;
    vrna_param_s* params = scale_parameters();

    if (!restricted.empty()) {
        hotspots.emplace_back(1, restricted.length(), restricted.length() + 1);
        hotspots.back().set_structure(restricted);
    }
    if (static_cast<int>(hotspots.size()) < suboptCount) {
        get_hotspots(seq, hotspots, suboptCount, params);
    }

    free(params);
    return hotspots;
}

std::vector<Result> fold_hotspots(
    const std::string& seq, 
    const std::vector<Hotspot>& hotspots,
    bool pk_free, bool pk_only, int dangles,
	bool input_structure_given
) {
    std::vector<Result> results;
    results.reserve(hotspots.size()); // Pre-allocate memory for results

    for (const Hotspot& hs : hotspots) {
        double energy = 0.0;
        sparse_tree tree(hs.get_structure(), seq.size());
        std::string final_structure = hfold(seq, hs.get_structure(), energy, tree, pk_free, pk_only, dangles);
		if (!input_structure_given && energy > 0.0) {
			energy = 0.0;
			final_structure = std::string(seq.length(), '.');
		}
        results.emplace_back(seq, hs.get_structure(), hs.get_energy(), final_structure, energy);
    }

    std::sort(results.begin(), results.end(), Result::Result_comp{});
    return results;
}

void output_results(
    const std::string& seq,
    const std::vector<Result>& results,
    const std::string& fileO,
    int suboptCount,
    const std::string& name,
    const std::size_t input_count,
    const bool skip_duplicates
) {
    int number_of_output = std::min(results.size(), static_cast<std::size_t>(suboptCount));

    if (!fileO.empty()) { //output to file
        std::ofstream out(fileO, std::ios::app);
        if (!name.empty()){
            out << ">" << name << std::endl; 
        }
        out << seq << std::endl;
        if (number_of_output == 1) {
            out << "Restricted" << ": " << results[0].get_restricted() << std::endl;
            out << "Result" << ":     " << results[0].get_final_structure()
                << " (" << results[0].get_final_energy() << ")" << std::endl;
        } else {
            for (int i = 0; i < number_of_output; i++) {
                if (skip_duplicates && i && results[i].get_final_structure() == results[i - 1].get_final_structure()){ // skip duplicates
                    continue;
                }
                out << "Restricted_" << i << ": " << results[i].get_restricted() << std::endl;
                out << "Result_" << i << ":     " << results[i].get_final_structure()
                    << " (" << results[i].get_final_energy() << ")" << std::endl;
            }
        }
    } else { // output to console
        if (!name.empty() && input_count > 1){
            std::cout << ">" << name << std::endl; 
        }
        
        if (results.size() == 1) { 
            std::cout << "Sequence:        " << seq << std::endl;
			std::cout << "Restricted:      " << results[0].get_restricted() << std::endl;
            std::cout << "Final Structure: " << results[0].get_final_structure()
                      << " (" << results[0].get_final_energy() << ")" << std::endl;
        } else {
            int alignment = std::floor(std::log10(number_of_output));
            std::cout << "Sequence:     " << std::string(alignment, ' ') << seq << std::endl;
            for (int i = 0; i < number_of_output; i++) {
                alignment = std::floor(std::log10(number_of_output)) - (std::floor(std::log10(i !=0 ? i : 1)));
                if (skip_duplicates && i && results[i].get_final_structure() == results[i - 1].get_final_structure()){ // skip duplicates
                    continue;
                }
                std::cout << "Restricted_" << i << ": " << std::string(alignment, ' ') << results[i].get_restricted() << std::endl;
                std::cout << "Result_" << i << ":     " << std::string(alignment, ' ') << results[i].get_final_structure()
                          << " (" << results[i].get_final_energy() << ")" << std::endl;
            }
        }
    }
}