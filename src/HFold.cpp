// Iterative HFold files
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
#include <filesystem>

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
    if(!std::filesystem::exists(file)){
        std::cerr << "Error: Input file not found: " << file << std::endl;
        std::exit(EXIT_FAILURE);
    }
    std::ifstream in(file.c_str());
    std::string line;
    RNAEntry current;
    std::vector<RNAEntry> entries;
    int state = 0;
    while(getline(in, line)){
        trim(line);
        if (line.empty()) continue;
        if (line[0] == '>'){
            if (!current.sequence.empty() && !current.structure.empty()) {
                entries.push_back(current);
                current = RNAEntry();
            }
            current.name = line.substr(1);
            state = 1;
        } else if (state == 1){
            current.sequence = line;
            state = 2;
        } else if (state == 2) {
            current.structure = line;
            state = 0;
        }
    }

    if (!current.sequence.empty() && !current.structure.empty()) {
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

void handle_output_file(std::string& fileO) {
    /**
     * @brief Handles the output file path.
     * 
     * This function checks if the output file already exists and prompts the user
     * for confirmation to overwrite it. If the user chooses to cancel, the program exits.
     * If the user provides an empty string, it will keep prompting until a valid path is given.
     * 
     * @param fileO The output file path to be checked and possibly modified.
     * @return void
     * 
     * @deprecated (Keeping this function in case we decide to overwrite the output file instead of append)
     */

    if (fileO.empty()) return; // User didn't provide an output file
    
    while (fileO.empty() || std::filesystem::exists(fileO)){
        if (fileO.empty()){
            std::cout << "Enter output file path (or 'c' to cancel): ";
            std::getline(std::cin, fileO);
            if (fileO == "c" or fileO == "C") std::exit(0);
        } else {
            std::cout << "Output file already exists. Overwrite? (y/n/c): ";
            char choice;
            std::cin >> choice;
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // clear the input buffer
            choice = std::tolower(choice);
            
            if (choice == 'y') return;
            if (choice == 'c') std::exit(0);
            if (choice == 'n') fileO.clear();
        }
    }
}

void preprocess_sequence(std::string& seq, std::string& restricted, bool noConv) {
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper); // convert sequence to uppercase
    if (!noConv) seqtoRNA(seq);
    validateSequence(seq);
    if (!restricted.empty()) validateStructure(seq, restricted);
}

void load_energy_parameters(const std::string& paramFile, const std::string& seq) {
    if (std::filesystem::exists(paramFile)) {
        vrna_params_load(paramFile.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
    } else if (seq.find('T') != std::string::npos) { // if T is present, load DNA parameters
        vrna_params_load_DNA_Mathews2004();
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
    const std::string& name = "",
    const std::size_t input_count = 1
) {
    int number_of_output = std::min(results.size(), static_cast<std::size_t>(suboptCount));

    if (!fileO.empty()) { //output to file
        std::ofstream out(fileO, std::ios::app);
        if (!name.empty()){
            out << ">" << name << std::endl; 
        }
        out << seq << std::endl;
        for (int i = 0; i < number_of_output; i++) {
            out << "Restricted_" << i << ": " << results[i].get_restricted() << std::endl;
            out << "Result_" << i << ":     " << results[i].get_final_structure()
                << " (" << results[i].get_final_energy() << ")" << std::endl;
        }
    } else { // output to console
        if (!name.empty() && input_count > 1){
            std::cout << ">" << name << std::endl; 
        }
        std::cout << seq << std::endl;
        if (results.size() == 1) { 
			// std::cout << "Restricted_" << 0 << ": " << results[0].get_restricted() << std::endl;
            std::cout << results[0].get_final_structure()
                      << " (" << results[0].get_final_energy() << ")" << std::endl;
        } else {
            for (int i = 0; i < number_of_output; i++) {
                if (i && results[i].get_final_structure() == results[i - 1].get_final_structure()) // skip duplicates
                    continue;
                std::cout << "Restricted_" << i << ": " << results[i].get_restricted() << std::endl;
                std::cout << "Result_" << i << ":     " << results[i].get_final_structure()
                          << " (" << results[i].get_final_energy() << ")" << std::endl;
            }
        }
    }
}