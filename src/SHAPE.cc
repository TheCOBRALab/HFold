#include "SHAPE.hh"
#include "ViennaRNA/utils.hh"
#include <sys/stat.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <stdexcept>


SHAPEData::SHAPEData(const std::string &filename, cand_pos_t n, double slope, double intercept): slope(slope), intercept(intercept), n(n), calculated(n + 1, 0.0){
    if (!exists(filename)){
        if(filename != ""){
            vrna_message_warning("Warning, invalid Shape file name: %s\n", filename.c_str());
        }
        return;
    }
    std::ifstream in(filename);
    std::string line;
    if (!std::getline(in, line)) return;

    std::istringstream ss(line);
    auto nextToken = [&](std::string &tok) {
        return static_cast<bool>(ss >> tok);
    };
    std::string tok;
    if (!nextToken(tok)) return;
    // Determine if the first field is a name or the length.
    // If it parses as an integer, there is no name field.
    if (!std::isdigit(static_cast<unsigned char>(tok.front()))) {
        // First token was the name; advance to the length field
        if (!nextToken(tok)) return;
    }
    cand_pos_t length = static_cast<cand_pos_t>(std::stoi(tok));
    if (length > n) {
        std::cerr << "ERROR: SHAPE file length (" << length << ") exceeds sequence length (" << n << ")\n";
        std::exit(EXIT_FAILURE);
    }
    // Skip the energy field
    if (!nextToken(tok)) return;
    // Read the reactivity portion
    cand_pos_t i = 0;
    while (nextToken(tok) && i < n) {
        ++i;
        if (tok == "NULL") continue;
        double reactivity = std::stod(tok);
        calculated[i] = 100.0*calculate(reactivity);
    }
    in.close();
}

bool SHAPEData::exists(const std::string &filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

double SHAPEData::get_calculated(cand_pos_t index)
{
    // An empty filename creates a zero-length SHAPE data set. This is the
    // normal "SHAPE disabled" state, so every contribution is zero.
    if (n == 0) {
        return 0.0;
    }

    if (index < 0 || index > n) {
        throw std::out_of_range(
            "SHAPEData index " + std::to_string(index) +
            " is outside [0, " + std::to_string(n) + "]"
        );
    }

    return calculated.at(static_cast<std::size_t>(index));
}

double SHAPEData::calculate(double reactivity){
    return slope*std::log(reactivity+1) + intercept;
}
