#ifndef SHAPE
#define SHAPE

#include "base_types.hh"
#include <vector>
#include <string>

// Consider floats instead of doubles
class SHAPEData {

    private:
        double  slope;         // m parameter for linear conversion to energy bonus
        double  intercept;     // b parameter
        cand_pos_t n;
        std::vector<double> calculated;  // pseudo-energy calculated SHAPE values, indexed by nucleotide position

        double calculate(double reactivity);
        bool exists(const std::string &filename);
    public:
        SHAPEData(const std::string &filename, cand_pos_t n, double slope = 1.8, double intercept = -.6);
        double get_calculated(cand_pos_t i);
};


#endif