// hfold_bindings.cpp
//  python3 setup.py build_ext --inplace

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>            // std::pair ⇆ Python tuple
#include <HFold.hpp>
#include <sparse_tree.hpp>
#include <W_final.hpp>
#include <filesystem>

namespace py = pybind11;


/*───────────────────────────────────────────────────────────
  Internal helpers for Python bindings (Functions that are not exposed to Python)
───────────────────────────────────────────────────────────*/
namespace bindings::helpers
{
    vrna_param_s* load_parameters(std::string param_file){
        if (std::filesystem::exists(param_file)) {
            vrna_params_load(param_file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
        } else {
            throw std::runtime_error("Parameter file not found: " + param_file);
        }
    
        // Scale parameters (matches main behavior)
        vrna_param_s* params = scale_parameters();
        return params;
    }
}


/*-───────────────────────────────────────────────────────────
Python Bindings (Functions that are exposed to Python)
───────────────────────────────────────────────────────────*/
namespace bindings
{   
    std::pair<std::string, double>
    hfold_simple(std::string sequence,
                 std::string constraint = "",
                 const std::string &param_file = "")
    {   
        /** 
        * @brief Predict a single RNA secondary structure using HFold with optional constraints.
        * 
        * This is a simplified wrapper around the HFold core algorithm, designed for use from Python.
        * It takes an RNA or DNA sequence, optionally applies a structure constraint (in dot-bracket notation),
        * and returns the predicted minimum-free-energy structure and its associated energy.
        * 
        * @param sequence The nucleotide sequence to fold. Can include A, C, G, U (or T, which will be converted to U).
        * @param constraint Optional dot-bracket structure constraint. Must be the same length as the sequence. 
        *                   If omitted, all positions are considered unpaired.
        * @param param_file Optional path to a ViennaRNA parameter file.
        */
        validateSequence(sequence);
        seqtoRNA(sequence);
        
        if (!constraint.empty()) {
            validateStructure(sequence, constraint);
        } else {
            // If the user gave no restriction, make an all‑dots string
            constraint = std::string(sequence.length(), '.');
        }

        if (!param_file.empty()) {
            helpers::load_parameters(param_file);
        }

        // get size of the sequence to allocate the sparse tree
        const int n = static_cast<int>(sequence.size());
        sparse_tree tree(constraint, n);

        double energy = 0.0;
        bool pk_free = false;
        bool pk_only = false;
        int dangles = 2;

        std::string final_structure =
            hfold(sequence,
                  constraint,
                  energy,
                  tree,
                  pk_free,
                  pk_only,
                  dangles);
        
        // If it takes energy to fold, it's best to not fold
        if (energy > 0.0) {
            energy = 0.0;
            constraint = std::string(sequence.length(), '.');
        }

        return {final_structure, energy};
    }
} 


/*───────────────────────────────────────────────────────────
  Pybind11 module definition
───────────────────────────────────────────────────────────*/
PYBIND11_MODULE(hfold, m)
{
    m.doc() = "HFold minimal Python bindings";
    m.def("hfold",
          &bindings::hfold_simple,
          py::arg("sequence"),
          py::arg("structure") = "",
          py::arg("param_file") = "",
          R"pbdoc(
              Fold an RNA/DNA sequence.

              Parameters
              ----------
              sequence : str
                  Nucleotide sequence (A,C,G,U/T).
              structure : str, optional
                  Dot-bracket restriction; must match sequence length.
              param_file : str, optional
                  Path to parameter file

              Returns
              -------
              tuple(str, float)
                  (final_structure, minimum_free_energy)
          )pbdoc");
}
