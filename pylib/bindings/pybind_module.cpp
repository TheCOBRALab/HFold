// hfold_bindings.cpp
//  python3 setup.py build_ext --inplace

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>            // std::pair ⇆ Python tuple
#include "../../src/HFold.hpp"
#include "../../src/sparse_tree.hpp"
#include "../../src/W_final.hpp"

namespace py = pybind11;


/*───────────────────────────────────────────────────────────
  Internal helpers for Python bindings (Functions that are not exposed to Python)
───────────────────────────────────────────────────────────*/
namespace bindings::helpers
{
    vrna_param_s* load_parameters(std::string param_file){
        if (exists(param_file)) {
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
                 const std::string &param_file = "./params/rna_DirksPierce09.par")
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
        * @param param_file Path to a ViennaRNA parameter file. Default is "./params/rna_DirksPierce09.par".
        */
        validateSequence(sequence);
        seqtoRNA(sequence);
        
        if (!constraint.empty()) {
            validateStructure(sequence, constraint);
        } else {
            // If the user gave no restriction, make an all‑dots string
            constraint = std::string(sequence.length(), '.');
        }

        helpers::load_parameters(param_file);

        // get size of the sequence
        const int n = static_cast<int>(sequence.size());

        sparse_tree tree(constraint, n);

        double energy = 0.0;
        std::string final_structure =
            hfold(sequence, constraint, energy,
                  tree,
                  /*pk_free =*/false,
                  /*pk_only =*/false,
                  /*dangles =*/2);

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
          py::arg("param_file") = "./params/rna_DirksPierce09.par",
          R"pbdoc(
              Fold an RNA/DNA sequence.

              Parameters
              ----------
              sequence : str
                  Nucleotide sequence (A,C,G,U/T).
              structure : str, optional
                  Dot‑bracket restriction; must match sequence length.
              param_file : str, optional
                  Path to parameter file (default: "../../params/rna_DirksPierce09.par").

              Returns
              -------
              tuple(str, float)
                  (final_structure, minimum_free_energy)
          )pbdoc");
}
