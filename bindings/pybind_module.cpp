// hfold_bindings.cpp
//  python setup.py build_ext --inplace

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>           // std::pair ⇆ Python tuple
#include "HFold.hpp"                 // real core routine
#include "sparse_tree.hpp"           // whatever defines sparse_tree

namespace py = pybind11;

/*───────────────────────────────────────────────────────────
  Helper that hides every C++‑level knob from Python
───────────────────────────────────────────────────────────*/
namespace   // ← keeps the symbol local to this TU
{
    std::pair<std::string,double>
    hfold_simple(const std::string &sequence,
                 const std::string &restricted = "")
    {
        const int n = static_cast<int>(sequence.size());

        // If the user gave no restriction, make an all‑dots string
        std::string res_for_tree =
            restricted.empty() ? std::string(n, '.') : restricted;

        sparse_tree tree(res_for_tree, n);

        double energy = 0.0;
        std::string final_structure =
            hfold(sequence, res_for_tree, energy,
                  tree,
                  /*pk_free =*/false,
                  /*pk_only =*/false,
                  /*dangles =*/2);

        return {final_structure, energy};
    }
} // anonymous namespace


/*───────────────────────────────────────────────────────────
  Pybind11 module definition
───────────────────────────────────────────────────────────*/
PYBIND11_MODULE(hfold, m)
{
    m.doc() = "HFold minimal Python bindings";

    m.def("hfold",
          &hfold_simple,
          py::arg("sequence"),
          py::arg("structure") = "",
          R"pbdoc(
              Fold an RNA/DNA sequence.

              Parameters
              ----------
              sequence : str
                  Nucleotide sequence (A,C,G,U/T).
              structure : str, optional
                  Dot‑bracket restriction; must match sequence length.

              Returns
              -------
              tuple(str, float)
                  (final_structure, minimum_free_energy)
          )pbdoc");
}
