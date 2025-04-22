// hfold_bindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <filesystem>

#include <HFold.hpp>
#include <Result.hpp>
#include <W_final.hpp>

namespace py = pybind11;


/*───────────────────────────────────────────────────────────
  Python Bindings (Functions that are exposed to Python)
───────────────────────────────────────────────────────────*/
namespace bindings
{
    std::vector<std::vector<Result>> hfold_simple(
        std::string sequence = "",
        std::string restricted = "",
        double energy = 0.0,
        bool pk_free = false,
        bool pk_only = false,
        int dangles = 2,
        int suboptCount = 1,
        std::string fileI = "",
        std::string fileO = "",
        std::string paramFile = "",
        bool noConv_given = false
    ) {
        bool input_structure_given = !restricted.empty();
        std::vector<RNAEntry> inputs = get_all_inputs(fileI, sequence, restricted);
        std::vector<std::vector<Result>> all_results;

        for (RNAEntry& current: inputs) {
            preprocess_sequence(current.sequence, current.structure, noConv_given);
            load_energy_parameters(paramFile, current.sequence);
            std::vector<Hotspot> hotspots = build_hotspots(current.sequence, current.structure, suboptCount);
            std::vector<Result> results = fold_hotspots(current.sequence, hotspots, pk_free, pk_only, dangles, input_structure_given);
            output_results(current.sequence, results, fileO, suboptCount, current.name, inputs.size());
            all_results.push_back(results);
        }

        return all_results;
    }
}


/*───────────────────────────────────────────────────────────
  Pybind11 module definition
───────────────────────────────────────────────────────────*/
PYBIND11_MODULE(hfold, m)
{
    m.doc() = "HFold minimal Python bindings";

    // Result class
    py::class_<Result>(m, "Result")
        .def(py::init<std::string, std::string, double, std::string, double>())
        .def("get_sequence", &Result::get_sequence)
        .def("get_restricted", &Result::get_restricted)
        .def("get_final_structure", &Result::get_final_structure)
        .def("get_restricted_energy", &Result::get_restricted_energy)
        .def("get_final_energy", &Result::get_final_energy)
        .def("__repr__", [](const Result& r) {
            return "<Result sequence='" + r.get_sequence() +
                   "', restricted='" + r.get_restricted() +
                   "', restricted_energy=" + std::to_string(r.get_restricted_energy()) +
                   ", final_structure='" + r.get_final_structure() +
                   "', final_energy=" + std::to_string(r.get_final_energy()) + ">";
        });

    // hfold function
    m.def(
        "hfold",
        &bindings::hfold_simple,
        py::arg("sequence") = "",
        py::arg("structure") = "",
        py::arg("energy") = 0.0,
        py::arg("pk_free") = false,
        py::arg("pk_only") = false,
        py::arg("dangles") = 2,
        py::arg("suboptCount") = 1,
        py::arg("fileI") = "",
        py::arg("fileO") = "",
        py::arg("paramFile") = "",
        py::arg("noConv_given") = false,
        R"pbdoc(
              Fold an RNA/DNA sequence.

              Parameters
              ----------
              sequence : str, optional
              structure : str, optional
              energy : float, optional
              pk_free : bool, optional
              pk_only : bool, optional
              dangles : int, optional
              suboptCount : int, optional
              fileI : str, optional
              fileO : str, optional
              paramFile : str, optional
              noConv_given : bool, optional

              Returns
              -------
              list[list[Result]]
        )pbdoc");
}
