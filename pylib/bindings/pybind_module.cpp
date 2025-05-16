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
    std::vector<std::vector<py::dict>> hfold_simple(
        std::string sequence = "",
        std::string restricted = "",
        bool pk_free = false,
        bool pk_only = false,
        int dangles = 2,
        int suboptCount = 1,
        std::string fileI = "",
        std::string fileO = "",
        std::string param_file = "",
        bool noConv_given = false,
        bool suppress_output = false
    ) {
        bool input_structure_given = !restricted.empty();
        std::vector<RNAEntry> inputs = get_all_inputs(fileI, sequence, restricted);
        std::vector<std::vector<py::dict>> all_results;

        for (RNAEntry& current: inputs) {
            preprocess_sequence(current.sequence, current.structure, noConv_given);
            load_energy_parameters(param_file, current.sequence);
            std::vector<Hotspot> hotspots = build_hotspots(current.sequence, current.structure, suboptCount);
            std::vector<Result> results = fold_hotspots(current.sequence, hotspots, pk_free, pk_only, dangles, input_structure_given);
            if (!suppress_output) {
                output_results(current.sequence, results, fileO, suboptCount, current.name, inputs.size());
            }

            std::vector<py::dict> result_dicts;
            for (const Result& r : results) {
                py::dict d;
                d["sequence"] = r.get_sequence();
                d["restricted"] = r.get_restricted();
                d["restricted_energy"] = r.get_restricted_energy();
                d["final_structure"] = r.get_final_structure();
                d["final_energy"] = r.get_final_energy();
                result_dicts.push_back(d);
            }

            all_results.push_back(result_dicts);
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
        py::arg("restriction") = "",
        py::arg("pk_free") = false,
        py::arg("pk_only") = false,
        py::arg("dangles") = 2,
        py::arg("subopt") = 1,
        py::arg("fileI") = "",
        py::arg("fileO") = "",
        py::arg("param_file") = "",
        py::arg("noConv_given") = false,
        py::arg("suppress_output") = false,
        R"pbdoc(
              Fold an RNA/DNA sequence.

              Parameters
              ----------
              sequence : str, optional
              restriction : str, optional
              pk_free : bool, optional
              pk_only : bool, optional
              dangles : int, optional
              subopt : int, optional
              fileI : str, optional
              fileO : str, optional
              param_file : str, optional
              noConv_given : bool, optional
              suppress_output : bool, optional

              Returns
              -------
              list[list[Result]]
        )pbdoc");
}
