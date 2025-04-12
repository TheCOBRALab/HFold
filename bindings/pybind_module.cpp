#include <pybind11/pybind11.h>

namespace py = pybind11;

// A simple function to be exposed
int add(int i, int j) {
    return i + j;
}

// Create the module
PYBIND11_MODULE(hfold, m) {
    m.doc() = "Example pybind11 module"; // optional module docstring

    m.def("add", &add, "A function that adds two integers");
}
