#include <pybind11/pybind11.h>

#include <goofit/PDFs/basic/ArgusPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;
using namespace pybind11::literals;

struct HelpPrinter {
    std::string help_str;
    std::string getHelp() const { return help_str; }
    HelpPrinter(std::string input)
        : help_str(input) {}
};

void init_ArgusPdf(py::module &m) {
    py::class_<HelpPrinter>(m, "HelpPrinter")
        .def(py::init<std::string>())
        .def("_repr_markdown_", &HelpPrinter::getHelp)
        .def("__repr__", &HelpPrinter::getHelp);

    m.attr("ArgusHelp") = HelpPrinter(ArgusHelp);

    py::class_<ArgusPdf, GooPdf>(m, "ArgusPdf")
        .def(py::init<std::string, Observable, Variable, Variable, bool>(),
             "Standard Argus",
             "n"_a,
             "x"_a,
             "m"_a,
             "s"_a,
             "upper"_a)

        .def(py::init<std::string, Observable, Variable, Variable, bool, Variable>(),
             "Power version of Argus",
             "n"_a,
             "x"_a,
             "m"_a,
             "s"_a,
             "upper"_a,
             "power"_a)
        .def_static("help", []() { return HelpPrinter(ArgusHelp); });
}
