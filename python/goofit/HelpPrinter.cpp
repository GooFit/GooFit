#include <goofit/Python.h>

void init_HelpPrinter(py::module &m) {
    py::class_<HelpPrinter>(m, "HelpPrinter")
        .def(py::init<std::string>())
        .def("_repr_markdown_", &HelpPrinter::getHelp)
        .def("__repr__", &HelpPrinter::getHelp);
}
