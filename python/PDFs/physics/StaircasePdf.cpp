#include <goofit/Python.h>

#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/physics/StaircasePdf.h>
#include <goofit/Variable.h>


using namespace GooFit;

void init_StaircasePdf(py::module &m) {
  py::class_<StaircasePdf, GooPdf>(m, "StaircasePdf")
    .def(py::init<std::string, Observable, std::vector<Variable*>>(),
	 "name"_a,
	 "_x"_a,
	 "x0list"_a);
}
