#include <pybind11/pybind11.h>

#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/PDFs/physics/TddpPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_DalitzVetoPdf(py::module &m) {
    /*

        py::class_<VetoInfo>(m, "VetoInfo")

                .def(py::init<const std::string&, Variable>())
                .def(py::init<const std::string&, Variable>())
                ;

     */
    py::class_<DalitzVetoPdf, GooPdf>(m, "DalitzVetoPdf");
}
