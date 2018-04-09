#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <thrust/copy.h>

#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

using namespace pybind11::literals;

void init_DalitzPlotPdf(py::module &m) {
    py::class_<DalitzPlotPdf, GooPdf>(m, "DalitzPlotPdf")
        .def(py::init<std::string, Observable, Observable, EventNumber, DecayInfo3, GooPdf *>(),
             "n",
             "m12",
             "m13",
             "eventNumber",
             "decay",
             "eff",
             py::keep_alive<1, 6>(), // Important to keep decay alive, to keep PDFs alive
             py::keep_alive<1, 7>())
        .def("setDataSize", &DalitzPlotPdf::setDataSize, "dataSize"_a, "evtSize"_a = 3)
        .def("getCachedWave",
             [](DalitzPlotPdf &self, size_t i) {
                 auto ret_thrust = self.getCachedWave(i);
                 std::vector<std::complex<fptype>> ret(ret_thrust.size());
                 thrust::copy(ret_thrust.begin(), ret_thrust.end(), ret.begin());
                 return ret;
             },
             "i"_a)
        .def("sumCachedWave",
             [](DalitzPlotPdf &self, size_t i) { return std::complex<fptype>(self.sumCachedWave(i)); },
             "i"_a)
        .def("fit_fractions",
             &DalitzPlotPdf::fit_fractions,
             "Using the current dataset, return the cached fit fraction values");
}
