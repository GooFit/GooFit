#include <pybind11/pybind11.h>

#include <pybind11/stl.h>

#include <Minuit2/FCNBase.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ROOT::Minuit2;

class PyFCNBase : public FCNBase {
  public:
    using FCNBase::FCNBase;

    double operator()(const std::vector<double> &v) const override {
        PYBIND11_OVERLOAD_PURE_NAME(double, FCNBase, "__call__", operator(), v);
    }

    double Up() const override { PYBIND11_OVERLOAD_PURE(double, FCNBase, Up, ); }

    void SetErrorDef(double v) override { PYBIND11_OVERLOAD(void, FCNBase, SetErrorDef, v); }

    double ErrorDef() const override { PYBIND11_OVERLOAD(double, FCNBase, ErrorDef, ); }
};

void init_FCNBase(py::module &m) {
    py::class_<FCNBase, PyFCNBase>(m, "FCNBase")
        .def(py::init<>())
        .def("__call__", &FCNBase::operator())
        .def("SetErrorDef", &FCNBase::SetErrorDef)
        .def("ErrorDef", &FCNBase::ErrorDef)
        .def("Up", &FCNBase::Up);
}
