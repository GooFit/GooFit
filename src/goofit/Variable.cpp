#include <goofit/Error.h>
#include <goofit/Variable.h>
#include <goofit/utilities/Uncertain.h>

#include <algorithm>

namespace GooFit {

bool Indexable::isEqualNameValLimits(const Indexable &other) const {
    return this->name == other.name && this->value.get() == other.value.get()
           && this->lowerlimit.get() == other.lowerlimit.get() && this->upperlimit.get() == other.upperlimit.get();
}

bool Variable::isEqualNameValLimitsErrBlindFixed(const Variable &other) const {
    return this->isEqualNameValLimits(other) && this->error.get() == other.error.get()
           && this->blind.get() == other.blind.get() && this->fixed.get() == other.fixed.get();
}

auto max_fitter_index(const std::vector<Variable> &vars) -> int {
    const Variable max_ind_ptr
        = *std::max_element(std::begin(vars), std::end(vars), [](const Variable &a, const Variable &b) {
              return a.getFitterIndex() < b.getFitterIndex();
          });
    return max_ind_ptr.getFitterIndex();
}

auto operator<<(std::ostream &o, const GooFit::Variable &var) -> std::ostream & {
    o << var.getName() << ": ";

    o << Uncertain(var.getValue(), var.getError());

    if(!var.IsFixed())
        o << " [" << var.getLowerLimit() << ", " << var.getUpperLimit() << "]";
    if(var.getIndex() >= 0)
        o << " GooFit index: " << var.getIndex();
    if(var.getFitterIndex() >= 0)
        o << " Fitter index: " << var.getFitterIndex();
    if(*var.blind != 0)
        o << " Blinded";

    return o;
}

auto operator<<(std::ostream &o, const GooFit::Observable &var) -> std::ostream & {
    o << var.getName() << ": " << var.getValue() << " (" << var.getNumBins() << " bins)";
    o << " [" << var.getLowerLimit() << ", " << var.getUpperLimit() << "]";
    if(var.getIndex() >= 0)
        o << " GooFit index: " << var.getIndex();

    return o;
}

auto operator>>(std::istream &i, GooFit::Observable &var) -> std::istream & { return i >> *var.value; }

} // namespace GooFit
