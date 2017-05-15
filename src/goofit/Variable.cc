#include "goofit/Variable.h"
#include "goofit/Error.h"

#include <algorithm>

void Variable::setNumBins(size_t num) {
    if(locked_)
        throw GooFit::GeneralError("BinnedDataSet does not allow the number of bins to be changed after creation");
    numbins = num;
}

std::ostream& operator<< (std::ostream& o, const Variable& var) {
    o << var.getName() << ": " << var.getValue() << " +/- " << var.getError();
    if(!var.fixed)
        o << " [" << var.getLowerLimit() << ", " << var.getUpperLimit() << "]";
    if(var.getIndex() >= 0)
        o << " GooFit index: " << var.getIndex();
    if(var.getFitterIndex() >= 0)
        o << " Fitter index: " << var.getFitterIndex();
    if(var.blind != 0)
        o << " Blinded";
    
    return o;
}

std::istream& operator>> (std::istream& i, Variable& var) {
    return i >> var.value;
}

/// Get the max index of a variable from a list
int max_index(const std::vector<Variable*> &vars) {
    const Variable* max_ind_ptr = *std::max_element(std::begin(vars),
                                                    std::end(vars),
                                                    [](const Variable *a, const Variable *b)
                                                    {return a->getIndex() < b->getIndex();});
    return max_ind_ptr->getIndex();
}
