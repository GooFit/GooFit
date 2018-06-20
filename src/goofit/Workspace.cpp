#include <goofit/Workspace.h>

#include <goofit/DataSet.h>
#include <goofit/PdfBase.h>

#include <algorithm>
#include <sstream>
#include <string>

namespace GooFit {

void Workspace::add_pdf(PdfBase *item) {
    if(pdf(item->getName()) == nullptr)
        pdfs_.emplace_back(item);
}

void Workspace::add_var(Variable item) {
    if(pdf(item.getName()) == nullptr)
        vars_.push_back(item);
}

void Workspace::add_obs(Observable item) {
    if(pdf(item.getName()) == nullptr)
        obs_.push_back(item);
}

void Workspace::add_obs(EventNumber item) {
    if(pdf(item.getName()) == nullptr)
        obs_.emplace_back(item);
}

void Workspace::add_dataset(DataSet *item) { datasets_.emplace_back(item); }

PdfBase *Workspace::pdf(std::string name) {
    auto item = std::find_if(
        pdfs_.begin(), pdfs_.end(), [name](const std::unique_ptr<PdfBase> &item) { return item->getName() == name; });
    if(item == pdfs_.end())
        return nullptr;
    else
        return item->get();
}

Variable *Workspace::var(std::string name) {
    auto item
        = std::find_if(vars_.begin(), vars_.end(), [name](const Variable &item) { return item.getName() == name; });
    if(item == vars_.end())
        return nullptr;
    else
        return &*item;
}

Observable *Workspace::obs(std::string name) {
    auto item
        = std::find_if(obs_.begin(), obs_.end(), [name](const Observable &item) { return item.getName() == name; });
    if(item == obs_.end())
        return nullptr;
    else
        return &*item;
}

DataSet *Workspace::dataset(int i) { return datasets_.at(i).get(); }

std::string Workspace::Print() const {
    std::stringstream out;

    if(!vars_.empty()) {
        out << "Variables\n"
               "---------\n";
        for(const Variable &item : vars_)
            out << item << "\n";

        out << "\n";
    }

    if(!obs_.empty()) {
        out << "Observables\n"
               "-----------\n";
        for(const Observable &item : obs_)
            out << item << "\n";

        out << "\n";
    }

    if(!pdfs_.empty()) {
        out << "PDFs\n"
               "----\n";
        for(const auto &item : pdfs_)
            out << *item << "\n";

        out << "\n";
    }

    if(!datasets_.empty()) {
        out << "DataSets\n"
               "--------\n";
        for(const auto &item : datasets_)
            out << *item << "\n";

        out << "\n";
    }

    return out.str();
}

} // namespace GooFit
