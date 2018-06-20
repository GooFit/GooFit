#pragma once

#include <memory>
#include <vector>

#include <goofit/Variable.h>

namespace GooFit {

class PdfBase;
class DataSet;

class Workspace {
  private:
    std::vector<std::unique_ptr<PdfBase>> pdfs_;
    std::vector<Variable> vars_;
    std::vector<Observable> obs_;
    std::vector<std::unique_ptr<DataSet>> datasets_;

  public:
    Workspace()             = default;
    Workspace(Workspace &)  = delete;
    Workspace(Workspace &&) = default;

    void add_pdf(PdfBase *);
    void add_var(Variable);
    void add_obs(Observable);
    void add_obs(EventNumber);
    void add_dataset(DataSet *);

    // add "import" - it should make copies, so ensure that works first
    // Does writeToFile make sense?
    // Adding factory will be very challenging (though rewarding)
    // Could a workspace be exported to Python somehow?

    PdfBase *pdf(std::string name);
    Variable *var(std::string name);
    Observable *obs(std::string name);
    DataSet *dataset(int i);

    /// Print the contents of a Workspace to a string
    std::string Print() const;
};

} // namespace GooFit
