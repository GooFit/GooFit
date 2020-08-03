#include <goofit/PDFs/GooPdf.h> 
#include <goofit/PdfBase.h>
#include <goofit/Variable.h>

#include <vector>
namespace GooFit {

class StaircasePdf : public GooPdf {
 public:
  StaircasePdf (std::string n, Observable _x, std::vector<Variable> x0list);

 private:
};
}
