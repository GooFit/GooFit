#include "goofit/FunctorWriter.h"
#include "goofit/PdfBase.h"
#include "goofit/Variable.h"
#include <fstream>
#include <map>

namespace GooFit {

void writeToFile(PdfBase *pdf, const char *fname) {
    std::vector<Variable *> params = pdf->getParameters();

    std::ofstream writer;
    writer.open(fname);

    for(Variable *p : params) {
        writer << p->getName() << " " << p->getValue() << " " << p->getError() << " " << p->getNumBins() << " "
               << p->getLowerLimit() << " " << p->getUpperLimit() << std::endl;
    }

    writer.close();
}

void readFromFile(PdfBase *pdf, const char *fname) {
    std::vector<Variable *> params = pdf->getParameters();

    std::map<std::string, Variable *> tempMap;

    for(Variable *p : params) {
        tempMap[p->getName()] = p;
    }

    std::ifstream reader;
    reader.open(fname);
    std::string buffer;
    char discard[1000];
    int numSet = 0;

    while(true) {
        reader >> buffer;

        if(reader.eof())
            break;

        Variable *var = tempMap[buffer];

        fptype value, error, lowerlimit, upperlimit;
        size_t numbins;

        if(var) {
            reader >> value >> error >> numbins >> lowerlimit >> upperlimit;

            var->setValue(value);
            var->setError(error);
            var->setNumBins(numbins);
            var->setLowerLimit(lowerlimit);
            var->setUpperLimit(upperlimit);

            if(++numSet == tempMap.size())
                break;
        } else {
            reader.getline(discard, 1000);
        }
    }

    reader.close();
}

void readListOfNumbers(thrust::host_vector<fptype> &target, const char *fname) {
    std::ifstream reader;
    reader.open(fname);
    fptype buffer = 0;

    while(true) {
        reader >> buffer;

        if(reader.eof())
            break;

        target.push_back(buffer);
    }

    reader.close();
}

void writeListOfNumbers(thrust::host_vector<fptype> &target, const char *fname) {
    std::ofstream writer;
    writer.open(fname);

    for(fptype t : target) {
        writer << t << " ";
    }

    writer.close();
}
} // namespace GooFit
