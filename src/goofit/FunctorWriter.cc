#include "goofit/FunctorWriter.h"
#include <fstream>
#include <map>
#include "goofit/PdfBase.h"
#include "goofit/Variable.h"

void writeToFile(PdfBase* pdf, const char* fname) {
    std::vector<Variable*> params = pdf->getParameters();

    std::ofstream writer;
    writer.open(fname);

    for(Variable* p : params) {
        writer << p->GetName() << " "
               << p->GetValue() << " "
               << p->GetError() << " "
               << p->GetNumBins() << " "
               << p->GetLowerLimit() << " "
               << p->GetUpperLimit()
               << std::endl;
    }

    writer.close();
}


void readFromFile(PdfBase* pdf, const char* fname) {
    std::vector<Variable*> params = pdf->getParameters();

    std::map<std::string, Variable*> tempMap;

    for(Variable* p : params) {
        tempMap[p->name] = p;
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

        Variable* var = tempMap[buffer];
        
        fptype value, error, lowerlimit, upperlimit;
        size_t numbins;

        if(var) {
            reader >> value
                   >> error
                   >> numbins
                   >> lowerlimit
                   >> upperlimit;
            
            var->SetValue(value);
            var->SetError(error);
            var->SetNumBins(numbins);
            var->SetLowerLimit(lowerlimit);
            var->SetUpperLimit(upperlimit);

            if(++numSet == tempMap.size())
                break;
        } else {
            reader.getline(discard, 1000);
        }
    }

    reader.close();
}

void readListOfNumbers(thrust::host_vector<fptype>& target, const char* fname) {
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

void writeListOfNumbers(thrust::host_vector<fptype>& target, const char* fname) {
    std::ofstream writer;
    writer.open(fname);

    for(fptype t : target) {
        writer << t << " ";
    }

    writer.close();
}
