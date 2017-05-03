#include "goofit/FunctorWriter.h"
#include <fstream>
#include <map>
#include "goofit/PdfBase.h"
#include "goofit/Variable.h"

void writeToFile(PdfBase* pdf, const char* fname) {
    Variable_v params = pdf->getParameters();

    std::ofstream writer;
    writer.open(fname);

    for(Variable* p : params) {
        writer << p->name << " "
               << p->value << " "
               << p->error << " "
               << p->numbins << " "
               << p->lowerlimit << " "
               << p->upperlimit
               << std::endl;
    }

    writer.close();
}


void readFromFile(PdfBase* pdf, const char* fname) {
    Variable_v params = pdf->getParameters();

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

        if(var) {
            reader >> var->value
                   >> var->error
                   >> var->numbins
                   >> var->lowerlimit
                   >> var->upperlimit;

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
