#include <iostream>
#include <fstream>
#include <string>

#include "ToyModel.h"

// GooFit stuff
#include "goofit/Application.h"
#include <algorithm>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <mcbooster/functors/FlagAcceptReject.h>
#include <numeric>
#include <random>
#include <thrust/count.h>

using namespace GooFit;

int main(int argc, char **argv) {
    // Set some constants.
    fptype MODEL_MAX_WEIGHT   = 35.0;
    uint MODEL_MC_EVENTS_NORM = 1;
    uint batchSize            = 20000000;
    uint genEvts              = 20000;
    uint seed                 = 123456;
    fptype xMixingValue       = 0.0042;
    fptype yMixingValue       = 0.0064;
    printf("Generating %i events with x=%f and y=%f\n", genEvts, xMixingValue, yMixingValue);

    // Build the model.
    ToyModel model(xMixingValue, yMixingValue, MODEL_MC_EVENTS_NORM);
    model.setModelMaxWeight(MODEL_MAX_WEIGHT);

    // Setup the output file.
    std::string fname = "tddp4_data.txt";
    std::ofstream outputfile;
    outputfile.open(fname);
    double tm12, tm34, tc12, tc34, tphi, tdtime;

    // Generate events in batches.
    uint generatedEvents  = 0;
    uint runNum           = 0;
    uint generationOffset = 0;
    while(generatedEvents < genEvts) {
        // Update the generation offset.
        model.setGenerationOffset(generationOffset);

        // Generate events for batch.
        unsigned int keptEvts = 0;
        std::cout << "Generating events for this run using random seed " << seed << " and generator offset "
                  << generationOffset << "..." << std::endl;
        auto tuple     = model.generateSig(batchSize, seed);
        auto particles = std::get<0>(tuple);
        auto variables = std::get<1>(tuple);
        auto flags     = std::get<3>(tuple);
        for(int i = 0; i < flags.size(); ++i) {
            if(generatedEvents < genEvts && flags[i] == 1) {
                ++generatedEvents;
                ++keptEvts;
                tm12   = (*(variables[0]))[i];
                tm34   = (*(variables[1]))[i];
                tc12   = (*(variables[2]))[i];
                tc34   = (*(variables[3]))[i];
                tphi   = (*(variables[4]))[i];
                tdtime = (*(variables[5]))[i];
                outputfile << tm12 << " " << tm34 << " " << tc12 << " " << tc34 << " " << tphi << " " << tdtime
                           << std::endl;
            }
        }

        generationOffset += batchSize;
        fmt::print("Run # {}: x={:.6} y={:.6} Using accept-reject method leaves "
                   "you with {} out of {} events. {:.4} percent of total.\n",
                   runNum,
                   xMixingValue,
                   yMixingValue,
                   keptEvts,
                   batchSize,
                   generatedEvents * 100.0 / genEvts);

        delete variables[0];
        delete variables[1];
        delete variables[2];
        delete variables[3];
        delete variables[4];
        delete variables[5];
        delete particles[0];
        delete particles[1];
        delete particles[2];
        delete particles[3];
    }

    outputfile.close();
    return 0;
}
