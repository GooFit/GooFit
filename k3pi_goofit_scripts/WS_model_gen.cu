//Generate WS signal only only toys
// GooFit stuff
#include <fstream>
#include <goofit/Application.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/physics/Amp4Body_TD.h>
#include <goofit/PDFs/physics/Lineshapes.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/TruthResolution.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <goofit/FitManager.h>
#include <thrust/count.h>
#include "ChristophModel.h"

using namespace std;
using namespace GooFit;

int main(int argc, char **argv) {
  // set constants
  const fptype MODEL_MAX_WEIGHT = 35.0;
  const unsigned int MODEL_MC_EVENTS_NORM = 1;

  // read in user input
  Application app(argc, argv);
  int BatchSize;
  app.add_option("batch,-b,--batch", BatchSize, "Batch size")->required();
  unsigned int genEvts;
  app.add_option("gen,-g,--gen", genEvts, "Number of events to generate")
      ->required();
  unsigned int seed;
  app.add_option("seed,-s,--seed", seed, "Generator random seed)")->required();
  std::string root_filename; // 4
  app.add_option("file,-f,--file", root_filename, "Output filename")
      ->required()
      ->check(CLI::NonexistentPath);
  fptype xmixing_value; // 5
  app.add_option("xmix,-x,--xmix", xmixing_value, "x mixing value")->required();
  fptype ymixing_value; // 6
  app.add_option("ymix,-y,--ymix", ymixing_value, "y mixing value")->required();
  GOOFIT_PARSE(app, argc, argv);

  // print status
  std::cout << "Generating toy signal MC using Christoph's K3PI model with ..."
            << std::endl;
  std::cout << "                    Batch size: " << BatchSize << ::std::endl;
  std::cout << "     Num of events to generate: " << genEvts << std::endl;
  std::cout << "         Generator random seed: " << seed << std::endl;
  std::cout << "                   Output file: " << root_filename << std::endl;
  std::cout << "                x mixing value: " << xmixing_value << std::endl;
  std::cout << "                y mixing value: " << ymixing_value << std::endl;

    // build model
  ChristophModel::ChristophModel model(xmixing_value, ymixing_value, MODEL_MC_EVENTS_NORM);
  model.setModelMaxWeight(MODEL_MAX_WEIGHT);

  // generate events in batches
  int generatedEvents = 0;
  int RunNum = 0;
  unsigned int generationOffset = 0;    

  while (generatedEvents < genEvts) {
    // update generation offset
    model.setGenerationOffset(generationOffset);

    // generate events for batch
    unsigned int keptEvts = 0;
    std::cout << "Generating events for this run using random seed " << seed
              << " and generator offset " << generationOffset << "..."
              << std::endl;
    auto tuple = model.generateSig(BatchSize, seed);
    auto particles = std::get<0>(tuple);
    auto variables = std::get<1>(tuple);
    // add accepted events to ntuple
    auto flags = std::get<3>(tuple);
    int accepted =
        thrust::count_if(flags.begin(), flags.end(), thrust::identity<bool>());
    ++RunNum;
    // print status and update generation offset to be used for next batch
    generationOffset += BatchSize;
    fmt::print("Run # {}: x={:.6} y={:.6} Using accept-reject method leaves "
               "you with {} out of {} events. {:.4} percent of total.\n",
               RunNum, xmixing_value, ymixing_value, keptEvts, BatchSize,
               generatedEvents * 100.0 / genEvts);
  }
return 0;
}