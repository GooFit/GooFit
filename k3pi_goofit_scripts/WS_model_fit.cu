// This file fits the signal toy MC samples generated using FinalGen0.cu with
// Christoph's model

#include <algorithm>
#include <cctype>
#include <ctime>
#include <functional>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <thrust/count.h>
#include <vector>

#include "ChristophModel.h"

// GooFit stuff
#include <goofit/Application.h>

using namespace GooFit;

// param filesString: comma separated list of input root files
// returns a vector of the input root files
std::vector<std::string> buildRootFileList(const std::string &filesString) {
    // make a copy so we don't modify the input string
    std::string filesStringCopy = filesString;
  
    // first remove any whitespace
    filesStringCopy.erase(
        std::remove_if(filesStringCopy.begin(), filesStringCopy.end(), ::isspace),
        filesStringCopy.end());
  
    // now separate on commas and build result vector
    std::stringstream ss(filesStringCopy);
    std::vector<std::string> result;
    while (ss.good()) {
      std::string substr;
      getline(ss, substr, ',');
      result.push_back(substr);
    }
  
    return result;
  }

  int main(int argc, char **argv) {
    // set constants
    // cudaSetDevice(1);
  
    const fptype START_MIXING_ERROR = 0.00001;
    const fptype MIXING_LOWER_LIMIT = -0.15;
    const fptype MIXING_UPPER_LIMIT = 0.15;
    const unsigned int MODEL_MC_EVENTS_NORM = 15e6;
  
    // read in user input
    Application app(argc, argv);
    unsigned int sampleNum;
    app.add_option("sampleNum, -n", sampleNum, "Sample number of input data")
        ->required();
    fptype xmixingStart;
    app.add_option("xmixingStart, --xStart", xmixingStart,
                   "Starting value for x mixing fit variable")
        ->required();
    fptype ymixingStart;
    app.add_option("ymixingStart, --yStart", ymixingStart,
                   "Starting value for y mixing fit variable")
        ->required();
    bool special_integral;
    app.add_flag("--specInt", special_integral,"Perform fit with 6D numerical integral and per-event efficiencies instead of 5D + decay time");
    std::string inputFilenames;
    app.add_option("inFiles, -f", inputFilenames,
                   "Comma separated list of root files with generated data")
        ->required();
    std::string outputFilename;
    app.add_option("outFile, -o", outputFilename,
                   "Fit log file to append fit results to")
        ->required();
    bool floatDCSCoeffs;
    app.add_flag("--floatDCS", floatDCSCoeffs,
                 "Use this flag if want to float DCS amplitude coefficents");
    unsigned int normalisation_seed;
    app.add_option("generationOffset, --genOff", normalisation_seed,"Set the seed used to generate normalisation events (particularly for 6D integration)");
    // parse user input
    GOOFIT_PARSE(app, argc, argv);
    std::vector<std::string> inputFiles = buildRootFileList(inputFilenames);
  
    // print status
    std::cout << "Fitting toy signal MC using Christoph's K3PI model with ..."
              << std::endl;
    std::cout << "                                Sample #: " << sampleNum
              << ::std::endl;
    std::cout << "Starting value for x mixing fit variable: " << xmixingStart
              << std::endl;
    std::cout << "Starting value for y mixing fit variable: " << ymixingStart
              << std::endl;
    std::cout << "                             Floating DCS amplitude coefficients: " << floatDCSCoeffs
              << std::endl;
    std::cout << "                             Calculating normalisation completely numerically: " << special_integral
              << std::endl;
    std::cout << "           Root files with generated data: " << std::endl;
    std::cout << "           Calculating normalisation events with seed: " << normalisation_seed << std::endl;
    for (auto const &f : inputFiles) {
      std::cout << f << std::endl;
    }
    std::cout << "       Appending fit results to log file: " << outputFilename
              << std::endl;
  
    // build model
    std::cout << "Building model..." << std::endl;
    ChristophModel model(xmixingStart, ymixingStart, MODEL_MC_EVENTS_NORM,special_integral,normalisation_seed);
    model.setXMixingRangeForFit(START_MIXING_ERROR, MIXING_LOWER_LIMIT,
                                MIXING_UPPER_LIMIT);
    model.setYMixingRangeForFit(START_MIXING_ERROR, MIXING_LOWER_LIMIT,
                                MIXING_UPPER_LIMIT);
    if (floatDCSCoeffs == true)
    {
        //std::vector<fptype> startValuesForDCSCoeffs = K3PiUtilities::genRandomStartValsForDCSCoeffs(randStartValueGen);
        model.allowDCSCoeffFloat(0.00001);
    }
    // read input files
    unsigned int fitEvts = 0;
    for (auto const &inputFile : inputFiles) {
      // build reader for root file
      std::cout << "Reading File" << std::endl;
      std::ifstream infile(inputFile.c_str());
      double m12,m34,c12,c34,phi,dtime;
      // read in events
      std::cout << "Reading in events from file " << inputFile << "..."
                << std::endl;
      while (infile >> m12 >> m34 >> c12 >> c34 >> phi >> dtime) {

        model.addEventToCurrentDataToFit(m12, m34, c12,
                                         c34, phi, dtime, 0.0,
                                         fitEvts++);
  
        if (fitEvts % 10000 == 0) {
          std::cout << "Read in " << fitEvts << " events so far." << std::endl;
        }
      } // end while loop for event reader
  
      // close root file
      std::cout << "Closing " << inputFile << "..." << std::endl;
      infile.close();
    } // end loop over input files
    std::cout << "Read in " << fitEvts << " events from " << inputFiles.size()
              << " files." << std::endl;
  
    // do fitting
    model.fitCurrentData(sampleNum, outputFilename);
  
    std::cout << "Done." << std::endl;
    return 0;
  }
  