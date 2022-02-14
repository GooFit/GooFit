

// System stuff
#include <fstream>
#include <sys/time.h>
#include <sys/times.h>

// GooFit stuff
#include <goofit/Application.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/PDFs/physics/Amp3Body.h>
#include <goofit/PDFs/physics/DalitzPlotter.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <goofit/utilities/Style.h>

using namespace std;
using namespace GooFit;

Variable fixedRhoMass("rho_mass", 0.7758, 0.01, 0.7, 0.8);
Variable fixedRhoWidth("rho_width", 0.1503, 0.01, 0.1, 0.2);

const fptype _mD0       = 1.86484;
const fptype _mD02      = _mD0 * _mD0;
const fptype _mD02inv   = 1. / _mD02;
const fptype piPlusMass = 0.13957018;
const fptype piZeroMass = 0.1349766;

// Constants used in more than one PDF component.
Variable motherM("motherM", _mD0);
Variable chargeM("chargeM", piPlusMass);
Variable neutrlM("neutrlM", piZeroMass);
Variable massSum("massSum", _mD0 *_mD0 + 2 * piPlusMass * piPlusMass + piZeroMass * piZeroMass); // = 3.53481
Variable constantOne("constantOne", 1);
Variable constantZero("constantZero", 0);

fptype cpuGetM23(fptype massPZ, fptype massPM) {
    return (_mD02 + piZeroMass * piZeroMass + piPlusMass * piPlusMass + piPlusMass * piPlusMass - massPZ - massPM);
}


Amp3Body *makeSignalPdf(Observable m12, Observable m13, EventNumber eventNumber, GooPdf *eff = 0) {
    DecayInfo3 dtop0pp;
    dtop0pp.motherMass   = _mD0;
    dtop0pp.daug1Mass    = piZeroMass;
    dtop0pp.daug2Mass    = piPlusMass;
    dtop0pp.daug3Mass    = piPlusMass;
    dtop0pp.meson_radius = 1.5;

    ResonancePdf *rhop = new Resonances::RBW(
        "rhop", Variable("rhop_amp_real", 1), Variable("rhop_amp_imag", 0), fixedRhoMass, fixedRhoWidth, 1, PAIR_12);

    bool fixAmps = false; // Takes ~400x longer

    ResonancePdf *rhom = new Resonances::RBW(
        "rhom",
        fixAmps ? Variable("rhom_amp_real", 0.714) : Variable("rhom_amp_real", 0.714, 0.001, 0, 0),
        fixAmps ? Variable("rhom_amp_imag", -0.025) : Variable("rhom_amp_imag", -0.025, 0.1, 0, 0),
        fixedRhoMass,
        fixedRhoWidth,
        1,
        PAIR_13);

    ResonancePdf *rho0 = new Resonances::RBW(
        "rho0",
        fixAmps ? Variable("rho0_amp_real", 0.565) : Variable("rho0_amp_real", 0.565, 0.001, 0, 0),
        fixAmps ? Variable("rho0_amp_imag", 0.164) : Variable("rho0_amp_imag", 0.164, 0.1, 0, 0),
        fixedRhoMass,
        fixedRhoWidth,
        1,
        PAIR_23);

    Variable sharedMass("rhop_1450_mass", 1.465, 0.01, 1.0, 2.0);
    Variable shareWidth("rhop_1450_width", 0.400, 0.01, 0.01, 5.0);

    ResonancePdf *rhop_1450 = new Resonances::RBW(
        "rhop_1450",
        fixAmps ? Variable("rhop_1450_amp_real", -0.174) : Variable("rhop_1450_amp_real", -0.174, 0.001, 0, 0),
        fixAmps ? Variable("rhop_1450_amp_imag", -0.117) : Variable("rhop_1450_amp_imag", -0.117, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_12);

    ResonancePdf *rho0_1450 = new Resonances::RBW(
        "rho0_1450",
        fixAmps ? Variable("rho0_1450_amp_real", 0.325) : Variable("rho0_1450_amp_real", 0.325, 0.001, 0, 0),
        fixAmps ? Variable("rho0_1450_amp_imag", 0.057) : Variable("rho0_1450_amp_imag", 0.057, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_23);

    ResonancePdf *rhom_1450 = new Resonances::RBW(
        "rhom_1450",
        fixAmps ? Variable("rhom_1450_amp_real", 0.788) : Variable("rhom_1450_amp_real", 0.788, 0.001, 0, 0),
        fixAmps ? Variable("rhom_1450_amp_imag", 0.226) : Variable("rhom_1450_amp_imag", 0.226, 0.1, 0, 0),
        sharedMass,
        shareWidth,
        1,
        PAIR_13);

    Variable sharedMass2("rhop_1700_mass", 1.720, 0.01, 1.6, 1.9);
    Variable shareWidth2("rhop_1700_width", 0.250, 0.01, 0.1, 1.0);

    ResonancePdf *rhop_1700 = new Resonances::RBW(
        "rhop_1700",
        fixAmps ? Variable("rhop_1700_amp_real", 2.151) : Variable("rhop_1700_amp_real", 2.151, 0.001, 0, 0),
        fixAmps ? Variable("rhop_1700_amp_imag", -0.658) : Variable("rhop_1700_amp_imag", -0.658, 0.1, 0, 0),
        sharedMass2,
        shareWidth2,
        1,
        PAIR_12);

    ResonancePdf *rho0_1700 = new Resonances::RBW(
        "rho0_1700",
        fixAmps ? Variable("rho0_1700_amp_real", 2.400) : Variable("rho0_1700_amp_real", 2.400, 0.001, 0, 0),
        fixAmps ? Variable("rho0_1700_amp_imag", -0.734) : Variable("rho0_1700_amp_imag", -0.734, 0.1, 0, 0),
        sharedMass2,
        shareWidth2,
        1,
        PAIR_23);

    ResonancePdf *rhom_1700 = new Resonances::RBW(
        "rhom_1700",
        fixAmps ? Variable("rhom_1700_amp_real", 1.286) : Variable("rhom_1700_amp_real", 1.286, 0.001, 0, 0),
        fixAmps ? Variable("rhom_1700_amp_imag", -1.532) : Variable("rhom_1700_amp_imag", -1.532, 0.1, 0, 0),
        sharedMass2,
        shareWidth2,
        1,
        PAIR_13);

    ResonancePdf *f0_980 = new Resonances::RBW("f0_980",
                                               fixAmps ? Variable("f0_980_amp_real", 0.008 * (-_mD02))
                                                       : Variable("f0_980_amp_real", 0.008 * (-_mD02), 0.001, 0, 0),
                                               fixAmps ? Variable("f0_980_amp_imag", -0.013 * (-_mD02))
                                                       : Variable("f0_980_amp_imag", -0.013 * (-_mD02), 0.1, 0, 0),
                                               Variable("f0_980_mass", 0.980, 0.01, 0.8, 1.2),
                                               Variable("f0_980_width", 0.044, 0.001, 0.001, 0.08),
                                               0,
                                               PAIR_23);

    ResonancePdf *f0_1370 = new Resonances::RBW("f0_1370",
                                                fixAmps ? Variable("f0_1370_amp_real", -0.058 * (-_mD02))
                                                        : Variable("f0_1370_amp_real", -0.058 * (-_mD02), 0.001, 0, 0),
                                                fixAmps ? Variable("f0_1370_amp_imag", 0.026 * (-_mD02))
                                                        : Variable("f0_1370_amp_imag", 0.026 * (-_mD02), 0.1, 0, 0),
                                                Variable("f0_1370_mass", 1.434, 0.01, 1.2, 1.6),
                                                Variable("f0_1370_width", 0.173, 0.01, 0.01, 0.4),
                                                0,
                                                PAIR_23);

    ResonancePdf *f0_1500 = new Resonances::RBW("f0_1500",
                                                fixAmps ? Variable("f0_1500_amp_real", 0.057 * (-_mD02))
                                                        : Variable("f0_1500_amp_real", 0.057 * (-_mD02), 0.001, 0, 0),
                                                fixAmps ? Variable("f0_1500_amp_imag", 0.012 * (-_mD02))
                                                        : Variable("f0_1500_amp_imag", 0.012 * (-_mD02), 0.1, 0, 0),
                                                Variable("f0_1500_mass", 1.507, 0.01, 1.3, 1.7),
                                                Variable("f0_1500_width", 0.109, 0.01, 0.01, 0.3),
                                                0,
                                                PAIR_23);

    ResonancePdf *f0_1710 = new Resonances::RBW("f0_1710",
                                                fixAmps ? Variable("f0_1710_amp_real", 0.070 * (-_mD02))
                                                        : Variable("f0_1710_amp_real", 0.070 * (-_mD02), 0.001, 0, 0),
                                                fixAmps ? Variable("f0_1710_amp_imag", 0.087 * (-_mD02))
                                                        : Variable("f0_1710_amp_imag", 0.087 * (-_mD02), 0.1, 0, 0),
                                                Variable("f0_1710_mass", 1.714, 0.01, 1.5, 2.9),
                                                Variable("f0_1710_width", 0.140, 0.01, 0.01, 0.5),
                                                0,
                                                PAIR_23);

    ResonancePdf *f2_1270
        = new Resonances::RBW("f2_1270",
                              fixAmps ? Variable("f2_1270_amp_real", -1.027 * (-_mD02inv))
                                      : Variable("f2_1270_amp_real", -1.027 * (-_mD02inv), 0.001, 0, 0),
                              fixAmps ? Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv))
                                      : Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv), 0.1, 0, 0),
                              Variable("f2_1270_mass", 1.2754, 0.01, 1.0, 1.5),
                              Variable("f2_1270_width", 0.1851, 0.01, 0.01, 0.4),
                              2,
                              PAIR_23);

    ResonancePdf *f0_600 = new Resonances::RBW("f0_600",
                                               fixAmps ? Variable("f0_600_amp_real", 0.068 * (-_mD02))
                                                       : Variable("f0_600_amp_real", 0.068 * (-_mD02), 0.001, 0, 0),
                                               fixAmps ? Variable("f0_600_amp_imag", 0.010 * (-_mD02))
                                                       : Variable("f0_600_amp_imag", 0.010 * (-_mD02), 0.1, 0, 0),
                                               Variable("f0_600_mass", 0.500, 0.01, 0.3, 0.7),
                                               Variable("f0_600_width", 0.400, 0.01, 0.2, 0.6),
                                               0,
                                               PAIR_23);

    ResonancePdf *nonr = new Resonances::NonRes(
        "nonr",
        fixAmps ? Variable("nonr_amp_real", 0.5595 * (-1)) : Variable("nonr_amp_real", 0.5595 * (-1), 0.001, 0, 0),
        fixAmps ? Variable("nonr_amp_imag", -0.108761 * (-1)) : Variable("nonr_amp_imag", -0.108761 * (-1), 0.1, 0, 0));



    std::vector<Variable> beta_r = {Variable("beta1_r", 9.3*cos(-78.7/180)),
            Variable("beta2_r", 10.89*cos(-159.1/180)),
            Variable("beta3_r", 24.20*cos(168.0/180)),
            Variable("beta4_r", 9.15*cos(90.5/180)),
            Variable("beta5_r", 0)};

     std::vector<Variable> beta_i =  {Variable("beta1_i", 9.3*sin(-78.7/180)),
              Variable("beta2_i", 10.89*sin(-159.1180)),
              Variable("beta3_i", 24.20*sin(168.0180)),
              Variable("beta4_i", 9.15*sin(90.5180)),
              Variable("beta5_i", 0)};
     std::vector<Variable> fprod_r =  {Variable("f_prod1_r", 7.94*cos(73.9/180)),
            Variable("f_prod2_r", 2.0*cos(-18.0/180)),
            Variable("f_prod3_r", 5.1*cos(33.0/180)),
            Variable("f_prod4_r", 3.23*cos(4.8/180)),
            Variable("f_prod5_r", 0)};
     std::vector<Variable> fprod_i =  {Variable("f_prod1_i", 7.94*sin(73.9/180)),
            Variable("f_prod2_i", 2.0*sin(-18.0/180)),
            Variable("f_prod3_i", 5.1*sin(33.0/180)),
            Variable("f_prod4_i", 3.23*sin(4.8/180)), 
            Variable("f_prod5_i", 0)};
     std::vector<Variable> poles =  {Variable("g00", 0.22889),
            Variable("g01", -0.55377),
            Variable("g02", 0),
            Variable("g03", -0.39899),
            Variable("g04", -0.34639),
            Variable("m0", 0.65100),
            Variable("g10", 0.94128),
            Variable("g11", 0.55095),
            Variable("g12", 0),
            Variable("g13",  0.39065),
            Variable("g14", 0.31503),
            Variable("m1", 1.20360),
            Variable("g20", 0.36856),
            Variable("g21", 0.23888),
            Variable("g22", 0.55639),
            Variable("g23", 0.18340),
            Variable("g24", 0.18681),
            Variable("m2", 1.55817),
            Variable("g30", 0.33650),
            Variable("g31", 0.40907),
            Variable("g32", 0.85679),
            Variable("g33",  0.19906),
            Variable("g34", -0.00984),
            Variable("m3", 1.21000), 
            Variable("g40", 0.18171),
            Variable("g41", -0.17558),
            Variable("g42", -0.79658),
            Variable("g43", -0.00355),
            Variable("g44", 0.22358),
            Variable("m4", 1.82206)};
 std::vector<Variable> fscat =  { Variable("f11_scatt", 0.23399),
             Variable("f12_scatt", 0.15044),
             Variable("f13_scatt", -0.20545),
             Variable("f14_scatt", 0.32825),
             Variable("f15_scatt", 0.35412)};
    

    ResonancePdf * kMatrix = new Resonances::kMatrix("kMatrix",
                                 Variable("a_r", 8.5, 0.0001, 0, 0),
                                 Variable("a_i", 68.5, 0.0001, 0, 0),
                                 Variable("sA0", -0.15),
                                 Variable("sA", 1),
                                 Variable("s0_prod", -0.07),
                                 Variable("s0_scatt", -3.92637),
                                 beta_r,
                                 beta_i,
                                 fprod_r,
                                 fprod_i,
                                 fscat,                                 
                                 poles,
                                 0, 
                                 PAIR_23);





    dtop0pp.resonances.push_back(kMatrix);

    bool fitMasses = false;

    if(!fitMasses) {
        for(vector<ResonancePdf *>::iterator res = dtop0pp.resonances.begin(); res != dtop0pp.resonances.end(); ++res) {
            (*res)->setParameterConstantness(true);
        }
    }

    if(!eff) {
        // By default create a constant efficiency.
        vector<Variable> offsets;
        vector<Observable> observables;
        vector<Variable> coefficients;

        observables.push_back(m12);
        observables.push_back(m13);
        offsets.push_back(constantZero);
        offsets.push_back(constantZero);
        coefficients.push_back(constantOne);
        eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);
    }

    return new Amp3Body("signalPDF", m12, m13, eventNumber, dtop0pp, eff);
}

int runToyFit(Amp3Body *signal, UnbinnedDataSet *data) {
    // EXERCISE 1 (real part): Create a PolynomialPdf which models
    // the efficiency you imposed in the preliminary, and use it in constructing
    // the signal PDF.

    // EXERCISE 2: Create a K0 veto function and use it as the efficiency.

    // EXERCISE 3: Make the efficiency a product of the two functions
    // from the previous exercises.

    signal->setData(data);
    signal->setDataSize(data->getNumEvents());
    FitManager datapdf(signal);

    datapdf.fit();

    ProdPdf prodpdf{"prodpdf", {signal}};

    DalitzPlotter plotter(&prodpdf, signal);



    return datapdf;
}

int main(int argc, char **argv) {
    GooFit::Application app("Dalitz example", argc, argv);

    print_goofit_info();

    // Observables setup
    Observable m12("m12", 0, 3);
    Observable m13("m13", 0, 3);
    EventNumber eventNumber("eventNumber");
    m12.setNumBins(1);
    m13.setNumBins(2);

    // Prepare the data
    UnbinnedDataSet data({m12, m13, eventNumber});

    for(int i=0; i<m12.getNumBins(); i++) {
        for(int j=0; j<m13.getNumBins(); j++) {
           m12.setValue(m12.getLowerLimit() + m12.getBinSize() * (i + 0.5));
           m13.setValue(m13.getLowerLimit() + m13.getBinSize() * (j + 0.5));
           data.addEvent();
            eventNumber.setValue(eventNumber.getValue() + 1);
        }
    }

                

    // Set up the model
    Amp3Body *signal = makeSignalPdf(m12, m13, eventNumber);

    // A wrapper for plotting without complex number segfault
    ProdPdf prodpdf{"prodpdf", {signal}};
    prodpdf.setData(&data);
    signal->setDataSize(data.getNumEvents());
    prodpdf.normalise();


}
