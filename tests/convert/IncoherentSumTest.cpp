#include <goofit/Catch.h>

#include "testhelpers.h"
#include <goofit/PDFs/basic/GaussianPdf.h>
#include <goofit/PDFs/physics/Amp3Body_IS.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/UnbinnedDataSet.h>

#include <goofit/Variable.h>

#include <random>

using namespace GooFit;

TEST_CASE("Amp3Body_IS", "[convert][fit]") {
    const double _mD0     = 1.86484;
    const double _mD02    = _mD0 * _mD0;
    const double _mD02inv = 1.0 / _mD02;

    // Random number generation
    std::mt19937 gen(137);
    std::exponential_distribution<> d(1.5);

    // Independent variable.
    Observable m12{"m12", 0, 10};
    Observable m13{"m13", 0, 10};
    EventNumber evtNum{"evtNum", 0, INT_MAX};

    std::vector<Observable> var_list;
    var_list.push_back(m12);
    var_list.push_back(m13);
    var_list.push_back(evtNum);

    // Data set
    UnbinnedDataSet data(var_list);

    // Generate toy events.
    for(int i = 0; i < 100000; ++i) {
        double m12_val = d(gen);
        double m13_val = d(gen);

        if(m12_val > 0 && m12_val < 10 && m13_val > 0 && m13_val < 10) {
            m12.setValue(m12_val);
            m13.setValue(m13_val);
            evtNum.setValue(data.getNumEvents());
            data.addEvent();
        }
    }

    // Fit parameter
    DecayInfo3 di;
    di.motherMass   = 1.86484;
    di.daug1Mass    = 0.1349766;
    di.daug2Mass    = 0.13957018;
    di.daug3Mass    = 0.13957018;
    di.meson_radius = 1.5;

    // All variables
    Variable fixedRhoMass{"rho_mass", 0.7758, 0.01, 0.7, 0.8};
    Variable fixedRhoWidth{"rho_width", 0.1503, 0.01, 0.1, 0.2};

    Variable rhop_amp_real{"rhop_amp_real", 1};
    Variable rhop_amp_imag{"rhop_amp_imag", 0};

    Variable rhom_amp_real{"rhom_amp_real", 0.714, 0.001, 0, 0};
    Variable rhom_amp_imag{"rhom_amp_imag", -0.025, 0.1, 0, 0};

    Variable rho0_amp_real{"rho0_amp_real", 0.565, 0.001, 0, 0};
    Variable rho0_amp_imag{"rho0_amp_imag", 0.164, 0.1, 0, 0};

    Variable shared_1450_mass{"rhop_1450_mass", 1.465, 0.01, 1.0, 2.0};
    Variable shared_1450_width{"rhop_1450_width", 0.4, 0.01, 0.01, 5.0};

    Variable rhop_1450_real{"rhop_1450_real", -0.174, 0.001, 0, 0};
    Variable rhop_1450_imag{"rhop_1450_imag", -0.117, 0.1, 0, 0};

    Variable rho0_1450_real{"rhop_1450_real", 0.325, 0.001, 0, 0};
    Variable rho0_1450_imag{"rhop_1450_imag", 0.057, 0.1, 0, 0};

    Variable rhom_1450_real{"rhop_1450_real", 0.788, 0.001, 0, 0};
    Variable rhom_1450_imag{"rhop_1450_imag", 0.226, 0.1, 0, 0};

    Variable shared_1700_mass{"rhop_1700_mass", 1.720, 0.01, 1.6, 1.9};
    Variable shared_1700_width{"rhop_1700_width", 0.25, 0.01, 0.1, 1.0};

    Variable rhop_1700_real{"rhop_1700_real", 2.151, 0.001, 0, 0};
    Variable rhop_1700_imag{"rhop_1700_imag", -0.658, 0.1, 0, 0};

    Variable rho0_1700_real{"rhop_1700_real", 2.4, 0.001, 0, 0};
    Variable rho0_1700_imag{"rhop_1700_imag", -0.734, 0.1, 0, 0};

    Variable rhom_1700_real{"rhop_1700_real", 1.286, 0.001, 0, 0};
    Variable rhom_1700_imag{"rhop_1700_imag", -1.532, 0.1, 0, 0};

    Variable f0_980_real{"f0_980_amp_real", 0.008 * (-_mD02), 0.001, 0, 0};
    Variable f0_980_imag{"f0_980_amp_imag", -0.013 * (-_mD02), 0.1, 0, 0};

    Variable f0_980_mass{"f0_980_mass", 0.98, 0.01, 0.8, 1.2};
    Variable f0_980_width{"f0_980_width", 0.044, 0.001, 0.001, 0.08};

    Variable f0_1370_real{"f0_1370_amp_real", -0.058 * (-_mD02), 0.001, 0, 0};
    Variable f0_1370_imag{"f0_1370_amp_imag", 0.026 * (-_mD02), 0.1, 0, 0};

    Variable f0_1370_mass{"f0_1370_mass", 1.434, 0.01, 1.2, 1.6};
    Variable f0_1370_width{"f0_1370_width", 0.173, 0.01, 0.01, 0.4};

    Variable f0_1500_real{"f0_1500_amp_real", 0.057 * (-_mD02), 0.001, 0, 0};
    Variable f0_1500_imag{"f0_1500_amp_imag", 0.012 * (-_mD02), 0.1, 0, 0};

    Variable f0_1500_mass{"f0_1500_mass", 1.507, 0.01, 1.3, 1.7};
    Variable f0_1500_width{"f0_1500_width", 0.109, 0.01, 0.01, 0.3};

    Variable f0_1710_real{"f0_1710_amp_real", 0.07 * (-_mD02), 0.001, 0, 0};
    Variable f0_1710_imag{"f0_1710_amp_imag", 0.087 * (-_mD02), 0.1, 0, 0};

    Variable f0_1710_mass{"f0_1710_mass", 1.714, 0.01, 1.5, 2.9};
    Variable f0_1710_width{"f0_1710_width", 0.14, 0.01, 0.01, 0.5};

    Variable f0_1270_real{"f0_1270_amp_real", -1.027 * (-_mD02), 0.001, 0, 0};
    Variable f0_1270_imag{"f0_1270_amp_imag", -0.162 * (-_mD02), 0.1, 0, 0};

    Variable f0_1270_mass{"f0_1270_mass", 1.2754, 0.01, 1.0, 1.5};
    Variable f0_1270_width{"f0_1270_width", 0.1851, 0.01, 0.01, 0.4};

    Variable f0_600_real{"f0_600_amp_real", -1.027 * (-_mD02), 0.001, 0, 0};
    Variable f0_600_imag{"f0_600_amp_imag", -0.162 * (-_mD02), 0.1, 0, 0};

    Variable f0_600_mass{"f0_600_mass", 1.2754, 0.01, 1.0, 1.5};
    Variable f0_600_width{"f0_600_width", 0.1851, 0.01, 0.01, 0.4};

    Variable nonr_amp_real{"nonr_amp_real", 0.5595 * (-1), 0.001, 0, 0};
    Variable nonr_amp_imag{"nonr_amp_imag", -0.108761 * (-1), 0.1, 0, 0};

    // Create resonances here
    ResonancePdf *rhop
        = new Resonances::RBW("rhop", rhop_amp_real, rhop_amp_imag, fixedRhoMass, fixedRhoWidth, 1, PAIR_12);
    ResonancePdf *rhom
        = new Resonances::RBW("rhom", rhom_amp_real, rhom_amp_imag, fixedRhoMass, fixedRhoWidth, 1, PAIR_13);
    ResonancePdf *rho0
        = new Resonances::RBW("rho0", rho0_amp_real, rho0_amp_imag, fixedRhoMass, fixedRhoWidth, 1, PAIR_23);

    ResonancePdf *rhop_1450 = new Resonances::RBW(
        "rhop_1450", rhop_1450_real, rhop_1450_imag, shared_1450_mass, shared_1450_width, 1, PAIR_12);
    ResonancePdf *rho0_1450 = new Resonances::RBW(
        "rho0_1450", rho0_1450_real, rho0_1450_imag, shared_1450_mass, shared_1450_width, 1, PAIR_23);
    ResonancePdf *rhom_1450 = new Resonances::RBW(
        "rhom_1450", rhom_1450_real, rhom_1450_imag, shared_1450_mass, shared_1450_width, 1, PAIR_13);

    ResonancePdf *rhop_1700 = new Resonances::RBW(
        "rhop_1700", rhop_1700_real, rhop_1700_imag, shared_1700_mass, shared_1700_width, 1, PAIR_12);
    ResonancePdf *rho0_1700 = new Resonances::RBW(
        "rho0_1700", rho0_1700_real, rho0_1700_imag, shared_1700_mass, shared_1700_width, 1, PAIR_23);
    ResonancePdf *rhom_1700 = new Resonances::RBW(
        "rhom_1700", rhom_1700_real, rhom_1700_imag, shared_1700_mass, shared_1700_width, 1, PAIR_13);

    ResonancePdf *f0_980
        = new Resonances::RBW("f0_980", f0_980_real, f0_980_imag, f0_980_mass, f0_980_width, 0, PAIR_23);

    ResonancePdf *f0_1370
        = new Resonances::RBW("f0_1370", f0_1370_real, f0_1370_imag, f0_1370_mass, f0_1370_width, 0, PAIR_23);

    ResonancePdf *f0_1500
        = new Resonances::RBW("f0_1500", f0_1500_real, f0_1500_imag, f0_1500_mass, f0_1500_width, 0, PAIR_23);

    ResonancePdf *f0_1710
        = new Resonances::RBW("f0_1710", f0_1710_real, f0_1710_imag, f0_1710_mass, f0_1710_width, 0, PAIR_23);

    ResonancePdf *f0_1270
        = new Resonances::RBW("f0_1270", f0_1270_real, f0_1270_imag, f0_1270_mass, f0_1270_width, 2, PAIR_23);

    ResonancePdf *f0_600
        = new Resonances::RBW("f0_600", f0_600_real, f0_600_imag, f0_600_mass, f0_600_width, 0, PAIR_23);

    ResonancePdf *nonr = new Resonances::NonRes("nonr", nonr_amp_real, nonr_amp_imag);

    // Add resonances to Decay structure
    di.resonances.push_back(nonr);
    di.resonances.push_back(rhop);
    di.resonances.push_back(rho0);
    di.resonances.push_back(rhom);

    di.resonances.push_back(rhop_1450);
    di.resonances.push_back(rho0_1450);
    di.resonances.push_back(rhom_1450);

    di.resonances.push_back(rhop_1700);
    di.resonances.push_back(rho0_1700);
    di.resonances.push_back(rhom_1700);

    di.resonances.push_back(f0_980);

    di.resonances.push_back(f0_1370);

    di.resonances.push_back(f0_1500);

    di.resonances.push_back(f0_1710);

    di.resonances.push_back(f0_1270);

    di.resonances.push_back(f0_600);

    Variable alpha{"alpha", 1, 0.1, -10, 10};
    Variable sigma{"sigma", 1, 0.1, -10, 10};
    GaussianPdf eff{"gauss", m12, alpha, sigma};

    // GooPdf object
    Amp3Body_IS pdf{"pdf", m12, m13, evtNum, di, &eff};
    pdf.setData(&data);
    pdf.setDataSize(data.getNumEvents(), 3);

    bool fitter = test_fitter(&pdf);

    CHECK(fitter);
    // CHECK(alpha.getError() ==  .1);
    // CHECK(alpha.getValue() == Approx(0.665178392).margin(alpha.getError() * 3));
}
