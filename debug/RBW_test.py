import sys
from goofit import *
from math import sqrt, pi, cos, sin

print_goofit_info()

def makeSignalPdf(m12, m13, eventNumber):
    # Now the amplitudes are fixed and the masses are floating within 3 sigma of PDG value
    d0toks0pipi = DecayInfo3()
    d0toks0pipi.motherMass = 1.86484
    d0toks0pipi.daug1Mass = 0.497614
    d0toks0pipi.daug2Mass = 0.139570
    d0toks0pipi.daug3Mass = 0.139570

    f2_1270_width = Variable("f2_1270_width", 0.1867, 0.0001, 0.1867-3*0.0025, 0.1867+3*0.0022)
    f2_1270_mass = Variable("f2_1270_mass", 1.2755, 0.0001, 1.250, 1.300)


    f2_1270_amp_real = Variable("f2_1270_amp_real", 1.43*cos(-36.3*pi/180)) 
    f2_1270_amp_imag = Variable("f2_1270_amp_imag", 1.43*sin(-36.3*pi/180)) 


    myRBW= Resonances.RBW("f2_1270",
                                 f2_1270_amp_real,
                             f2_1270_amp_imag,
                             f2_1270_mass,
                             f2_1270_width,
                             2,
                             PAIR_23)
    print('Add to list of resonances.')
    d0toks0pipi.resonances = ([myRBW])



    print('Define efficiency and Amp3Body Pdf.')
    # Constant efficiency
    constantOne = Variable('One', 1)
    constantZero = Variable('Zero', 0)
    eff = PolynomialPdf('constantEff', [m12, m13], [constantOne], [constantZero,constantZero], 0)
    d = Amp3Body('signalPDF', m12, m13, eventNumber, d0toks0pipi, eff)

    print(d.getName())
    print(d.getParameters())
    print(d.getObservables())

    return d

def inDalitz(m12, m13):
    mD0, mKS, mPi = 1.86486, 0.497614, 0.13957

    if m12 < pow(mKS + mPi, 2) or m12 > pow(mD0 - mPi, 2):
        return False

    # Calculate energies of particles 1 and 3 in m12 rest frame
    E1 = 0.5 * (m12 - mPi*mPi + mKS*mKS) / sqrt(m12)
    E3 = 0.5 * (mD0*mD0 - m12 - mPi*mPi) / sqrt(m12)    

    minimum = pow(E1 + E3, 2) - pow(sqrt(E1*E1 - mKS*mKS) + sqrt(E3*E3 - mPi*mPi), 2)
    maximum = pow(E1 + E3, 2) - pow(sqrt(E1*E1 - mKS*mKS) - sqrt(E3*E3 - mPi*mPi), 2)
    if m13 < minimum or m13 > maximum:
        return False

    return True

def make_toy_data(dp, m12, m13, eventNumber):
    data = UnbinnedDataSet(m12, m13, eventNumber)

    print('ProdPdf')
    prod = ProdPdf('prod', [dp])

    print(m12.getNumBins())
    print(m13.getNumBins())

    xbins, ybins = [], []
    for i in range(m12.getNumBins()):
        m12.setValue(m12.getLowerLimit() + m12.getBinSize() * (i + 0.5))
        for j in range(m13.getNumBins()):
            m13.setValue(m13.getLowerLimit() + m13.getBinSize() * (j + 0.5))
            if inDalitz(m12.getValue(), m13.getValue()):
                xbins.append(i)
                ybins.append(j)
                data.addEvent()
                eventNumber.setValue(eventNumber.getValue() + 1)
    prod.setData(data)
    dp.setDataSize(data.getNumEvents())

    print('Normalise')
    prod.normalize()
    pdfvals = prod.getCompProbsAtDataPoints() # Memory error
    print (pdfvals)

    return data

def main():

    nBins = 10
    m12 = Observable('m12', 0, 3)
    m13 = Observable('m13', 0, 3)
    m12.setNumBins(nBins)
    m13.setNumBins(nBins)
    eventNumber = EventNumber('eventNumber')

    signal = makeSignalPdf(m12, m13, eventNumber)
    dataset = make_toy_data(signal, m12, m13, eventNumber)

if __name__ == '__main__':
    sys.exit(main())

