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

    # beta coeeficients
    beta_r = [Variable("beta1_r", 9.3*cos(-78.7*pi/180)),
              Variable("beta2_r", 10.89*cos(-159.1*pi/180)),
              Variable("beta3_r", 24.20*cos(168.0*pi/180)),
              Variable("beta4_r", 9.15*cos(90.5*pi/180)),
              Variable("beta5_r", 0)
             ]

    beta_i = [Variable("beta1_i", 9.3*sin(-78.7*pi/180)),
              Variable("beta2_i", 10.89*sin(-159.1*pi/180)),
              Variable("beta3_i", 24.20*sin(168.0*pi/180)),
              Variable("beta4_i", 9.15*sin(90.5*pi/180)),
              Variable("beta5_i", 0)
             ]
    
    # f_prod coefficients
    f_prod_r = [Variable("f_prod1_r", 7.94*cos(73.9*pi/180)),
                Variable("f_prod2_r", 2.0*cos(-18.0*pi/180)),
                Variable("f_prod3_r", 5.1*cos(33.0*pi/180)),
                Variable("f_prod4_r", 3.23*cos(4.8*pi/180)),
                Variable("f_prod5_r", 0)
               ]

    f_prod_i = [Variable("f_prod1_i", 7.94*sin(73.9*pi/180)),
                Variable("f_prod2_i", 2.0*sin(-18.0*pi/180)),
                Variable("f_prod3_i", 5.1*sin(33.0*pi/180)),
                Variable("f_prod4_i", 3.23*sin(4.8*pi/180)), 
                Variable("f_prod5_i", 0)                     
               ]
    
    

    # kMatrix
    fscat = [Variable("f11_scatt", 0.23399),
             Variable("f12_scatt", 0.15044),
             Variable("f13_scatt", -0.20545),
             Variable("f14_scatt", 0.32825),
             Variable("f15_scatt", 0.35412)]
    poles = [Variable("g00", 0.22889),
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
             Variable("m4", 1.82206)]
    print('Define K Matrix Resonance.')
    print("beta_r =", beta_r)
    print("beta_i =", beta_i)
    print("f_prod_r =", f_prod_r)
    print("f_prod_i =", f_prod_i)

    kMatrix = Resonances.kMatrix("kMatrix",
                                 Variable("a_r", 8.5, 0.0001, 0, 0),
                                 Variable("a_i", 68.5, 0.0001, 0, 0),
                                 Variable("sA0", -0.15),
                                 Variable("sA", 1),
                                 Variable("s0_prod", -0.07),
                                 Variable("s0_scatt", -3.92637),
                                 beta_r,
                                 beta_i,
                                 f_prod_r,
                                 f_prod_i,
                                 fscat,
                                 poles,
                                 0, # L (spin)
                                 PAIR_23
    )
    print('Add to list of resonances.')
    d0toks0pipi.resonances = ([kMatrix])

    print(kMatrix.get_amp_real())
    print(kMatrix.get_amp_img())
    print(kMatrix.getParameters())
    print(kMatrix.getObservables())

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


