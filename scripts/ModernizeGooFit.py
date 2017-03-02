#!/usr/bin/env python

from __future__ import print_function
import re
try:
    from plumbum import local, cli
except ImportError:
    print("This file uses the plumbum library. Install with pip or conda (user directory or virtual environment OK).")
    raise

conversion = [
    ('["<]cuda_runtime_api.hh?[">]', '"fakecuda/cuda_runtime_api.h"'),
    ('["<]driver_types.hh?[">]', '"fakecuda/driver_types.h"'),
    ('["<]host_defines.hh?[">]', '"fakecuda/host_defines.h"'),
    ('["<]Application.hh?[">]', '"goofit/Application.h"'),
    ('["<]BinnedDataSet.hh?[">]', '"goofit/BinnedDataSet.h"'),
    ('["<]DataSet.hh?[">]', '"goofit/DataSet.h"'),
    ('["<]Faddeeva.hh?[">]', '"goofit/Faddeeva.h"'),
    ('["<]FitControl.hh?[">]', '"goofit/FitControl.h"'),
    ('["<]FitManager.hh?[">]', '"goofit/FitManager.h"'),
    ('["<]FitManagerMinuit1.hh?[">]', '"goofit/FitManagerMinuit1.h"'),
    ('["<]FitManagerMinuit2.hh?[">]', '"goofit/FitManagerMinuit2.h"'),
    ('["<]FitManagerMinuit3.hh?[">]', '"goofit/FitManagerMinuit3.h"'),
    ('["<]FunctorWriter.hh?[">]', '"goofit/FunctorWriter.h"'),
    ('["<]GlobalCudaDefines.hh?[">]', '"goofit/GlobalCudaDefines.h"'),
    ('["<]PdfBase.hh?[">]', '"goofit/PdfBase.h"'),
    ('["<]UnbinnedDataSet.hh?[">]', '"goofit/UnbinnedDataSet.h"'),
    ('["<]Variable.hh?[">]', '"goofit/Variable.h"'),
    ('["<]AddPdf.hh?[">]', '"goofit/PDFs/AddPdf.h"'),
    ('["<]ArgusPdf.hh?[">]', '"goofit/PDFs/ArgusPdf.h"'),
    ('["<]BifurGaussPdf.hh?[">]', '"goofit/PDFs/BifurGaussPdf.h"'),
    ('["<]BinTransformPdf.hh?[">]', '"goofit/PDFs/BinTransformPdf.h"'),
    ('["<]BWPdf.hh?[">]', '"goofit/PDFs/BWPdf.h"'),
    ('["<]CompositePdf.hh?[">]', '"goofit/PDFs/CompositePdf.h"'),
    ('["<]ConvolutionPdf.hh?[">]', '"goofit/PDFs/ConvolutionPdf.h"'),
    ('["<]CorrGaussianPdf.hh?[">]', '"goofit/PDFs/CorrGaussianPdf.h"'),
    ('["<]CrystalBallPdf.hh?[">]', '"goofit/PDFs/CrystalBallPdf.h"'),
    ('["<]DalitzPlotHelpers.hh?[">]', '"goofit/PDFs/DalitzPlotHelpers.h"'),
    ('["<]DalitzPlotPdf.hh?[">]', '"goofit/PDFs/DalitzPlotPdf.h"'),
    ('["<]DalitzVetoPdf.hh?[">]', '"goofit/PDFs/DalitzVetoPdf.h"'),
    ('["<]devcomplex.hh?[">]', '"goofit/PDFs/devcomplex.h"'),
    ('["<]DP4Pdf.hh?[">]', '"goofit/PDFs/DP4Pdf.h"'),
    ('["<]EvalVar.hh?[">]', '"goofit/PDFs/EvalVar.h"'),
    ('["<]EventWeightedAddPdf.hh?[">]', '"goofit/PDFs/EventWeightedAddPdf.h"'),
    ('["<]ExpGausPdf.hh?[">]', '"goofit/PDFs/ExpGausPdf.h"'),
    ('["<]ExpPdf.hh?[">]', '"goofit/PDFs/ExpPdf.h"'),
    ('["<]GaussianPdf.hh?[">]', '"goofit/PDFs/GaussianPdf.h"'),
    ('["<]GooPdf.hh?[">]', '"goofit/PDFs/GooPdf.h"'),
    ('["<]IncoherentSumPdf.hh?[">]', '"goofit/PDFs/IncoherentSumPdf.h"'),
    ('["<]InterHistPdf.hh?[">]', '"goofit/PDFs/InterHistPdf.h"'),
    ('["<]JohnsonSUPdf.hh?[">]', '"goofit/PDFs/JohnsonSUPdf.h"'),
    ('["<]KinLimitBWPdf.hh?[">]', '"goofit/PDFs/KinLimitBWPdf.h"'),
    ('["<]LandauPdf.hh?[">]', '"goofit/PDFs/LandauPdf.h"'),
    ('["<]LineshapesPdf.hh?[">]', '"goofit/PDFs/LineshapesPdf.h"'),
    ('["<]MappedPdf.hh?[">]', '"goofit/PDFs/MappedPdf.h"'),
    ('["<]MixingTimeResolution_Aux.hh?[">]', '"goofit/PDFs/MixingTimeResolution_Aux.h"'),
    ('["<]NovosibirskPdf.hh?[">]', '"goofit/PDFs/NovosibirskPdf.h"'),
    ('["<]PolynomialPdf.hh?[">]', '"goofit/PDFs/PolynomialPdf.h"'),
    ('["<]ProdPdf.hh?[">]', '"goofit/PDFs/ProdPdf.h"'),
    ('["<]ResonancePdf.hh?[">]', '"goofit/PDFs/ResonancePdf.h"'),
    ('["<]ScaledGaussianPdf.hh?[">]', '"goofit/PDFs/ScaledGaussianPdf.h"'),
    ('["<]SmoothHistogramPdf.hh?[">]', '"goofit/PDFs/SmoothHistogramPdf.h"'),
    ('["<]SpinFactors.hh?[">]', '"goofit/PDFs/SpinFactors.h"'),
    ('["<]SpinHelper.hh?[">]', '"goofit/PDFs/SpinHelper.h"'),
    ('["<]StepPdf.hh?[">]', '"goofit/PDFs/StepPdf.h"'),
    ('["<]Tddp4Pdf.hh?[">]', '"goofit/PDFs/Tddp4Pdf.h"'),
    ('["<]TddpPdf.hh?[">]', '"goofit/PDFs/TddpPdf.h"'),
    ('["<]ThreeGaussResolution_Aux.hh?[">]', '"goofit/PDFs/ThreeGaussResolution_Aux.h"'),
    ('["<]TrigThresholdPdf.hh?[">]', '"goofit/PDFs/TrigThresholdPdf.h"'),
    ('["<]TruthResolution_Aux.hh?[">]', '"goofit/PDFs/TruthResolution_Aux.h"'),
    ('["<]VoigtianPdf.hh?[">]', '"goofit/PDFs/VoigtianPdf.h"'),
    ('["<]fakeTH1F.hh?[">]', '"goofit/rootstuff/fakeTH1F.h"'),
    ('["<]TMinuit.hh?[">]', '"goofit/rootstuff/TMinuit.h"'),
    ('["<]TRandom.hh?[">]', '"goofit/rootstuff/TRandom.h"'),
    ('["<]TRandom3.hh?[">]', '"goofit/rootstuff/TRandom3.h"'),
]

def fix_text(contents):
    """
    >>> text = '''
    ... #include "Variable.h"
    ... #include "GlobalCudaDefines.h"
    ... #include "goofit/FitControl.h"
    ... #include <set>
    ... #include <goofit/BinnedDataSet.h>
    ... #include <UnbinnedDataSet.h>
    ... '''

    >>> corrected = fix_text(text)
      Converting ["<]GlobalCudaDefines.hh?[">] -> "goofit/GlobalCudaDefines.h"
      Converting ["<]UnbinnedDataSet.hh?[">] -> "goofit/UnbinnedDataSet.h"
      Converting ["<]Variable.hh?[">] -> "goofit/Variable.h"
    >>> print(corrected.strip())
    #include "goofit/Variable.h"
    #include "goofit/GlobalCudaDefines.h"
    #include "goofit/FitControl.h"
    #include <set>
    #include <goofit/BinnedDataSet.h>
    #include "goofit/UnbinnedDataSet.h"
    """

    for r in conversion:
        after = r[1]
        before = re.compile(r[0])
        if before.search(contents):
            print('  Converting {0} -> {1}'.format(r[0], after))
            contents = before.sub(after,contents)
    return contents

def fix_files(src):
    for name in src:
        print('Converting: {0}'.format(name))
        with name.open('r') as f:
            contents = f.read()
        contents = fix_text(contents)
        with name.open('w') as f:
            f.write(contents)



class ModernizeGooFit(cli.Application):

    @cli.positional(cli.ExistingFile)
    def main(self, *src):
        fix_files(src)


if __name__ == '__main__':
    ModernizeGooFit()

