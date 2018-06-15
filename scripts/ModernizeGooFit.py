#!/usr/bin/env python

from __future__ import print_function
import re
try:
    from plumbum import local, cli
except ImportError:
    print("This file uses the plumbum library. Install with pip or conda (user directory or virtual environment OK).")
    raise

DIR = local.path(__file__).dirname

conversion = dict()
conversion['brackets'] = [
    (r'^#include\s+"(\S+)"', r'#include <\1>'),
]
conversion['2.0'] = [
    #(r'^#include\s+["<]cuda_runtime_api.hh?[">]\s*$', r'// Fake cuda has been removed (cuda_runtime_api.h requested)'),
    #(r'^#include\s+["<]driver_types.hh?[">]\s*$', r'// Fake cuda has been removed (driver_types.h requested)'),
    #(r'^#include\s+["<]host_defines.hh?[">]\s*$', r'// Fake cuda has been removed (host_defines.h requested)'),
    ('["<]Application.hh?[">]', '<goofit/Application.h>'),
    ('["<]BinnedDataSet.hh?[">]', '<goofit/BinnedDataSet.h>'),
    ('["<]DataSet.hh?[">]', '<goofit/DataSet.h>'),
    ('["<]Faddeeva.hh?[">]', '<goofit/Faddeeva.h>'),
    ('["<]FitControl.hh?[">]', '<goofit/FitControl.h>'),
    ('["<]FitManager.hh?[">]', '<goofit/FitManager.h>'),
    ('["<]FitManagerMinuit1.hh?[">]', '<goofit/fitting/FitManagerMinuit1.h>'),
    ('["<]FitManagerMinuit2.hh?[">]', '<goofit/fitting/FitManagerMinuit2.h>'),
    #(r'^#include\ ["<]FitManagerMinuit3.hh?[">]', r'// Fit Manager 3 removed'),
    ('["<]FunctorWriter.hh?[">]', '<goofit/FunctorWriter.h>'),
    ('["<]GlobalCudaDefines.hh?[">]', '<goofit/GlobalCudaDefines.h>'),
    ('["<]PdfBase.hh?[">]', '<goofit/PdfBase.h>'),
    ('["<]UnbinnedDataSet.hh?[">]', '<goofit/UnbinnedDataSet.h>'),
    ('["<]Variable.hh?[">]', '<goofit/Variable.h>'),
    ('["<]AddPdf.hh?[">]', '<goofit/PDFs/AddPdf.h>'),
    ('["<]ArgusPdf.hh?[">]', '<goofit/PDFs/ArgusPdf.h>'),
    ('["<]BifurGaussPdf.hh?[">]', '<goofit/PDFs/BifurGaussPdf.h>'),
    ('["<]BinTransformPdf.hh?[">]', '<goofit/PDFs/BinTransformPdf.h>'),
    ('["<]BWPdf.hh?[">]', '<goofit/PDFs/BWPdf.h>'),
    ('["<]CompositePdf.hh?[">]', '<goofit/PDFs/CompositePdf.h>'),
    ('["<]ConvolutionPdf.hh?[">]', '<goofit/PDFs/ConvolutionPdf.h>'),
    ('["<]CorrGaussianPdf.hh?[">]', '<goofit/PDFs/CorrGaussianPdf.h>'),
    ('["<]CrystalBallPdf.hh?[">]', '<goofit/PDFs/CrystalBallPdf.h>'),
    ('["<]DalitzPlotHelpers.hh?[">]', '<goofit/PDFs/DalitzPlotHelpers.h>'),
    ('["<]DalitzPlotPdf.hh?[">]', '<goofit/PDFs/DalitzPlotPdf.h>'),
    ('["<]DalitzVetoPdf.hh?[">]', '<goofit/PDFs/DalitzVetoPdf.h>'),
    ('["<]devcomplex.hh?[">]', '<goofit/PDFs/devcomplex.h>'),
    ('["<]DP4Pdf.hh?[">]', '<goofit/PDFs/DP4Pdf.h>'),
    ('["<]detail/EvalVar.hh?[">]', '<goofit/PDFs/detail/EvalVar.h>'),
    ('["<]EventWeightedAddPdf.hh?[">]', '<goofit/PDFs/EventWeightedAddPdf.h>'),
    ('["<]ExpGausPdf.hh?[">]', '<goofit/PDFs/ExpGausPdf.h>'),
    ('["<]ExpPdf.hh?[">]', '<goofit/PDFs/ExpPdf.h>'),
    ('["<]GaussianPdf.hh?[">]', '<goofit/PDFs/GaussianPdf.h>'),
    ('["<]GooPdf.hh?[">]', '<goofit/PDFs/GooPdf.h>'),
    ('["<]IncoherentSumPdf.hh?[">]', '<goofit/PDFs/IncoherentSumPdf.h>'),
    ('["<]InterHistPdf.hh?[">]', '<goofit/PDFs/InterHistPdf.h>'),
    ('["<]JohnsonSUPdf.hh?[">]', '<goofit/PDFs/JohnsonSUPdf.h>'),
    ('["<]KinLimitBWPdf.hh?[">]', '<goofit/PDFs/KinLimitBWPdf.h>'),
    ('["<]LandauPdf.hh?[">]', '<goofit/PDFs/LandauPdf.h>'),
    ('["<]LineshapesPdf.hh?[">]', '<goofit/PDFs/LineshapesPdf.h>'),
    ('["<]MappedPdf.hh?[">]', '<goofit/PDFs/MappedPdf.h>'),
    ('["<]MixingTimeResolution.hh?[">]', '<goofit/PDFs/MixingTimeResolution.h>'),
    ('["<]NovosibirskPdf.hh?[">]', '<goofit/PDFs/NovosibirskPdf.h>'),
    ('["<]PolynomialPdf.hh?[">]', '<goofit/PDFs/PolynomialPdf.h>'),
    ('["<]ProdPdf.hh?[">]', '<goofit/PDFs/ProdPdf.h>'),
    ('["<]ResonancePdf.hh?[">]', '<goofit/PDFs/ResonancePdf.h>'),
    ('["<]ScaledGaussianPdf.hh?[">]', '<goofit/PDFs/ScaledGaussianPdf.h>'),
    ('["<]SmoothHistogramPdf.hh?[">]', '<goofit/PDFs/SmoothHistogramPdf.h>'),
    ('["<]SpinFactors.hh?[">]', '<goofit/PDFs/SpinFactors.h>'),
    ('["<]detail/SpinHelper.hh?[">]', '<goofit/PDFs/detail/SpinHelper.h>'),
    ('["<]StepPdf.hh?[">]', '<goofit/PDFs/StepPdf.h>'),
    ('["<]Tddp4Pdf.hh?[">]', '<goofit/PDFs/Tddp4Pdf.h>'),
    ('["<]TddpPdf.hh?[">]', '<goofit/PDFs/TddpPdf.h>'),
    ('["<]ThreeGaussResolution.hh?[">]', '<goofit/PDFs/ThreeGaussResolution.h>'),
    ('["<]TrigThresholdPdf.hh?[">]', '<goofit/PDFs/TrigThresholdPdf.h>'),
    ('["<]TruthResolution.hh?[">]', '<goofit/PDFs/TruthResolution.h>'),
    ('["<]VoigtianPdf.hh?[">]', '<goofit/PDFs/VoigtianPdf.h>'),
    ('["<]fakeTH1F.hh?[">]', '<TH1F.h>'),
    (r'\bALIGN\b', '__align__'),
    (r'\bMEM_CONSTANT\b', '__constant__'),
    (r'\bMEM_DEVICE\b', '__device__'),
    (r'\bEXEC_TARGET\b', '__device__'),
    (r'\bMEM_SHARED\b', '__shared__'),
    (r'\bSYNCH\b', 'cudaDeviceSynchronize'),
    (r'\bdevcomplex\b', 'thrust::complex'),
    (r'\bstd::complex\b', 'thrust::complex'),
    (r'\bDEVICE_VECTOR\b', 'thrust::device_vector'),
    (r'\bATAN2\b', 'atan2'),
    (r'\bCOS\b', 'cos'),
    (r'\bCOSH\b', 'cosh'),
    (r'\bSINH\b', 'sinh'),
    (r'\bERF\b', 'erf'),
    (r'\bERFC\b', 'erfc'),
    (r'\bEXP\b', 'exp'),
    (r'\bFABS\b', 'fabs'),
    (r'\bFMOD\b', 'fmod'),
    (r'\bLOG\b', 'log'),
    (r'\bMODF\b', 'modf'),
    (r'\bSIN\b', 'sin'),
    (r'\bSQRT\b', 'sqrt'),
    (r'\bRSQRT\b', 'rsqrt'),
    (r'\bFLOOR\b', 'floor'),
    (r'\bCONST_PI\b', 'M_PI'),
    #(r'normalise', 'normalize'),
    (r'Normalise', 'Normalize'),
    (r'initialise', 'initialize'),
    (r'Initialize', 'Initialize'),
    (r'\bPdfBase::parCont\b', 'std::vector<Variable*>'),
    (r'\bPdfBase::obsCont\b', 'std::vector<Variable*>'),
    (r'\bobsCont\b', 'std::vector<Variable*>'),
    (r'\bparCont\b', 'std::vector<Variable*>'),
    (r'\babortWithCudaPrintFlush\b', 'GooFit::abort'),
    (r'goofit/PDFs/ArgusPdf.h', 'goofit/PDFs/basic/ArgusPdf.h'),
    (r'goofit/PDFs/BifurGaussPdf.h', 'goofit/PDFs/basic/BifurGaussPdf.h'),
    (r'goofit/PDFs/BWPdf.h', 'goofit/PDFs/basic/BWPdf.h'),
    (r'goofit/PDFs/BinTransformPdf.h', 'goofit/PDFs/basic/BinTransformPdf.h'),
    (r'goofit/PDFs/CorrGaussianPdf.h', 'goofit/PDFs/basic/CorrGaussianPdf.h'),
    (r'goofit/PDFs/CrystalBallPdf.h', 'goofit/PDFs/basic/CrystalBallPdf.h'),
    (r'goofit/PDFs/ExpGausPdf.h', 'goofit/PDFs/basic/ExpGausPdf.h'),
    (r'goofit/PDFs/ExpPdf.h', 'goofit/PDFs/basic/ExpPdf.h'),
    (r'goofit/PDFs/GaussianPdf.h', 'goofit/PDFs/basic/GaussianPdf.h'),
    (r'goofit/PDFs/InterHistPdf.h', 'goofit/PDFs/basic/InterHistPdf.h'),
    (r'goofit/PDFs/JohnsonSUPdf.h', 'goofit/PDFs/basic/JohnsonSUPdf.h'),
    (r'goofit/PDFs/KinLimitBWPdf.h', 'goofit/PDFs/basic/KinLimitBWPdf.h'),
    (r'goofit/PDFs/LandauPdf.h', 'goofit/PDFs/basic/LandauPdf.h'),
    (r'goofit/PDFs/NovosibirskPdf.h', 'goofit/PDFs/basic/NovosibirskPdf.h'),
    (r'goofit/PDFs/PolynomialPdf.h', 'goofit/PDFs/basic/PolynomialPdf.h'),
    (r'goofit/PDFs/ScaledGaussianPdf.h', 'goofit/PDFs/basic/ScaledGaussianPdf.h'),
    (r'goofit/PDFs/SmoothHistogramPdf.h', 'goofit/PDFs/basic/SmoothHistogramPdf.h'),
    (r'goofit/PDFs/StepPdf.h', 'goofit/PDFs/basic/StepPdf.h'),
    (r'goofit/PDFs/VoigtianPdf.h', 'goofit/PDFs/basic/VoigtianPdf.h'),
    (r'goofit/PDFs/TrigThresholdPdf.h', 'goofit/PDFs/basic/TrigThresholdPdf.h'),
    (r'goofit/PDFs/AddPdf.h', 'goofit/PDFs/combine/AddPdf.h'),
    (r'goofit/PDFs/CompositePdf.h', 'goofit/PDFs/combine/CompositePdf.h'),
    (r'goofit/PDFs/ConvolutionPdf.h', 'goofit/PDFs/combine/ConvolutionPdf.h'),
    (r'goofit/PDFs/EventWeightedAddPdf.h', 'goofit/PDFs/combine/EventWeightedAddPdf.h'),
    (r'goofit/PDFs/MappedPdf.h', 'goofit/PDFs/combine/MappedPdf.h'),
    (r'goofit/PDFs/ProdPdf.h', 'goofit/PDFs/combine/ProdPdf.h'),
    (r'goofit/PDFs/DalitzPlotHelpers.h', 'goofit/PDFs/physics/DalitzPlotHelpers.h'),
    (r'goofit/PDFs/DalitzPlotPdf.h', 'goofit/PDFs/physics/DalitzPlotPdf.h'),
    (r'goofit/PDFs/DalitzVetoPdf.h', 'goofit/PDFs/physics/DalitzVetoPdf.h'),
    (r'goofit/PDFs/DP4Pdf.h', 'goofit/PDFs/physics/DP4Pdf.h'),
    (r'goofit/PDFs/detail/EvalVar.h', 'goofit/PDFs/physics/detail/EvalVar.h'),
    (r'goofit/PDFs/IncoherentSumPdf.h', 'goofit/PDFs/physics/IncoherentSumPdf.h'),
    (r'goofit/PDFs/LineshapesPdf.h', 'goofit/PDFs/physics/LineshapesPdf.h'),
    (r'goofit/PDFs/MixingTimeResolution.h', 'goofit/PDFs/physics/MixingTimeResolution.h'),
    (r'goofit/PDFs/ResonancePdf.h', 'goofit/PDFs/physics/ResonancePdf.h'),
    (r'goofit/PDFs/SpinFactors.h', 'goofit/PDFs/physics/SpinFactors.h'),
    (r'goofit/PDFs/detail/SpinHelper.h', 'goofit/PDFs/physics/detail/SpinHelper.h'),
    (r'goofit/PDFs/Tddp4Pdf.h', 'goofit/PDFs/physics/Tddp4Pdf.h'),
    (r'goofit/PDFs/TddpPdf.h', 'goofit/PDFs/physics/TddpPdf.h'),
    (r'goofit/PDFs/ThreeGaussResolution.h', 'goofit/PDFs/physics/ThreeGaussResolution.h'),
    (r'goofit/PDFs/TruthResolution.h', 'goofit/PDFs/physics/TruthResolution.h')
]
conversion['2.1'] = [
    (r'thrust::complex<fptype>', r'fpcomplex'),
    (r'CountingVariable', r'EventNumber'),
    (r'\bDecayInfo\b\s*\*', r'DecayInfo3'), # Might need a t if time dependent
    (r'\bDecayInfo_DP\b\s*\*', r'DecayInfo4'), # Might need a t if time dependent
    (r'Variable\s*\*\s*(\w+)\s*=\s*new\s+Variable', r'Variable \1'),
    (r'new\s+Variable', r'Variable'),
    (r'Variable(\s*)\*', r'Variable\1'),
    (r'ResonancePDF\((.*?)$', r'Resonances::RBW(\1 // Could also be LASS, GS, FLATTE, Gauss, NonRes'),
    (r'LS::BW', r'LS::RBW'),
    (r'LS::ONE', r'LS::One'),
    (r'LS::nonRes', r'LS::NonRes'),
    (r'LS::Lass', r'LS::LASS'),
    (r'LS::Lass_M3', r'LS::GLASS'), # No change to Bugg, Bugg3, Flatte, SBW
    (r'''Variable\s+\w+\s*\(
        ([^,^)]+),
        ([^,^)]+),
        ([^,^)]+)
        \)''', r'Observable(\1,\2,\3)'),
    (r'''Lineshape\s*\(
        ([^,]+),
        ([^,^)]+),
        ([^,^)]+),
        ([^,^)]+),
        ([^,^)]+),
        \s*LS::(\w+)\s*,
        ([^)]+)     \)
    ''', r'Lineshapes::\6(\1,\2,\3,\4,\5,\7)'),
    (r'''Lineshape\s*\(
        ([^,]+),
        (\s*Variable\([^)]+\)),
        (\s*Variable\([^)]+\)),
        ([^,^)]+),
        ([^,^)]+),
        \s*LS::(\w+)\s*,
        ([^)]+)     \)
    ''', r'Lineshapes::\6(\1,\2,\3,\4,\5,\7)'),
    (r'vector<Variable>\s*observable', r'vector<Observable> observable'), # Common pattern
]

conversion['2.2'] = [
    (r'DalitzPlotPdf', r'Amp3Body'),
    (r'TddpPdf', r'Amp3Body_TD'),
    (r'DPPdf', r'Amp4Body'),
    (r'DP4Pdf', r'Amp4Body'),
    (r'TDDP4', r'Amp4Body_TD'),
    (r'LineshapesPdf', r'Lineshapes'),
    (r'IncoherentSumPdf', r'Amp3Body_IS'),
]

def fix_text(contents, version='2.0'):
    r"""
    >>> text = '''
    ... #include "Variable.h"
    ... #include "GlobalCudaDefines.h"
    ... #include "goofit/FitControl.h"
    ... #include <set>
    ... #include <goofit/BinnedDataSet.h>
    ... #include <UnbinnedDataSet.h>
    ... #include <FitManagerMinuit3.hh>
    ... SYNCH();
    ... THREAD_SYNCH();
    ... '''

    >>> corrected = fix_text(text)
      Converting #include ["<]FitManagerMinuit3.hh?[">] -> \\\\ Fit Manager 3 removed
      Converting ["<]GlobalCudaDefines.hh?[">] -> "goofit/GlobalCudaDefines.h"
      Converting ["<]UnbinnedDataSet.hh?[">] -> "goofit/UnbinnedDataSet.h"
      Converting ["<]Variable.hh?[">] -> "goofit/Variable.h"
      Converting \bSYNCH\b -> cudaDeviceSynchronize
    >>> print(corrected.strip())
    #include "goofit/Variable.h"
    #include "goofit/GlobalCudaDefines.h"
    #include "goofit/FitControl.h"
    #include <set>
    #include <goofit/BinnedDataSet.h>
    #include "goofit/UnbinnedDataSet.h"
    \\ Fit Manager 3 removed
    cudaDeviceSynchronize();
    THREAD_SYNCH();
    """

    for r in conversion[version]:
        after = r[1]
        before = re.compile(r[0], re.MULTILINE | re.VERBOSE)
        if before.search(contents):
            print('  Converting {0} -> {1}'.format(r[0], after))
            contents = before.sub(after,contents)
    return contents

def fix_files(src, version):
    for name in src:
        with name.open('r') as f:
            contents = f.read()
        new_contents = fix_text(contents, version)
        if contents != new_contents:
            print('Converted: {0}'.format(name))
            with name.open('w') as f:
                f.write(new_contents)



class ModernizeGooFit(cli.Application):
    set_version = cli.SwitchAttr(['-v','--version'], cli.Set(*conversion), default='2.0', help='The version to convert')
    source = cli.Flag(['--source'], help="Run on the GooFit Source")

    @cli.positional(cli.ExistingFile)
    def main(self, *src):
        if not src:
            assert self.source, "You must use the --source flag to run over GooFit's source"
            git = local['git']
            with local.cwd(DIR / '../'):
                src = [local.path(n) for n in
                        git('ls-files', '--', '*.cpp', '*.h', '*.cu').splitlines()]
        fix_files(src, self.set_version)



if __name__ == '__main__':
    ModernizeGooFit()

