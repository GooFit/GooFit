# This is only needed for IDE header display
file(GLOB_RECURSE GOOFIT_PHYSICS_HEADERS "${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/*.h")

if(GOOFIT_KMATRIX)
  set(GOOFIT_KMATRIX_SOURCE lineshapes/kMatrix.cu resonances/kMatrix.cu kMatrixUtils.cu)
else()
  set(GOOFIT_KMATRIX_SOURCE)
endif()

goofit_add_pdf_library(
  PDFPhysics
  Amp3Body.cu
  Amp3BodyBase.cu
  Amp3Body_TD.cu
  Amp3Body_IS.cu
  Amp4Body.cu
  Amp4BodyBase.cu
  Amp4Body_TD.cu
  Amp4BodyGlobals.cu
  Amplitude.cu
  DalitzPlotHelpers.cu
  DalitzVetoPdf.cu
  SpinFactors.cu
  MixingTimeResolution.cu
  ThreeGaussResolution.cu
  TruthResolution.cu
  SquareDalitzEffPdf.cu
  detail/AmpCalc.cu
  detail/AmpCalc_TD.cu
  detail/LSCalculator.cu
  detail/LSCalculator_TD.cu
  detail/NormIntegrator.cu
  detail/NormIntegrator_TD.cu
  detail/NormLSCalculator.cu
  detail/NormLSCalculator_TD.cu
  detail/NormSpinCalculator.cu
  detail/NormSpinCalculator_TD.cu
  detail/SFCalculator.cu
  detail/SFCalculator_TD.cu
  detail/FourDblTupleAdd.cu
  detail/SpecialIncoherentIntegrator.cu
  detail/SpecialIncoherentResonanceCalculator.cu
  detail/SpecialResonanceIntegrator.cu
  detail/SpecialResonanceCalculator.cu
  detail/SpecialDalitzIntegrator.cu
  detail/SpecialWaveCalculator.cu
  lineshapes/Lineshape.cu
  lineshapes/BW.cu
  lineshapes/Common.cu
  lineshapes/FOCUS.cu
  lineshapes/Bugg.cu
  lineshapes/Flatte.cu
  lineshapes/GLASS.cu
  lineshapes/GSpline.cu
  lineshapes/NonRes.cu
  lineshapes/One.cu
  ${GOOFIT_KMATRIX_SOURCE}
  resonances/Resonance.cu
  resonances/FLATTE.cu
  resonances/Gauss.cu
  resonances/GS.cu
  resonances/LASS.cu
  resonances/NonRes.cu
  resonances/RBW.cu
  resonances/RhoOmegaMix.cu
  resonances/Spline.cu
  ${GOOFIT_PHYSICS_HEADERS}
  detail/NormEvents_4Body_Base.cu
  detail/NormEvents_4Body_DeviceCached.cu
  detail/NormEvents_4Body_HostCached.cu)

if(IS_NOT_CUDA
   AND UNIX
   AND NOT APPLE)
  target_link_libraries(PDFPhysics PRIVATE rt)
endif()
