goofit_add_library(_Physics STATIC
    DalitzPlotPdf.cpp
    DalitzVetoPdf.cpp
    DalitzPlotHelpers.cpp
    DalitzPlotter.cpp
    DP4Pdf.cu
    IncoherentSumPdf.cpp
    LineshapesPdf.cpp
    MixingTimeResolution_Aux.cpp
    ResonancePdf.cpp
    SpinFactors.cpp
    Tddp4Pdf.cu
    TddpPdf.cpp
    ThreeGaussResolution_Aux.cpp
    TruthResolution_Aux.cpp
)
target_link_libraries(_Physics PRIVATE _goofit_python)
set_target_properties(_Physics PROPERTIES FOLDER python/PDFs)

goofit_add_python_docs(_Physics PDFs/physics/DalitzPlotPdf.h)
goofit_add_python_docs(_Physics PDFs/physics/DalitzVetoPdf.h)
goofit_add_python_docs(_Physics PDFs/physics/IncoherentSumPdf.h)
goofit_add_python_docs(_Physics PDFs/physics/MixingTimeResolution_Aux.h)
goofit_add_python_docs(_Physics PDFs/physics/ResonancePdf.h)
goofit_add_python_docs(_Physics PDFs/physics/TddpPdf.h)
goofit_add_python_docs(_Physics PDFs/physics/ThreeGaussResolution_Aux.h)
goofit_add_python_docs(_Physics PDFs/physics/TruthResolution_Aux.h)
