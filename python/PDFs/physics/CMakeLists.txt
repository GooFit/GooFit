goofit_add_library(_Physics STATIC
    AmpNBodyBase.cpp
    Amp3BodyBase.cpp
    Amp3Body.cpp
    Amp3Body_TD.cpp
    Amp3Body_IS.cpp
    Amp4BodyBase.cpp
    Amp4Body.cpp
    Amp4Body_TD.cpp
    DalitzVetoPdf.cpp
    DalitzPlotHelpers.cpp
    DalitzPlotter.cpp
    Lineshapes.cpp
    MixingTimeResolution.cpp
    ResonancePdf.cpp
    SpinFactors.cpp
    ThreeGaussResolution.cpp
    TruthResolution.cpp
)
target_link_libraries(_Physics PRIVATE _goofit_python)
set_target_properties(_Physics PROPERTIES FOLDER python/PDFs)

goofit_add_python_docs(_Physics PDFs/physics/Amp3Body.h)
goofit_add_python_docs(_Physics PDFs/physics/Amp3Body_IS.h)
goofit_add_python_docs(_Physics PDFs/physics/Amp3Body_TD.h)
goofit_add_python_docs(_Physics PDFs/physics/DalitzVetoPdf.h)
goofit_add_python_docs(_Physics PDFs/physics/MixingTimeResolution.h)
goofit_add_python_docs(_Physics PDFs/physics/resonances/Resonance.h)
goofit_add_python_docs(_Physics PDFs/physics/ThreeGaussResolution.h)
goofit_add_python_docs(_Physics PDFs/physics/TruthResolution.h)
