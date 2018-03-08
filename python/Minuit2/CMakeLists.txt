
add_library(minuit2 MODULE
    main.cpp
    FCNBase.cpp
    FunctionMinimum.cpp
    MinuitParameter.cpp
    MnApplication.cpp
    MnMigrad.cpp
    MnPrint.cpp
    MnUserCovariance.cpp
    MnUserParameterState.cpp
    MnUserParameters.cpp
)

process_pybind_module(minuit2)
target_link_libraries(minuit2 PUBLIC ROOT::Minuit2)
