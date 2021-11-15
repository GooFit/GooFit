target_sources(_goofit PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/VariableBinTransform1DPdf.cpp")
set_target_properties(_goofit PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)