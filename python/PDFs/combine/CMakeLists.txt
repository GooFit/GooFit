add_library(_Combine STATIC AddPdf.cpp CompositePdf.cpp ConvolutionPdf.cpp EventWeightedAddPdf.cpp
                            MappedPdf.cpp ProdPdf.cpp)

target_link_libraries(_Combine PRIVATE _goofit_python)
set_target_properties(_Combine PROPERTIES FOLDER python/PDFs CUDA_RESOLVE_DEVICE_SYMBOLS ON)

goofit_add_python_docs(_Combine PDFs/combine/AddPdf.h)
goofit_add_python_docs(_Combine PDFs/combine/CompositePdf.h)
goofit_add_python_docs(_Combine PDFs/combine/ConvolutionPdf.h)
goofit_add_python_docs(_Combine PDFs/combine/EventWeightedAddPdf.h)
goofit_add_python_docs(_Combine PDFs/combine/MappedPdf.h)
goofit_add_python_docs(_Combine PDFs/combine/ProdPdf.h)
