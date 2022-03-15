add_library(
  _Basic STATIC
  ArgusPdf.cpp
  BifurGaussPdf.cpp
  BinTransformPdf.cpp
  BernsteinPdf.cpp
  BWPdf.cpp
  CorrGaussianPdf.cpp
  CrystalBallPdf.cpp
  ExpGausPdf.cpp
  ExpPdf.cpp
  GaussianPdf.cpp
  InterHistPdf.cpp
  JohnsonSUPdf.cpp
  KinLimitBWPdf.cpp
  LandauPdf.cpp
  NovosibirskPdf.cpp
  PolynomialPdf.cpp
  ScaledGaussianPdf.cpp
  SmoothHistogramPdf.cpp
  StepPdf.cpp
  TrigThresholdPdf.cpp
  VoigtianPdf.cpp)

target_link_libraries(_Basic PRIVATE _goofit_python)
set_target_properties(_Basic PROPERTIES FOLDER python/PDFs CUDA_RESOLVE_DEVICE_SYMBOLS ON)

goofit_add_python_docs(_Basic PDFs/basic/ArgusPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/BifurGaussPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/BernsteinPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/BWPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/CorrGaussianPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/CrystalBallPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/ExpGausPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/ExpPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/GaussianPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/InterHistPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/JohnsonSUPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/KinLimitBWPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/LandauPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/NovosibirskPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/PolynomialPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/ScaledGaussianPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/SmoothHistogramPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/StepPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/VoigtianPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/TrigThresholdPdf.h)
goofit_add_python_docs(_Basic PDFs/basic/BinTransformPdf.h)
