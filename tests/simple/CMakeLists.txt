goofit_catch_test(VectorsTest)
goofit_catch_test(NormalizeTest)
goofit_catch_test(SimpleTest)
goofit_catch_test(BinningTest)

goofit_catch_test(BlindTest)
goofit_catch_test(MonteCarloTest)

if(ROOT_FOUND)
  goofit_catch_test(Minuit1Test)
endif()
