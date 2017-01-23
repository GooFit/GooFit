
set(SRC
    TMinuit.cc
    TRandom.cc
    TRandom3.cc
)

add_library(rootstuff STATIC ${SRC})
target_include_directories(rootstuff PUBLIC ${PROJECT_SOURCE_DIR}/include/goofit/rootstuff)
install(TARGETS rootstuff EXPORT GooFitTargets DESTINATION lib)
