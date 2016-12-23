
set(SRC
    TMinuit.cc
    TRandom.cc
    TRandom3.cc
)

add_library(rootstuff STATIC ${SRC})
target_include_directories(rootstuff PUBLIC ${PROJECT_SOURCE_DIR}/include/goofit/rootstuff)
export(TARGETS rootstuff APPEND FILE ${PROJECT_BINARY_DIR}/GooFitConfig.cmake)
