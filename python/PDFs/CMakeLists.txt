target_sources(_goofit PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/GooPdf.cpp")
target_sources(_goofit PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/CombinePdf.cpp")
target_sources(_goofit PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/GooPdfBkgSub.cpp")

add_subdirectory(basic)
add_subdirectory(combine)
add_subdirectory(physics)
add_subdirectory(utilities)
