
add_library(catch_main catch_main.cpp)
target_link_libraries(catch_main PUBLIC Catch2::Catch)
set_target_properties(catch_main PROPERTIES FOLDER tests)

macro(goofit_catch_test TESTNAME)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TESTNAME}.cu")
        goofit_add_executable(${TESTNAME} ${TESTNAME}.cu "${PROJECT_SOURCE_DIR}/tests/include/goofit/Catch.h" ${ARGN})
    else()
        goofit_add_executable(${TESTNAME} ${TESTNAME}.cpp "${PROJECT_SOURCE_DIR}/tests/include/goofit/Catch.h" ${ARGN})
    endif()
    target_link_libraries(${TESTNAME} PUBLIC catch_main)
    target_include_directories(${TESTNAME} PUBLIC "${PROJECT_SOURCE_DIR}/tests/include")
    add_test(${TESTNAME} ${TESTNAME})

    get_filename_component(CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    set_target_properties(${TESTNAME} PROPERTIES FOLDER "tests/${CURRENT_DIR_NAME}")
    set_property(TEST ${TESTNAME} PROPERTY COST 1000)
endmacro()

add_subdirectory(simple)
add_subdirectory(PDFs)
add_subdirectory(convert)

