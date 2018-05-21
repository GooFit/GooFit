
# Warning suppression only needed until next release
if(NOT DEFINED CMAKE_SUPPRESS_DEVELOPER_WARNINGS)
    set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE BOOL "")
    add_subdirectory("${PROJECT_SOURCE_DIR}/extern/googletest" "extern/googletest")
    unset(CMAKE_SUPPRESS_DEVELOPER_WARNINGS)
else()
    add_subdirectory("${PROJECT_SOURCE_DIR}/extern/googletest" "extern/googletest")
endif()


mark_as_advanced(
    BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
    gmock_build_tests gtest_build_samples gtest_build_tests
    gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)

set_target_properties(gtest PROPERTIES FOLDER extern)
set_target_properties(gtest_main PROPERTIES FOLDER extern)
set_target_properties(gmock PROPERTIES FOLDER extern)
set_target_properties(gmock_main PROPERTIES FOLDER extern)

macro(goofit_add_test TESTNAME)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TESTNAME}.cu")
        goofit_add_executable(${TESTNAME} ${TESTNAME}.cu ${ARGN})
    else()
        goofit_add_executable(${TESTNAME} ${TESTNAME}.cpp ${ARGN})
    endif()
    target_link_libraries(${TESTNAME} PUBLIC gtest gmock gtest_main)
    add_test(${TESTNAME} ${TESTNAME})
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
    set_property(TEST ${TESTNAME} PROPERTY COST 1000)
endmacro()

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
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
    set_property(TEST ${TESTNAME} PROPERTY COST 1000)
endmacro()

add_subdirectory(simple)
add_subdirectory(PDFs)
add_subdirectory(convert)

