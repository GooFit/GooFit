add_subdirectory(2d_plot)
add_subdirectory(addition)
add_subdirectory(chisquare)
add_subdirectory(convolution)
add_subdirectory(dalitz)
add_subdirectory(DP4)
add_subdirectory(exponential)
add_subdirectory(pipipi0DPFit)
add_subdirectory(product)
add_subdirectory(SigGen)
add_subdirectory(simpleFit)
add_subdirectory(TDDP4)
add_subdirectory(zachFit)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/int_testing.py
    COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/int_testing.py
    ${CMAKE_CURRENT_BINARY_DIR}/int_testing.py
    )

add_custom_target(
    int_testing ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/int_testing.py
    )
