if(ROOT_FOUND)
    add_subdirectory(2d_plot)
    add_subdirectory(addition)
    add_subdirectory(chisquare)
    add_subdirectory(convolution)
    add_subdirectory(dalitz)
    add_subdirectory(DP4)
    add_subdirectory(pipipi0DPFit)
    add_subdirectory(product)
    add_subdirectory(SigGen)
    add_subdirectory(simpleFit)
    add_subdirectory(TDDP4)
    add_subdirectory(zachFit)
endif()

# Root not required
add_subdirectory(exponential)

goofit_add_link(RunAll.py)
