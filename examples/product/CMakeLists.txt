goofit_add_directory()
goofit_add_executable(product product.cu)
set_target_properties(product PROPERTIES FOLDER examples)

add_test(NAME product_Example COMMAND product)
