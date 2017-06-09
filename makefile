.PHONY: default auto xcode omp cuda warning mpi docs

default: auto warning

auto: 
	@mkdir -p build
	cd build && cmake .. && $(MAKE) --no-print-directory
	cd build && ctest --output-on-failure

omp:
	@mkdir -p build-omp
	cd build-omp && cmake .. -DGOOFIT_DEVICE=OMP && $(MAKE) --no-print-directory
	cd build-omp && ctest --output-on-failure

cuda:
	@mkdir -p build-cuda
	cd build-cuda && cmake .. -DGOOFIT_DEVICE=CUDA && $(MAKE) --no-print-directory
	cd build-cuda && ctest --output-on-failure

xcode:
	@mkdir -p xbuild
	cd xbuild && cmake .. -GXcode
	open xbuild/GOOFIT.xcodeproj

mpi:
	@mkdir -p build-mpi
	cd build-mpi && cmake .. -DGOOFIT_MPI=ON && $(MAKE) --no-print-directory
	cd build-mpi && ctest --output-on-failure

docs:
	cd docs && doxygen Doxyfile

clang-format:
	git ls-files -- '*.cu' '*.cc' '*.h' '*.cpp' | xargs clang-format -i -style=file

warning:
	@echo "This project builds with CMake 3.4+."
	@echo "This makefile is just a shortcut to prepare your build with CMake."
	@echo "This is roughly equivelent to the standard CMake procedure::"
	@echo "  mkdir build"
	@echo "  cd build"
	@echo "  cmake .."
	@echo "  make"
	@echo ""
	@echo "You can use 'make omp cuda' to setup or rebuild the build-cuda and build-omp directories."
	@echo "On a Mac, you can also prepare the Xcode build with make xcode"
	@echo ""
