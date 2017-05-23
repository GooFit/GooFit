.PHONY: auto omp cuda warning

auto: warning
	@mkdir -p build
	cd build && cmake .. && $(MAKE) --no-print-directory

omp:
	@mkdir -p build-omp
	cd build-omp && cmake .. -DGOOFIT_DEVICE=OMP && $(MAKE) --no-print-directory

cuda:
	@mkdir -p build-cuda
	cd build-cuda && cmake .. -DGOOFIT_DEVICE=CUDA && $(MAKE) --no-print-directory

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
