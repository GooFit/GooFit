.PHONY: default auto xcode omp cuda info mpi docs cmake base

UNAME := $(shell uname)
CMAKE_VER := 3.11.2
CMAKE_BASE := $(basename $(CMAKE_VER))
ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
CMAKE_DIR := cmake-$(CMAKE_VER)

HAS_WGET := $(shell command -v wget 2> /dev/null)
ifdef HAS_WGET
    GET := wget -qO-
else
    GET := curl -s
endif

HAS_DIR := $(shell if [ -d $(CMAKE_DIR) ] ; then echo "yes" ; else echo "no" ; fi)

ifeq ($(HAS_DIR), yes)
    CMAKE_EXE := $(ROOT_DIR)/$(CMAKE_DIR)/bin/cmake
    CTEST_EXE := $(ROOT_DIR)/$(CMAKE_DIR)/bin/ctest
else
    CMAKE_EXE := cmake
    CTEST_EXE := ctest
endif

default: info

base:
	@echo "Using: $(CMAKE_EXE)"

auto: base
	@mkdir -p build
	cd build && $(CMAKE_EXE) .. && $(MAKE) --no-print-directory
	cd build && $(CTEST_EXE) --output-on-failure

omp: base
	@mkdir -p build-omp
	cd build-omp && $(CMAKE_EXE) .. -DGOOFIT_DEVICE=OMP && $(MAKE) --no-print-directory
	cd build-omp && $(CTEST_EXE) --output-on-failure

cuda: base
	@mkdir -p build-cuda
	cd build-cuda && $(CMAKE_EXE) .. -DGOOFIT_DEVICE=CUDA && $(MAKE) --no-print-directory
	cd build-cuda && $(CTEST_EXE) --output-on-failure

mpi: base
	@mkdir -p build-mpi
	cd build-mpi && $(CMAKE_EXE) .. -DGOOFIT_MPI=ON && $(MAKE) --no-print-directory
	cd build-mpi && $(CTEST_EXE) --output-on-failure

docs:
	cd docs && doxygen Doxyfile

clang-format:
	git ls-files -- '*.cu' '*.cc' '*.h' '*.cpp' | xargs clang-format -i -style=file

ifeq ($(UNAME), Darwin)

xcode: base
	@mkdir -p build-xcode
	cd build-xcode && cmake .. -GXcode
	open build-xcode/GOOFIT.xcodeproj

cmake:
	@echo "CMake download on macOS not supported. Please use:"
	@echo "  brew install cmake"
	@echo "to keep your CMake install up to date."

else

cmake: $(CMAKE_DIR)/bin/cmake
	@echo "Run:"
	@echo '  export PATH=$(ROOT_DIR)/cmake-$(CMAKE_VER)/bin:$$PATH'
	@echo "Note: this makefile will use the new cmake from now on regardless of path settings"

$(CMAKE_DIR):
	mkdir -p cmake-$(CMAKE_VER)

$(CMAKE_DIR)/bin/cmake: | $(CMAKE_DIR)
	$(GET) "https://cmake.org/files/v$(CMAKE_BASE)/cmake-$(CMAKE_VER)-$(UNAME)-x86_64.tar.gz" | tar --strip-components=1 -xz -C cmake-$(CMAKE_VER)

endif

info:
	@echo "Shortcuts for the laziest of GooUsers:"
	@echo "  auto  - Create a  best-available build in the build       directory"
	@echo "  omp   - Create an     OpenMP     build in the build-omp   directory"
	@echo "  cuda  - Create a       CUDA      build in the build-cuda  directory"
	@echo "  cuda  - Create a        MPI      build in the build-mpi   directory"
	@echo "  xcode - Create an  Xcode (macOS) build in the build-xcode directory"
	@echo ""
	@echo "Other functionalities:"
	@echo "  cmake        - download CMake (on Linux systems)"
	@echo "  docs         - Build the documentation in docs/html/index.html"
	@echo "  clang-format - Run clang-format on the source code (modify inplace)"
	@echo ""
	@echo "Using the cmake target will cause subsequent make commands to use the"
	@echo "local version of CMake instead of the system default."
	@echo ""
	@echo "The recommended way to build this package is directly with CMake 3.4+"
	@echo "Using 'make auto' is equivalent to the standard CMake procedure:"
	@echo "  mkdir build"
	@echo "  cd build"
	@echo "  cmake .."
	@echo "  make"
	@echo ""
