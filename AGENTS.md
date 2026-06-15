# Agent notes

GooFit is a massively-parallel C++ framework for maximum-likelihood fits (mostly high-energy physics
amplitude analyses). The same source compiles to four Thrust backends selected at build time:
`CUDA`, `OMP` (OpenMP), `TBB`, or `CPP` (single-threaded). There are also pybind11-based Python bindings.

Most "device" source files are `.cu` (CUDA syntax) but compile for every backend — Thrust abstracts the
device. Files compile under `nvcc` for CUDA and as plain C++ otherwise.

## Building

The build is CMake-based. Convenience targets in the top-level `makefile` create a build dir, configure,
build, and run ctest:

```bash
make auto   # best-available backend  -> build/
make omp    # OpenMP                  -> build-omp/
make cuda   # CUDA                    -> build-cuda/
make mpi    # MPI on                  -> build-mpi/
```

Direct CMake (preferred for control):

```bash
cmake -S . -B build -DGOOFIT_DEVICE=OMP -DGOOFIT_TESTS=ON
cmake --build build -j$(nproc)
```

Key options: `-DGOOFIT_DEVICE=` (`Auto`/`CUDA`/`OMP`/`TBB`/`CPP`), `-DGOOFIT_EXAMPLES=ON`,
`-DGOOFIT_TESTS=ON`, `-DGOOFIT_PYTHON=ON`, `-DGOOFIT_MPI=ON`, `-DGOOFIT_DEBUG=ON`/`-DGOOFIT_TRACE=ON`
(enable the matching printout macros). `Auto` picks CUDA if found, else OpenMP, else CPP. On macOS the
default Clang gives a `CPP` build; OpenMP needs `brew install libomp`.

ROOT 6 is optional but recommended — without it, the bundled Minuit2 submodule is used and the Minuit1
fitter is unavailable. Submodules under `extern/` are required: clone with `--recursive` or run
`git submodule update --init --recursive`.

## Tests

C++ tests use Catch2 (`tests/`, enabled by `-DGOOFIT_TESTS=ON`):

```bash
cd build && ctest --output-on-failure        # all tests
./tests/GooFitTest "[simple]"                 # filter by Catch2 tag
./tests/GooFitTest "Adding values*"           # filter by test name
```

Python tests (when `-DGOOFIT_PYTHON=ON`) use pytest, run from the build dir; `pytest` config lives in
`pyproject.toml` (`testpaths = ["python/tests"]`). `nox -s test` does an isolated install + pytest.

Examples: `./examples/RunAll.py` (needs `plumbum`); individual examples build into the build dir from
sources in `examples/`.

## Linting

Use `prek -a --quiet` (not `pre-commit run -a`). Hooks: black + ruff (Python), clang-format (C/C++/CUDA,
`.clang-format`), cmake-format, codespell, shellcheck. C++ style can also be applied with `make clang-format`.

## Architecture

### PDF model — the core abstraction

Every distribution is a `GooPdf` (subclass of `PdfBase`). A PDF has two halves:

- **Host side** — a C++ class (e.g. `GaussianPdf`) constructed with `Observable`s and `Variable`s. The
  constructor calls `registerFunction("ptr_to_X", ptr_to_X)` to bind its device evaluator, then
  `initialize()`.
- **Device side** — a free `__device__` function (e.g. `device_Gaussian(fptype *evt, ParameterContainer &pc)`)
  exposed through a `device_function_ptr`. This runs once per event on the device.

Parameters and observables are **not** passed as arguments. They are packed into flat arrays and indexed
at runtime through a `ParameterContainer` (`pc.getParameter(i)`, `pc.getObservable(i)`). Each device
function ends with `pc.incrementIndex(...)` to advance the cursor to the next PDF in the tree — getting
these counts wrong silently corrupts evaluation of every following PDF. See `src/PDFs/basic/GaussianPdf.cu`
for the minimal template, and `docs/CreatingPDFs.md`.

`fptype` is the float/double scalar; use `RO_CACHE(...)` when reading event data on the device.

### Evaluation pipeline

`FitManager` (alias for `FitManagerMinuit2`; `FitManagerMinuit1` exists only with ROOT) drives a Minuit
minimizer. Minuit varies the `Variable`s; for each step `GooPdf::calculateNLL()` runs a Thrust
transform-reduce over the dataset using a `MetricTaker`/`MetricPointer` (the per-event metric, e.g. NLL),
after `normalize()` integrates the PDF over a grid. Data lives in `UnbinnedDataSet` / `BinnedDataSet`,
backed by `thrust::device_vector`.

### Source layout

- `include/goofit/`, `src/goofit/` — core: `PdfBase`, `Variable`/`Observable`, datasets, `Application`
  (CLI11-based, `GOOFIT_PARSE(app)`), `FitManager`, `FCN`, `Faddeeva`, `GlobalCudaDefines.h` (backend
  abstraction macros: `MEMCPY*`, `__constant__`, etc.).
- `src/PDFs/{basic,combine,physics,utilities}` — PDF implementations. `physics/` is the largest:
  Dalitz-plot / amplitude analysis (`Amp3Body*`, `Amp4Body*`, resonances, lineshapes, spin factors,
  K-matrix, time-dependent resolution models). `combine/` has `AddPdf`, `ProdPdf`, etc.
- `python/` — pybind11 bindings mirroring the C++ tree (`python/PDFs/...`); `examples/`, `tests/`;
  built via scikit-build (`setup.py`, `pyproject.toml`).
- `extern/` — git submodules: thrust, Minuit2, Eigen, fmt, CLI11, catch2, pybind11, MCBooster, etc.
- `examples/` — runnable analyses; each subdir has a `CMakeLists.txt` using `goofit_add_executable()` /
  `goofit_add_directory()` helpers defined in the top-level `CMakeLists.txt`.

### Adding a PDF or example

New example: make a dir with a `CMakeLists.txt` containing `goofit_add_directory()` and
`goofit_add_executable(Name Name.cu)`, then add the dir name to `examples/CMakeLists.txt`. New PDF:
add the `.cu` to the matching `src/PDFs/*/CMakeLists.txt` and the header under `include/goofit/PDFs/*`.

## Conventions

- C++11, `GooFit` namespace throughout. Includes use the prefixed form: `#include <goofit/Variable.h>`.
- Host functions are `__host__`, device evaluators `__device__`; trailing-return-type style
  (`-> fptype`) is used across the codebase.
- `Variable` (fit parameter) vs `Observable` (data dimension) are distinct types — don't access members
  directly, use getters/setters or treat the object as its value.
