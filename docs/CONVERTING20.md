# Converting from older code for GooFit 2.0


## Build system

The build system underwent a major upgrade in the move to CMake.  You should remove your old Makefiles and use the new `CMakeFiles.txt` files provided in examples - this should require
writing two lines of code instead of the 50 or so previously needed. 

An example of a `CMakeLists.txt` file:

```cmake
# Add links to all the data and source files in this directory
goofit_add_directory()

# Add an explicit link if you prefer finer control
goofit_add_link("myfile.txt")

# Add a GooFit program
goofit_add_executible(MyNewExample MyNewExample.cu)

# Add a GooFit PDF
include_directory(.)
goofit_add_pdf(MyPdf.cu MyPdf.h)
```

## Renames and moves

The folders that were introduced to keep the includes structured require modifications of source code, converting lines like `#include "Variable.hh"` to `#include "goofit/Variable.h"`. This modification can be done for you by running the provided script, `scripts/ModernizeGooFit.py` on your source files (requires Python and [Plumbum]). For example:

```bash
./scripts/ModernizeGooFit.py goofit_mypackage/*.cu
```

## Application

The new `GooFit::Application`, which is not required but provides GooFit options, like GPU selection and status, as well as MPI support and configurable command line options, is available by adding:

```cpp
#include "goofit/Application.h"

// Place this at the beginning of main
GooFit::Application app{"Optional discription", argc, argv};

// Command line options can be added here.

try {
    app.run();
} catch(const GooFit::ParseError &e) {
    return app.exit(e);
}
```

See [CLI11] for more details. The [pipipi0](../examples/pipipi0DPFit) example has an example of a complex set of options.

## Variable

The `Variable` class no longer allows direct access to members: please use getters and setters. For convenience, quick access to the value of a variable is also provided.

```cpp
Variable var {"var", 1.0, .01, -10, 10};

// Was fptype x = var.value;
fptype x = var.getValue();
fptype x = var; // Shortcut

// Was var.value = x;
var.SetValue(x);
var = x; // Shortcut

// Was istream >> var.value;
fptype x; istream >> x; val.setValue(x);
istream >> var; // Shortcut
```

The other values do not provide shortcut access. Copying a variable is explicitly disallowed; it was possible before but gave incorrect results.

## Namespaces

GooFit no longer leaks the `std::` namespace, and all of GooFit is now inside the GooFit namespace. You can add `using namespace GooFit;` at the top of your code.

[CLI11]:   https://github.com/CLIUtils/CLI11
[Plumbum]: https://github.com/tomerfiliba/plumbum
