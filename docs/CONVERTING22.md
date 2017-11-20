# Converting from older code for GooFit 2.2

## Writing PDFs

The indexing used internally has been rewritten for performance and better read-ability.  Writing a new PDF will be broken into multiple sections.  First, we have the distinction between `device` function and `class`.  There are numerous examples to utilize and meet the desired functionality.

We will start with creating a basic class.  Each class needs to inherit from GooPdf.  This exposes the appropriate methods to use in a constructor for initializing your new PDF.  Three types of information are utilized, which are constants, observables, and parameters.  Any value that will be change by MINUIT will be a `Parameter`.  A parameter can be added using the `registerParameter` function, which takes a single Parameter object.  `Observables` are passed through the constructor, and should be set at construction.  The 'index' into the observables list is stored, not the observable value itself.  Constants can be added using `registerConstant`.  Note, that there is not a distinction between integer and floating point, so currently all values are stored as an integer.  It may be helpful for the developer to track the indexes, which will be used later.

Here is an example for creating the constructor of the class:

```cpp
host__ GaussianPdf::GaussianPdf(std::string n, Observable _x, Variable mean, Variable sigma)
    : GooPdf(n, _x) {
    registerParameter(mean);
    registerParameter(sigma);

    initialize();
}

```

In the above code, GooPdf will appropriately `registerObservable(_x)`, so the developer does not need to do this registration.

Once everything has been registered, each constructor calls `initialize`.  In addition, each PDF needs to overload the function `recursiveSetIndices`, which provides the conversion of all information into a format used in device code.  Unless the PDF needs extra information, this should suffice:

```cpp
__host__ void GaussianPdf::recursiveSetIndices() 
{
    GET_FUNCTION_ADDR(ptr_to_Gaussian);

    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_Gaussian");
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx                               = num_device_functions++;

    populateArrays();
}

```

## Device functions

These are function implementations used in the minimization algorithm.  There are two passes currently performed.  First, is the normalisation phase.  This routine is run after updating all parameters from MINUIT.  The next function run is the summation of the negative log likeli-hood.  It is important to note the distinction, an example of this distinction is with the dalitz plot.  Caches are pre-calculated, and stored for use within the normalisation.

When writing the device function, all parameters, observables, etc are accessible through the ParameterContainer.  It is important as the developer to also increment the ParameterContainer using two helper functions, incrementIndex(func, parms, obs, cons, norms), and with incrementIndex().  If possible, please always use the first function, as the latter function requires four very slow memory reads.

Here is an example using a gaussian on how this is setup.

```
__device__ fptype device_Gaussian(fptype *evt, ParameterContainer &pc) {
    int id       = pc.getObservable(0);
    fptype mean  = pc.getParameter(0);
    fptype sigma = pc.getParameter(1);
    fptype x     = evt[id];

    pc.incrementIndex(1, 2, 0, 1, 1);
    fptype ret = exp(-0.5 * (x - mean) * (x - mean) / (sigma * sigma));

    return ret;
}
```

The observable index is fetched, then the two parameters are accessed.  


## Advanced PDFs

The above example PDF was quite simple, and didn't have anything complicated going on.  This section will discuss some of the problems and solutions encountered with integrating the DalitzPlotPdf.  The first problem encountered is with calculating and saving a cache.  Two operators work on the data, and store these values separately.  This is the Calculator and Integrator functions.  Since a resonance now inherits from a GooPdf, a resonance also contains two parameters, these 16 resonances need to exist within the PDF creation system.  When the Calculator and integrator functions are called, first we need to seek to the appropriate dalitz PDF, then we need to seek to the appropriate resonance(s) functions to be called.  Under `device_dalitz`, the same thing needs to happen.  The device function is executed, then the resonance PDFs need to be skipped in order to get to the efficiency function.  The c_motherMass, c_daug1Mass, etc. were moved to constant memory as an optimization.  Duplicating these at least 17x will be far too inefficient.


