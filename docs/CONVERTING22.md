# Converting from older code for GooFit 2.2

## Writing PDFs

The indexing used internally has been rewritten for performance and better read-ability.  This indexing rewrite is being done in order to remove a double-memory access required for parameter lookups and constant lookups for every PDF.  Instead of accessing specifically three buffers, now four buffers are used containing PDF-specific information.  If the PDF is thought as a tree, then each PDF has values that are written sequentially into these buffers and are accessed into a 1D buffer of information.  If the buffers are accessed directly, the very first value is the total number of items for this PDF.  This is not recommended to access, please use the getter functions available in ParameterContainer.  Any optimizations that will improve memory access will be hidden within the ParameterContainer for improvements to all PDFs.

Writing new PDFs will be broken into multiple sections.  First, we have the distinction between `device` function and `class`.  There are numerous examples to utilize in order to meet the desired functionality.  Please post any questions or issues if encountered to the issue tracker.

We will start with creating a basic gaussian distribution class.  Each class needs to inherit from GooPdf.  This exposes the appropriate methods to use in a constructor for initializing the new PDF.  Three types of information are utilized, which are constants, observables, and parameters.  Any value that will be change by MINUIT will be a `Parameter`.  A parameter can be added using the `registerParameter` function, which takes a single Parameter object.  `Observables` are passed through the GooPdf constructor, and few instances in which the developer needs to configure these manually.  Please note that the 'index' into the observables list is stored, not the observable value itself.  Constants can be added using `registerConstant`.  There currently is not a distinction between integer and floating point, so currently all values are stored as a fptype.  It may be helpful for the developer to track the indexes, which will be used later.  Registering any Parameters and Constants are placed into a 0-based index, and access from the device function uses the same offset.

Here is an example for creating the constructor of the class:

```cpp
host__ GaussianPdf::GaussianPdf(std::string n, Observable _x, Variable mean, Variable sigma)
    : GooPdf("GaussianPdf", n, _x, mean, sigma) {

    registerFunction("device_Gaussian", device_Gaussian);

    initialize();
}
```

In the above code, GooPdf will appropriately `registerObservable(_x)` and `registerParameter`, so the developer does not need to do this registration. The first parameter is the Pdf class name.

Once everything has been registered, each constructor calls `initialize`.  This routine will setup each PDF with `setMetrics`, which is required to be done.


## Device functions

These are function implementations used in the minimization algorithm.  There are two passes currently performed.  First, is the normalization phase.  This routine is run after updating all parameters from MINUIT.  The next function run is the summation of the negative log likeli-hood.  It is important to note the distinction, an example of this distinction is with the dalitz plot.  Caches are pre-calculated, and stored for use within the normalization.

When writing the device function, all parameters, observables, etc are accessible through the ParameterContainer.  It is important as the developer to also increment the ParameterContainer using two helper functions, incrementIndex(func, parms, obs, cons, norms), and with incrementIndex().  If possible, please always use the first function, as the latter function requires four very slow memory reads.

Here is an example using a gaussian on how this is setup.

```cpp
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

The observable index is fetched, then the two parameters are accessed.  Since we know how much to increment our container, we can increment by 1 function, 2 parameters, 0, constants, 1 observable, and 1 normalization.


## Advanced PDFs

The above example PDF was quite simple, and didn't have anything complicated going on.  This section will discuss some of the problems and solutions encountered with integrating the DalitzPlotPdf.  This section will build off the previous section.

The first problem encountered is with calculating and saving the cache.  Two operators work on the data, and store these values separately.  This is the Calculator and Integrator functions.  Since a resonance needs to inherit from a GooPdf because a resonance also contains two parameters, these 16 resonances need to exist within the PDF creation system.  When the Calculator and Integrator functions are called, first we need to seek to the appropriate dalitz PDF, then we need to seek to the appropriate resonance(s) functions to be called.

Under `device_dalitz`, the function needs to increment, skip over the resonance functions to be able to call the efficiency function.  Also, conditional code should be removed from this particular PDF for branching and for skipping.  need to be skipped in order to get to the efficiency function.  The c_motherMass, c_daug1Mass, etc. were moved to constant memory as an optimization.  Duplicating these at least 17x will be far too inefficient.

One potential optimization is to pre-compute the 'resonance jump' such that we avoid using the `incrementIndex();` function.

Another issue arises if the developer needs to track the `end of the efficiency function` or track the `index to start the efficiency function`.  The developer will need to overload the `recursiveSetIndices` to save these ID's, and utilize them appropriately.  In the case of the DalitzPlot, we need to increment our resonance functions, save our function ID to begin efficiency, then call recursiveSetIndices on the efficiency function.
