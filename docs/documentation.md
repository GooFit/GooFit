The GooFit Framework {#mainpage}
================================

Introduction
============

[GooFit](https://github.com/GooFit/GooFit) (\ref footnote1 "1")
is a framework for creating arbitrary probability density
functions (PDFs) and evaluating them over large datasets using nVidia
Graphics Processing Units (GPUs). New PDFs are written partly in
nVidia's CUDA programming language and partly in C++; however, no
expertise in CUDA is required to get started, because the
already-existing PDFs can be put together in plain C++.

Aside from the mass of unenlightened hominids who have not yet
discovered their need for a massively-parallel fitting framework, there
are three kinds of GooFit users:

-   Initiates, who write "user-level code" - that is, code which
    instantiates existing PDF classes in some combination. No knowledge
    of CUDA is required for this level. If your data can be described by
    a combination of not-too-esoteric functions, even if the combination
    is complicated, then user code is sufficient. Section
    [User code](@ref usercode) gives an example of how to write a simple fit.

-   Acolytes, or advanced users, who have grasped the art of creating
    new PDF classes. This involves some use of CUDA, but is mainly a
    question of understanding the variable-index organization that
    GooFit PDFs use. Section [New PDFs](@ref newpdfs) considers this organization
    in some depth.

-   Ascended Masters, or architects, who by extended meditation have
    acquired a full understanding of the core engine of GooFit, and can
    modify it to their desire (\ref footnote2 "2"). Section [Engine](@ref engine) gives a
    detailed narrative of the progress of a PDF evaluation through the
    engine core, thus elucidating its mysteries. It should only rarely
    be necessary to acquire this level of mastery; in principle only the
    developers of GooFit need to know its internal details.

Aside from considerations of the user's understanding, GooFit does
require a CUDA-capable graphics card to run on, with compute capability
at least 2.1. Further, you will need nVidia's CUDA SDK 7.0 or better, in particular
the `nvcc` compiler. Aside from this, GooFit is known to compile and run
on Fedora 14, Ubuntu 12.04, and OSX 10.8.4. It has been tested on the
Tesla, Fermi, and Kepler generations of nVidia GPUs.

Getting started
---------------

You will need to have a CUDA-capable device and to have the development
environment (also known as the software development kit or SDK) set up,
with access to the compiler `nvcc` and its libraries. If you have the
hardware, you can get the SDK from [nVidia's
website](https://developer.nvidia.com/gpu-computing-sdk).

With your CUDA environment set up, you can install GooFit thus:

-   Clone from the GitHub repository:

        git clone git://github.com/GooFit/GooFit.git
        cd GooFit

-   Compile with `cmake`:

        mkdir build
        cd build
        cmake ..
        make

    Do not be alarmed by warning
    messages saying that such-and-such a function's stack size could not
    be statically determined (currently hidden); this is an unavoidable (so far) side
    effect of the function-pointer implementation discussed in section
    [Engine](@ref engine).

-   Run the ‘simpleFitExample' program, which generates
    three distributions, fits them, and plots the results:

        cd examples/simpleFit
        ./simpleFit

    The expected output is a MINUIT log for three different fits, and
    three image files.

-   Run the Dalitz-plot tutorial, which fits a text file
    containing toy Monte Carlo data to a coherent sum of Breit-Wigner
    resonances:

        cd examples/dalitz
        ./dalitz dalitz_toyMC_000.txt

Quick troubleshooting: GooFit uses [FindCUDA](https://cmake.org/cmake/help/v3.7/module/FindCUDA.html), and expects
to find `root-config` in your path (ROOT optional, but used in most examples). Check the docs for FindCUDA if you need help locating your CUDA install.

The text file contains information about simulated decays of the \f$D^0\f$
particle to \f$\pi^+\pi^-\pi^0\f$; in particular, in each line, the second
and third numbers are the Dalitz-plot coordinates \f$m^2(pi^+\pi^0)\f$ and
\f$m^2(pi^-\pi^0)\f$. The `dalitz` program creates a PDF describing the
distribution of these two variables in terms of Breit-Wigner resonances,
reads the data, sends it to the GPU, and fits the PDF to the data - the
floating parameters are the complex coefficients of the resonances. The
expected output is a MINUIT fit log showing that the fit converged, with
such-and-such values for the real and imaginary parts of the resonance
coefficients.

User-level code {#usercode}
===============

From the outside, GooFit code should look like ordinary, object-oriented
C++ code: The CUDA parts are hidden away inside the engine, invisible to
the user. Thus, to construct a simple Gaussian fit, merely declare three
GooFit::Variable objects and a GooFit::GaussianPdf object that uses them, and create
an appropriate GooFit::UnbinnedDataSet to fit to:


Simple Gaussian fit {#listinggaussfit}
-------------------

```{.cpp}
int main (int argc, char** argv) {

  // Optional, but highly recommended. Based loosly
  // on TApplication.
  GooFit::Application app {"Simple Gaussian Fit", argc, argv};

  // Run the application parser, setup MPI if
  // needed, and exit if parsing failed
  GOOFIT_PARSE(app);


  // Create an object to represent the observable,
  // the number we have measured. Give it a name,
  // upper and lower bounds, and a number of bins
  // to use in numerical integration.
  Observable xvar {"xvar", -5, 5};
  xvar.setNumBins(1000);

  // A data set to store our observations in.
  UnbinnedDataSet data {xvar};

  // "Observe" ten thousand events and add them
  // to the data set, throwing out any events outside
  // the allowed range. In a real fit this step would
  // involve reading a previously created file of data
  // from an _actual_ experiment.
  TRandom donram(42);
  for (int i = 0; i < 10000; ++i) {
    fptype val = donram.Gaus(0.2, 1.1);
    if (fabs(val) > 5) {--i; continue;}
    data.addEvent(val);
  }

  // Variables to represent the mean and standard deviation
  // of the Gaussian PDF we're going to fit to the data.
  // They take a name, starting value, optional initial
  // step size and upper and lower bounds. Notice that
  // here only the mean is given a step size; the sigma
  // will use the default step of one-tenth of its range.
  Variable mean {"mean", 0, 1, -10, 10};
  Variable sigm {"sigm", 1, 0.5, 1.5};

  // The actual PDF. The Gaussian takes a name, an independent
  // (ie observed) variable, and a mean and width.
  GaussianPdf gauss {"gauss", xvar, mean, sigm};

  // Copy the data to the GPU.
  gauss.setData(&data);

  // A class that talks to MINUIT and GooFit. It needs
  // to know what PDF it should set up in MINUIT.
  FitManager fitter {&gauss};

  // The actual fit.
  fitter.fit();
  return 0;
}
```

Notice that, behind the scenes, GooFit assumes that there will be
exactly one top-level PDF and data set; it is not advised to break this
assumption unless you know what you are doing and exactly how you are
getting around it.

Data sets
---------

To create a data set with several dimensions, supply a `vector` of
`Observable`s:

```{.cpp}
Observable xvar {"xvar", -10, 10};
Observable yvar {"yvar", -10, 10};
vector<Observable> vars {xvar, yvar};
UnbinnedDataSet data(vars);
```

In this case, to fill the data set, set the GooFit::Variable values and call
the `addEvent` method without arguments:

```{.cpp}
xvar.setValue(3);
yvar.setValue(-2);
data.addEvent();
```

This will add an event with the current values of the GooFit::Variable list to
the data set. In general, where an unknown number of arguments are
wanted, GooFit prefers to use a `vector` of pointers.

Fit types
---------

By default, GooFit will do an unbinned maximum-likelihood fit, where the
goodness-of-fit metric that is minimized (\ref footnote3 "3") is the negative sum of
logarithms of probabilities, which is equivalent to maximizing the joint
overall probability:
\f{align}{
\cal P &=& -2\sum\limits_{events} \log(P_i)
\f}
where \f$P_i\f$
is the PDF value for event \f$i\f$.

To get a binned fit, you should create a GooFit::BinnedDataSet instead of the
GooFit::UnbinnedDataSet; the procedure is otherwise the same. Notice that the
GooFit::BinnedDataSet will use the number of bins that its constituent
GooFit::Variables have at the moment of its creation. Supplying a
GooFit::BinnedDataSet to a GooFit::GooPdf (which is the base class of all the
`FooPdf` classes such as `GaussianPdf`) will, by default, make it do a
binned negative-log-likelihood fit, in which the goodness-of-fit
criterion is the sum of the logarithms of the Poisson probability of
each bin:
\f{align}{
{\cal P} &=& -2*\sum\limits_{bins}(N * \log(E) - E)
\f}
where
\f$E\f$ is the expected number of events in a bin and \f$N\f$ is the observed
number.

There are two non-default variants of binned fits: A chisquare fit where
the error on a bin entry is taken as the square root of the number of
observed entries in it (or 1 if the bin is empty), and a "bin error" fit
where the error on each bin is supplied by the GooFit::BinnedDataSet. To do
such a fit, in addition to supplying the GooFit::BinnedDataSet (and providing
the errors through the `setBinError` method in the case of the bin error
fit), you should create a suitable `FitControl` object and send it to
the top-level GooFit::GooPdf:

```{.cpp}
Variable decayTime {"decayTime", 100, 0, 10};
BinnedDataSet* ratioData {&decayTime};
for (int i = 0; i < 100; ++i) {
  ratioData.SetBinContent(getRatio(i));
  ratioData.SetBinError(getError(i));
}

vector<Variable*> weights;
weights.push_back(new Variable("constaCoef", 0.03, 0.01, -1, 1));
weights.push_back(new Variable("linearCoef", 0, 0.01, -1, 1));
weights.push_back(new Variable("secondCoef", 0, 0.01, -1, 1));

PolynomialPdf poly {"poly", decayTime, weights};
poly.setFitControl(new BinnedErrorFit());
poly.setData(&ratioData);
```

The `FitControl` classes are `UnbinnedNLLFit` (the default),
`BinnedNLLFit` (the default for binned fits), `BinnedErrorFit` and
`BinnedChisqFit`.

Creating new PDF classes {#newpdfs}
========================

The simplest way to create a new PDF is to take the existing
GooFit::GaussianPdf class as a template. The existence of a `FooPdf.cu` file
in the `FPOINTER` directory is, because of Makefile magic, sufficient to
get the `Foo` PDF included in the GooFit library. However, a certain
amount of boilerplate is necessary to make the PDF actually work. First
of all, it needs a device-side function with a particular signature:


Signature of evaluation functions. {#listingfsign}
-------------------------------------------

```{.cpp}
__device__ fptype device_Gaussian (fptype* evt,
                                   fptype* p,
                                   unsigned int* indices);
```

Notice that this is a standalone function, not part of any class; `nvcc`
does not play entirely nicely with device-side polymorphism, which is
why we organize the code using a table of function pointers - a poor
man's implementation of the virtual-function lookups built into C++.
Second, we need a pointer to the evaluation function:

```{.cpp}
__device__ device_function_ptr ptr_to_Gaussian = device_Gaussian;
```

where `device_function_ptr` is defined (using `typedef`) as a pointer to
a function with the signature shown in the listing [here](@ref listingfsign):

```{.cpp}
typedef fptype (*device_function_ptr) (fptype*,
                                       fptype*,
                                       unsigned int*);
```

This pointer (\ref footnote4 "4") will be copied into the `device_function_table` array,
and its index in that array is the PDF's internal representation of "my
evaluation function".

Finally, the new PDF needs a bog-standard C++ class definition,
extending the GooFit::GooPdf superclass, which will allow it to be
instantiated and passed around in user-level code. [The indeces section](@ref subindexarray)
discusses what should happen in the constructor;
otherwise the class may have any supporting paraphernalia that are
necessary or useful to its evaluation - caches, lists of components,
pointers to device-side constants, whatever.

The indices array {#subindexarray}
-----------------

The heart of a PDF's organization is its index array, which appears in
the arguments to its device-side evaluation function as
`unsigned int* indices`. The index array stores the position of the
parameters of the PDF within the global parameter array; this allows
different PDFs to share the same parameters, as in two Gaussians with a
common mean. It also stores the position of the event variables,
sometimes called observables, within the event array passed to the
evaluation function; this is the argument `fptype* evt`.

The index array is created by the constructor of a PDF class; in
particular, the constructor should call `registerParameter` so as to
obtain the global indices of its parameters, store these numbers in a
`vector<unsigned int>` (conventionally called `pindices`), and pass this
`vector` to `initialize`. The PDF constructor should also call
`registerObservable` on each of the event variables it depends on.

The `initialize` method constructs the array that is used on the GPU
side, which consists of four parts. First is stored the number of
parameters, which is equal to the size of the `pindices vector`. Next
come the indices of the parameters, in the order they were put into
`pindices`. Then comes the number of observables, and finally the
indices of the observables, again in the order they were registered.

An example may be useful at this point. Consider the simple Gaussian PDF
constructor:

```{.cpp}
GaussianPdf::GaussianPdf (std::string name,
                          Observable _x,
                          Variable mean,
                          Variable sigma)
  : GooPdf(name, _x, mean, sigma) {
  registerFunction("ptr_to_Gaussian", ptr_to_Gaussian);
  initialize();
}
```

This is almost the simplest possible PDF: Two parameters, one
observable, no messing about! Notice that the call to
`registerFunction` sets the CUDA function this PDF is tied to. You can simply pass as many GooFit::Observable and GooFit::Variable parameters as you'd like to the parent constructor, or you can use `registerObservable` and `registerParameter` yourself.

The device-side index array for the Gaussian, assuming it is the only
PDF in the system, looks like this:

    index  0 1 2 3 4
    value  2 0 1 1 0

Here the initial 2 is the number of parameters - mean and sigma. Then
come their respective indices; since by assumption the Gaussian is the
only PDF we're constructing, these will simply be 0 and 1. Then comes
the number of observables, which is 1, and finally the index of the
observable - which, as it is the only observable registered, must be 0.
Now we can consider how the device-side code makes use of this:

```{.cpp}
__device__ fptype device_Gaussian (fptype* evt,
                                   fptype* p,
                                   unsigned int* indices) {
  fptype x = evt[indices[2 + indices[0]]];
  fptype mean = p[indices[1]];
  fptype sigma = p[indices[2]];

  fptype ret = exp(-0.5*(x-mean)*(x-mean)/(sigma*sigma));
  return ret;
}
```

The calculation of the Gaussian is straightforward enough, but let's
look at where the numbers `mean, sigma` and especially `x` come from.
The function is passed a pointer to the particular event it is to
calculate the value for, a global parameter array, and the index array.
The parameter array, in the case of a single Gaussian, consists simply
of the values for the mean and sigma in the current MINUIT iteration.
Let us replace the index lookups in those lines with their values from
above:

```{.cpp}
fptype mean = p[0];
fptype sigma = p[1];
```

which is exactly what we want. The fetching of `x` appears a little more
formidable, with its double `indices` lookup; it calls for some
explanation. First, `indices[0]` is the number of parameters of the
function; we want to skip ahead by this number to get to the ‘event'
part of the array. In the Gaussian, this is known at compile-time to be
2, but not every PDF writer is so fortunate; a polynomial PDF, for
example, could have an arbitrary number of parameters. (Or it might
specify a maximum number, say 10, and in most cases leave seven or eight
of them fixed at zero - but then there would be a lot of wasted
multiplications-by-zero and additions-of-zero.) Thus, as a convention,
lookups of event variables should always use `indices[0]` even if the
coder knows what that number is going to be. Then, 2 must be added to
this number to account for the space taken by the number-of-parameters
and number-of-observables entries in the array. So, replacing the first
level of lookup by the values, we have:

```{.cpp}
fptype x = evt[indices[4]];
```

and `indices[4]` is just 0; so in other words, `x` is the first
observable in the event. In the case of the single Gaussian, it is also
the *only* observable, so we've done quite a bit of work to arrive at a
zero that we knew from the start; but in more complex fits this would
not be true. The `x` variable could be observable number 5, for all we
know to the contrary in the general case. Likewise the mean and sigma
could be stored at positions 80 and 101 of the global parameter array.

Constants
---------

There are two ways of storing constants, or three if we count
registering a GooFit::Variable as a parameter and telling MINUIT to keep it
fixed. For integer constants, we may simply store them in the index
array; since it is up to the programmer to interpret the indices, there
is no rule that says it absolutely must be taken as an offset into the
global parameter array! An index can also store integers for any other
purpose - the maximum degree of a polynomial, flagging the use of an
optional parameter, or anything else you can think of. Indeed, this is
exactly what the framework does in enforcing the convention that the
first number in the index array is the number of parameters.

However, this will not serve for non-integer-valued constants. They must
either go through MINUIT as fixed parameters, or else go into the
`functorConstants` array. `functorConstants` works just like the global
parameters array, except that it does not update on every MINUIT
iteration since it is meant for storing constants. To use it, you should
first reserve some space in it using the `registerConstants` method,
which takes the number of constants you want as an argument and returns
the index of the first one. Usually you will want to put that index in
the `pindices` array. For example, suppose I want to store \f$\sqrt{2\pi}\f$
as a constant for use in the Gaussian. Then I would modify the
constructor thus:

```{.cpp}
__host__ GaussianPdf::GaussianPdf (std::string n,
                                   Variable* _x,
                                   Variable* mean,
                                   Variable* sigma)
  : GooPdf(_x, n)
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(mean));
  pindices.push_back(registerParameter(sigma));

  pindices.push_back(registerConstants(1));
  fptype sqrt2pi = sqrt(2*M_PI);
  MEMCPY_TO_SYMBOL(functorConstants, &sqrt2pi, sizeof(fptype),
                     cIndex*sizeof(fptype), cudaMemcpyHostToDevice);

  MEMCPY_FROM_SYMBOL((void**) &host_fcn_ptr, ptr_to_Gaussian, sizeof(void*));
  initialize(pindices);
}
```

Notice the member variable `cIndex`, which is set (and returned) by
`registerConstants`; it is the index of the first constant belonging to
this object. To extract my constant for use in the device function, I
look it up as though it were a parameter, but the target array is
`functorConstants` instead of the passed-in `p`:

```{.cpp}
__device__ fptype device_Gaussian (fptype* evt,
                                   fptype* p,
                                   unsigned int* indices) {
  fptype x = evt[indices[2 + indices[0]]];
  fptype mean = p[indices[1]];
  fptype sigma = p[indices[2]];
  fptype sqrt2pi = functorConstants[indices[3]];

  fptype ret = exp(-0.5*(x-mean)*(x-mean)/(sigma*sigma));
  ret /= sqrt2pi;
  return ret;
}
```

If I had registered two constants instead of one, the second one would
be looked up by `functorConstants[indices[3] + 1]`, not the
`functorConstants[indices[4]]` one might naively expect. This is because
the constant is stored next to the first one registered, but its *index*
is not stored at all; it has to be calculated from the index of the
first constant. Thus the `+1` must go outside the indices lookup, not
inside it! Keeping the levels of indirection straight when constructing
this sort of code calls for some care and attention.

Note that `functorConstants[0]` is reserved for the number of events in
the fit.

Program flow {#engine}
============

This section narrates the course of a fit after it is created, passing
through MINUIT and the core GooFit engine. In particular, we will
consider the example Gaussian fit shown in listing [Gauss fit](@ref listinggaussfit)
and look at what happens in these innocent-looking lines:

## Data transfer and fit invocation {#listingactualfit}

```{.cpp}
gauss.setData(&data);
FitManager fitter(&gauss);
fitter.fit();
```

Copying data
------------

The `setData` method copies the contents of the supplied GooFit::DataSet to
the GPU:

Internals of the setData method {#listingsetData}
--------------------------------------------------

```{.cpp}
setIndices();
int dimensions = observables.size();
numEntries = data->getNumEvents();
numEvents = numEntries;

fptype* host_array = new fptype[numEntries*dimensions];
for (int i = 0; i < numEntries; ++i) {
  for (obsIter v = obsBegin(); v != obsEnd(); ++v) {
    fptype currVal = data->getValue((*v), i);
    host_array[i*dimensions + (*v)->index] = currVal;
  }
}

gooMalloc((void**) &dev_event_array, dimensions*numEntries*sizeof(fptype));
cudaMemcpy(dev_event_array, host_array, dimensions*numEntries*sizeof(fptype), cudaMemcpyHostToDevice);
MEMCPY_TO_SYMBOL(functorConstants, &numEvents, sizeof(fptype), 0, cudaMemcpyHostToDevice);
delete[] host_array;
```

Notice the call to `setIndices`; this is where the indices of
observables passed to the PDF are decided and copied into the indices
array. This step cannot be done before all the subcomponents of the PDF
have had a chance to register their observables. Hence `setData` should
be called only after the creation of all PDF components, and only on the
top-level PDF.

The array thus created has the simple structure
`x1 y1 z1 x2 y2 z2 ... xN yN zN`, that is, the events are laid out
contiguously in memory, each event consisting simply of the observables,
in the same order every time. Notice that if the GooFit::DataSet contains
GooFit::Variables that have not been registered as observables, they are
ignored. If `setData` is called with an GooFit::BinnedDataSet object, the
procedure is similar except that each ‘event' consists of the
coordinates of the bin center, the number of events in the bin, and
either the bin error or the bin size. We will see later how the engine
uses the `dev_event_array` either as a list of events or a list of bins.

MINUIT setup
------------

Having copied the data to the GPU, the next task is to create the MINUIT
object that will do the actual fit; this is done by creating a
GooFit::FitManager object, with the top-level PDF as its argument, and calling
its `fit` method. The `fit` method does two things: First it calls the
`getParameters` method of the supplied PDF, which recursively gets the
registered parameters of all the component PDFs, and from the resulting
list of GooFit::Variables it creates MINUIT parameters by calling
`DefineParameter`. Second, it sets the method `FitFun` to be MINUIT's
function-to-minimize, and calls MINUIT's `mnmigr` method.

A few variants on the above procedure exist. Most obviously, ROOT
contains three implementations of the MINUIT algorithm, named `TMinuit`,
`TMinuit2`, and `TVirtualFitter` (\ref footnote5 "5"). GooFit provides the first two methods, and defaults to GooFit::FitManagerMinuit2. The interfaces differ, but the essential procedure is the
one described above: Define parameters, set function-to-minimize, run
MIGRAD. (NB: As of v0.2, GooFit has not recently been tested with
`MINUIT_VERSION` set to 2 or 3.) In the case of `TMinuit`, one can call
`setMaxCalls` to override the usual MINUIT limitation on the number of
iterations, although my experience is that this is not usually helpful
because running into the iteration limit tends to indicate a deeper
problem with the fit. Finally, the underlying `TMinuit` object is
available through the `getMinuitObject` method, allowing fine-grained
control of what MINUIT does, for example by calling `mnhess` in place of
`mnmigr`.

PDF evaluation
--------------

We have copied the data to the GPU, set up MINUIT, and invoked `mnmigr`.
Program flow now passes to MINUIT, which for purposes of this
documentation is a black box, for some time; it returns to GooFit by
calling the `FitFun` method with a list of parameters for which MINUIT
would like us to evaluate the NLL. `FitFun` translates MINUIT indices
into GooFit indices, and calls `copyParams`, which eponymously copies
the parameter array to `cudaArray` on the GPU. `FitFun` then returns the
value from `GooPdf::calculateNLL` to MINUIT, which absorbs the number
into its inner workings and eventually comes back with another set of
parameters to be evaluated. Control continues to pass back and forth in
this way until MINUIT converges or gives up, or until GooFit crashes.

The `calculateNLL` method does two things: First it calls the
`normalize` function of the PDF, which in turn will usually recursively
normalize the components; the results of the `normalize` call are copied
into the `normalizationFactors` array on the GPU. Next it calls
`sumOfNll` and returns the resulting value. Particular PDF
implementations may override `sumOfNll`; most notably `AddPdf` does so
in order to have the option of returning an ‘extended' likelihood, with
a term for the Poisson probability of the observed number of events in
addition to the event probabilities.

The `normalize` method, by default, simply evaluates the PDF at a grid
of points, returning the sum of all the values multiplied by the grid
fineness - a primitive algorithm for numerical integration, but one
which takes advantage of the GPU's massive parallelization. The fineness
of the grid usually depends on the `getNumBins` member of the observables;
in the case of the example Gaussian fit in listing [Gauss fit](@ref listinggaussfit),
the PDF will be evaluated at 1000 points, evenly spaced between -5 and 5.
However, this behavior can be overridden by calling the
`setIntegrationFineness` method of the PDF object, in which case the
number of bins (in each observable) will be equal to the supplied
fineness.

Stripped of complications, the essential part of the `normalize`
function is a call to `transform_reduce`:

normalization code. {#listingnormalization}
------------------------------------------

```{.cpp}
fptype dummy = 0;
static plus<fptype> cudaPlus;
constant_iterator<fptype*> arrayAddress(normRanges);
constant_iterator<int> eventSize(observables.size());
counting_iterator<int> binIndex(0);

fptype sum = transform_reduce(make_zip_iterator(
                               make_tuple(binIndex,
                                          eventSize,
                                          arrayAddress)),
                              make_zip_iterator(
                               make_tuple(binIndex + totalBins,
                                          eventSize,
                                          arrayAddress)),
              *logger, dummy, cudaPlus);
```

Here `normRanges` is an array of triplets `lower, upper, bins` for each
observable, created by the `generateNormRanges` method. The member
`logger` points to an instance of the `MetricTaker` class, which has an
operator method that Thrust will invoke on each bin index between the
initial value of zero and the final value of `totalBins-1`. This
operator method, which is invoked once per thread with a separate
(global) bin number for each invocation, calculates the bin center and
returns the value of the PDF at that point. The `dummy` and `cudaPlus`
variables merely indicate that Thrust should add (rather than, say,
multiply) all the returned values, and that it should start the sum at
zero. The `normalization` method returns this sum, but stores its
inverse in the `host_normalization` array that will eventually be copied
to `normalizationFactors` on the GPU; this is to allow the
micro-optimization of multiplying by the inverse rather than dividing in
every thread.

PDF implementations may override the `normalization` method, and among
the default PDFs, both `AddPdf` and `ProdPdf` do so to ensure that their
components are correctly normalized. Among the more specialized
implementations, `TddpPdf` overrides `normalize` so that it may cache
the slowly-changing Breit-Wigner calculations, and also because its time
dependence is analytically integrable and it is a good optimization to
do only the Dalitz-plot part numerically. This points to a more general
rule, that once a PDF depends on three or four observables, the
relatively primitive numerical integration outlined above may become
unmanageable because of the number of points it creates. Finally, note
that PDFs may, without overriding `normalize`, advertise an analytical
integral by overriding GooFit::GooPdf's `hasAnalyticIntegral` method to return
`true`, and then implementing an `integrate` method to be evaluated on
the CPU.

The `logger` object will appear again in the actual PDF evaluation,
performing a very similar function, so it is worth taking a moment to
consider in detail exactly what the `transform_reduce` call does. The
first two parameters (involving `make_tuple` calls) define the range of
evaluation: In this case, global bins (\ref footnote6 "6") 0 through \f$N-1\f$. They also
specify which `operator` method of `MetricTaker` should be called: It is
the one which takes as arguments two integers (the bin index and event
size) and an `fptype` array (holding the `normRanges` values), in that
order. Conceptually, Thrust will create one thread for each unique value
of the iterator range thus created - that is, one per global bin - and
have each thread invoke the indicated `operator` method. As a matter of
organization on the physical chip, it is likely that Thrust will
actually create a thousand or so threads and have each thread evaluate
as many bins as needed; but at any rate, the
`operator(int, int, fptype*)` method will be called once per global bin.
The last two arguments indicate that the return value should be
calculated as the sum of the return values from each `operator`
invocation, and that the sum should start at zero. Finally, the
`*logger` argument indicates the specific `MetricTaker` object to use,
which is important because this is where the function-pointer and
parameter indices are stored.

The `operator` does two things: First it calculates the bin centers, in
each observable, of the global bin:

Bin-center calculation {#listingbincenter}
----------------------

```{.cpp}
__shared__ fptype binCenters[1024*MAX_NUM_OBSERVABLES];

// To convert global bin number to (x,y,z...) coordinates:
// For each dimension, take the mod with the number of bins
// in that dimension. Then divide by the number of bins, in
// effect collapsing so the grid has one fewer dimension.
// Rinse and repeat.

int offset = threadIdx.x*MAX_NUM_OBSERVABLES;
unsigned int* indices = paramIndices + parameters;
for (int i = 0; i < evtSize; ++i) {
  fptype lowerBound = thrust::get<2>(t)[3*i+0];
  fptype upperBound = thrust::get<2>(t)[3*i+1];
  int numBins    = (int) floor(thrust::get<2>(t)[3*i+2] + 0.5);
  int localBin = binNumber % numBins;

  fptype x = upperBound - lowerBound;
  x /= numBins;
  x *= (localBin + 0.5);
  x += lowerBound;
  binCenters[indices[indices[0] + 2 + i]+offset] = x;
  binNumber /= numBins;
}
```

in the straightforward way, and stores the bin centers in a *fake
event*. Since events are just lists of observables, all that's necessary
is to keep track of which part of the `__shared__` (\ref footnote7 "7") `binCenters`
array is owned by this thread, look up the index-within-events of each
observable, and set the entries of the locally-owned part of
`binCenters` accordingly. This fake event is then sent to the PDF for
evaluation:

```{.cpp}
fptype ret = callFunction(binCenters+offset,
                          functionIdx,
                          parameters);
```

where `callFunction` is just a wrapper for looking up the function
referred to by `functionIdx` and calling it with the right part of the
parameter array:


Code to call device-side PDF implementations (some lines broken up for clarity) {#listingcallfunction}
----------------------------------


```{.cpp}
__device__ fptype callFunction (fptype* eventAddress,
                                unsigned int functionIdx,
                                unsigned int paramIdx) {
  void* rawPtr = device_function_table[functionIdx];
  device_function_ptr fcn;
  fcn = reinterpret_cast<device_function_ptr>(rawPtr);
  return (*fcn)(eventAddress,
                cudaArray,
                paramIndices + paramIdx);
}
```

This, finally, is where the `__device__` function from the PDF
definitions in section [New PDFs](@ref newpdfs) is called; we have now connected
all this engine code with the evaluation code for the Gaussian,
Breit-Wigner, polynomial, sum of functions, or whatever calculation we
happen to be doing today.

Having found the integral of the PDF, either using fake events as
outlined above or with an analytic calculation, we are now ready to find
the actual NLL, or sum of chi-squares, or other goodness-of-fit metric,
using the actual, observed events that we copied across in `setData`.
The procedure is similar to that for the normalization:

Goodness-of-fit evaluation {#listingnlleval}
-----------------


```{.cpp}
transform_reduce(make_zip_iterator(make_tuple(eventIndex,
                                              arrayAddress,
                                              eventSize)),
                 make_zip_iterator(make_tuple(eventIndex + numEntries,
                                              arrayAddress,
                                              eventSize)),
                 *logger, dummy, cudaPlus);   
```

Here the `*logger`, `dummy`, and `cudaPlus` arguments are doing the same
jobs as before. The tuple arguments, however, differ: In particular,
they are now indicating the range 0 to \f$N-1\f$ in *events*, not bins, and
`arrayAddress` this time points to the array of events, not to a set of
normalization triplets from which bin centers can be calculated. Since
the order of the arguments differs - it is now `int, fptype*, int` - a
different `operator` method is called:

Main evaluation operator (some lines broken up for clarity) {#listingmaineval}
----------------------------------------------------------

```{.cpp}
__device__ fptype MetricTaker::operator ()
  (thrust::tuple<int, fptype*, int> t) const {
  // Calculate event offset for this thread.
  int eventIndex = thrust::get<0>(t);
  int eventSize  = thrust::get<2>(t);
  fptype* eventAddress = thrust::get<1>(t);
  eventAddress += (eventIndex * abs(eventSize));

  // Causes stack size to be statically undeterminable.
  fptype ret = callFunction(eventAddress, functionIdx, parameters);

  // Notice assumption here! For unbinned fits the
  // eventAddress pointer won't be used in the metric,
  // so it doesn't matter what it is. For binned fits it
  // is assumed that the structure of the event is
  // (obs1 obs2... binentry binvolume), so that the array
  // passed to the metric consists of (binentry binvolume).
  void* fcnAddr = device_function_table[metricIndex];
  device_metric_ptr fcnPtr;
  fcnPtr = reinterpret_cast<device_metric_ptr>(fcnAddr);
  eventAddress += abs(eventSize)-2;
  ret = (*fcnPtr)(ret, eventAddress, parameters);
  return ret;
}
```

Observe that, for binned events, `eventSize` is negative; in this case
the event array looks like `x1 y1 n1 v1 x2 y2 n2 v2 ... xN yN nN vN`
where `x` and `y` are bin centers, `n` is the number of entries, and `v`
is the bin volume or error. This does not matter for the PDF evaluation
invoked by `callFunction`, which will just get a pointer to the start of
the event and read off the bin centers as event variables; hence the
`abs(eventSize)` in the calculation of the event address allows binned
and unbinned PDFs to be treated the same. However, it very much does
matter for the goodness-of-fit metric. Suppose the fit is the default
NLL: Then all the operator needs to do at this point is take the
logarithm of what the PDF returned, multiply by -2, and be on its way.
But if it is a chi-square fit, then it must calculate the expected
number of hits in the bin, which depends on the PDF value, the bin
volume, and the total number of events (\ref footnote8 "8"), subtract the observed
number, square, and divide by the observed number. Hence there is a
second function-pointer lookup, but now the `void*` stored in
`device_function_table` is to be interpreted as a different kind of
function - a "take the metric" function rather than a "calculate the
PDF" function. The `metricIndex` member of `MetricTaker` is set by the
`FitControl` object of the PDF; it points to one of the `calculateFoo`
functions:

Metric-taking functions {#listingmetrics}
-------

```{.cpp}
__device__ fptype calculateEval (fptype rawPdf,
                                 fptype* evtVal,
                                 unsigned int par) {
  // Just return the raw PDF value, for use
  // in (eg) normalization.
  return rawPdf;
}

__device__ fptype calculateNLL (fptype rawPdf,
                                 fptype* evtVal,
                                 unsigned int par) {
  rawPdf *= normalizationFactors[par];
  return rawPdf > 0 ? -log(rawPdf) : 0;
}

__device__ fptype calculateProb (fptype rawPdf,
                                 fptype* evtVal,
                                 unsigned int par) {
  // Return probability, ie normalized PDF value.
  return rawPdf * normalizationFactors[par];
}

__device__ fptype calculateBinAvg (fptype rawPdf,
                                 fptype* evtVal,
                                 unsigned int par) {
  rawPdf *= normalizationFactors[par];
  rawPdf *= evtVal[1]; // Bin volume
  // Log-likelihood of numEvents with expectation of exp
  // is (-exp + numEvents*ln(exp) - ln(numEvents!)).
  // The last is constant, so we drop it; and then multiply
  // by minus one to get the negative log-likelihood.
  if (rawPdf > 0) {
    fptype expEvents = functorConstants[0]*rawPdf;
    return (expEvents - evtVal[0]*log(expEvents));
  }
  return 0;
}

__device__ fptype calculateBinWithError (fptype rawPdf,
                                 fptype* evtVal,
                                 unsigned int par) {
  // In this case interpret the rawPdf as just a number,
  // not a number of events. Do not divide by integral over
  // phase space, do not multiply by bin volume, and do not
  // collect 200 dollars. evtVal should have the structure
  // (bin entry, bin error).
  rawPdf -= evtVal[0]; // Subtract observed value.
  rawPdf /= evtVal[1]; // Divide by error.
  rawPdf *= rawPdf;
  return rawPdf;
}

__device__ fptype calculateChisq (fptype rawPdf,
                                 fptype* evtVal,
                                 unsigned int par) {
  rawPdf *= normalizationFactors[par];
  rawPdf *= evtVal[1]; // Bin volume

  fptype ret = pow(rawPdf * functorConstants[0] - evtVal[0], 2);
  ret /= (evtVal[0] > 1 ? evtVal[0] : 1);
  return ret;
}
```

Notice the use of `normalizationFactors` in most of the metric
functions, and the special cases when the PDF or the observed number of
events is zero.

It is worth noting that the PDF evaluation function may itself call
other functions, either using `callFunction` or manually casting a
function index into other kinds of functions, as in the metric
calculation of listing [Main Eval](@ref listingmaineval). For example, in
`DalitzPlotPdf`, each resonance may be parametrized by a relativistic
Breit-Wigner, a Gaussian, a Flatte function, or more esoteric forms; so
the main function is supplied with a list of function indices and
parameter indices for them, and interprets the `void` pointer from
`device_function_table` as a specialized function type taking
Dalitz-plot location (rather than a generic event) as its argument. More
prosaically, `AddPdf` simply carries a list of PDF function indices and
indices of weights to assign them, and invokes `callFunction` several
times, multiplying the results by its weight parameters and returning
the sum.

We have now calculated the function value that we ask MINUIT to
minimize, for a single set of parameters; this value is passed back to
MINUIT, which does its thing and comes up with another set of parameters
for us, completing the loop. Ends here the saga of the fit iteration;
you now know the entire essential functionality of GooFit's core engine.

Existing PDF classes
====================

The GooFit PDFs, like ancient Gaul, are roughly divisible into three:

-   Basic functions, written because they are (expected to be)
    frequently used, such as the Gaussian and polynomial PDFs.

-   Combiners, functions that take other functions as arguments and spit
    out some combination of the inputs, for example sums and products.

-   Specialized PDFs, written for the \f$D^0\to\pi\pi\pi^0\f$ mixing
    analysis that is the driving test case for GooFit's capabilities.

In the lists below, note that all the constructors take pointers to
GooFit::Variable objects; rather than repetitively repeat "GooFit::Variable pointer"
in a redundantly recurring manner, we just say GooFit::Variable. Additionally,
the first argument in every constructor is the name of the object being
created; again this is not mentioned in every item. By convention,
constructors take observables first, then parameters.

Basic PDFs
----------

Basic PDFs are relatively straightforward: They take one or more
observables and one or more parameters, and implement operations that
are in some sense ‘atomic' - they do not combine several functions in
any way. Usually they have a reasonably well-known given name, for
example "the threshold function" or "a polynomial". The canonical
example is the Gaussian PDF.

- GooFit::ArgusPdf
- GooFit::BifurGaussPdf
- GooFit::BWPdf
- GooFit::CorrGaussianPdf
- GooFit::CrystalBallPdf
- GooFit::ExpGausPdf
- GooFit::ExpPdf

- GooFit::GaussianPdf
- GooFit::InterHistPdf
- GooFit::JohnsonSUPdf
- GooFit::KinLimitBWPdf
- GooFit::LandauPdf
- GooFit::NovosibirskPdf
- GooFit::PolynomialPdf
- GooFit::ScaledGaussianPdf
- GooFit::SmoothHistogramPdf
- GooFit::StepPdf
- GooFit::VoigtianPdf

-  GooFit::TrigThresholdPdf (was in physics)
-  GooFit::BinTransformPdf (was in combination)

Combination PDFs
----------------

These are the tools that allow GooFit to be more than a collection of
special cases. The most obvious example is a sum of PDFs - without a
class for this, you'd have to write a new PDF every time you added a
Gaussian to your fit.

- GooFit::AddPdf
- GooFit::CompositePdf
- GooFit::ConvolutionPdf
- GooFit::EventWeightedAddPdf
- GooFit::MappedPdf
- GooFit::ProdPdf

Specialized amplitude-analysis functions
----------------------------------------

These functions exist mainly for use in a specific physics analysis,
mixing in \f$D^0\to\pi\pi\pi^0\f$. Nonetheless, if you are doing a
Dalitz-plot analysis, you may find them, and conceivably even this
documentation, helpful.

-   GooFit::DalitzPlotPdf
-   GooFit::DalitzVetoPdf
-   GooFit::IncoherentSumPdf
-   GooFit::MixingTimeResolution (in `MixingTimeResolution_Aux.h`)
-   GooFit::ResonancePdf (subclasses not separately documented yet)
-   GooFit::TddpPdf
-   GooFit::ThreeGaussResolution
-   GooFit::TruthResolution

\anchor footnote1 1: Named in homage to RooFit, with the ‘G' standing for ‘GPU'.

\anchor footnote2 2: Although, if they are *Buddhist* masters, they don't even though
    they can, since they have transcended desire - and suffering with
    it.

\anchor footnote3 3: For historical reasons, MINUIT always minimizes rather than
    maximizing.

\anchor footnote4 4: You might ask, why not copy the function directly? The reason is
    that `cudaMemcpy` doesn't like to get the address of a function, but
    `nvcc` is perfectly happy to statically initialize a pointer. It's a
    workaround, in other words.

\anchor footnote5 5: These are, respectively, ancient FORTRAN code translated
    line-by-line into C++, almost literally by the addition of
    semicolons; someone's obsessively-detailed object-oriented
    implementation of the same algorithm, with the same spaghetti logic
    chopped into classes instead of lines of code; and what seems to be
    intended as a common interface for a large number of possible
    fitting backends, which falls a little flat since it has only the
    MINUIT backend to talk to. You pays your money and takes your
    choice.

\anchor footnote6 6: A global bin ranges from 0 to \f$n_1n_2\ldots n_N-1\f$ where \f$n_j\f$ is
    the number of bins in the \f$j\f$th variable and \f$N\f$ is the number of
    variables. In two dimensions, with three bins in each of \f$x\f$ and
    \f$y\f$, the global bin is given by \f$3b_y+b_x\f$, where \f$b_{x,y}\f$ is the
    bin number in \f$x\f$ or \f$y\f$ respectively, as shown here:

\f[
\begin{array}{l|ccc}
    2 & 6 & 7 & 8 \\
    1 & 3 & 4 & 5 \\
    0 & 0 & 1 & 2 \\
    \hline
      & 0 & 1 & 2
\end{array}
\f]
 where the leftmost column and bottom row indicate the
    \f$y\f$ and \f$x\f$ bin number.

\anchor footnote7 7: That is, `__shared__` for the default CUDA target.

\anchor footnote8 8: This is why `functorConstants[0]` is reserved for that value!
