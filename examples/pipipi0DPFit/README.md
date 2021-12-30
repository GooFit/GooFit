This example is to perform a Dalitz-plot for $D^0\to \pi^+\pi^-\pi^0$.


## Setup details:

All of the data files are downloaded automatically by the CMake build system.


## Running the fits

You can see the options with:

```
$ ./pipipi0DPFit -h
pipipi0 Dalitz fit example
Usage: ./examples/pipipi0DPFit/pipipi0DPFit [OPTIONS] SUBCOMMAND

GooFit:
  -q,--quiet                  Reduce the verbosity of the Application
  --nosplash Excludes: -q,--quiet
                              Do not print a splash
  --config TEXT=config.ini    An ini file with command line options in it

Options:
  -h,--help                   Print this help message and exit
  --minuit1                   Use Minuit 1 instead of Minuit 2

Subcommands:
  toy                         Toy MC Performance evaluation
  truth                       Truth Monte Carlo fit
  sigma                       Run sigma fit
  efficiency                  Run efficiency fit
  canonical                   Run the canonical fit
  background_dalitz           Run the background Dalitz fit
  background_sigma            Run background sigma fit
  background_histograms       Write background histograms
  run_gen_mc                  Run generated Monte Carlo fit
  make_time_plots             Make time plots
```

### Toy MC

All of the subcommands also have help. For example:

```
$ ./pipipi0DPFit toy -h
Toy MC Performance evaluation
Usage: ./examples/pipipi0DPFit/pipipi0DPFit toy [OPTIONS] [sample] [load]

Positionals:
  sample INT=0                Sample number to use
  load INT=1                  Number of times to load

Options:
  -h,--help                   Print this help message and exit
  -s,--sample INT=0           Sample number to use
  -l,--load INT=1             Number of times to load
  -p,--plot                   Also make plots
```

Toy MC study is added here for the purpose of performance evaluation. The command to execute toy fits is:

```
$./pipipi0DPFit toy 0 1
```

The first argument is to call for the toy fit function, and the second points to the zeroth toy sample file "dataFiles/toyPipipi0/dalitz_toyMC_000.txt". The third and last argument is the number of times this toy sample is to be loaded. That is to say, the larger this number is, the longer the fitting procedure takes to complete. Roughly, for N times the toy sample file is loaded, the mixing fit would last N x 20 seconds on one computer tested.

### Canonical

For the executable "pipipi0DPFit", it takes several command-line arguments, and will fail if the arguments are not properly given. Now only for the demonstration of Dalitz plot fits, you just run:

```
$ ./pipipi0DPFit canonical dataFiles/cocktail_pp_0.txt --blindSeed=0
```

The canonical subcommand calls for a "canonical" DP fit over the given "data" stored in "dataFiles/cocktail_pp_0.txt" as the third argument. The last argument is to disable the blinding of mixing parameters x and y which is only necessary for real data. The "data" in "dataFiles/cocktail_pp_0.txt" is completely from the simulation (MC cocktails of signals / different background sources). The included signals events are simulated with the input of x = y = +1%. The statistics of this sample is comparable to that of real data.
