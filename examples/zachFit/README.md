# ZachFit example

This example performs a staged fit measuring the mass difference
between the `D*(2010)+` and `D0` using `D*+ -> D0 pi+` events recorded by the
BaBar detector (approximately 477 inverse femtobarn). The D0 is
reconstructed in the `D0 -> K- pi+` and `D0 -> K- pi+ pi- pi+` decay
channels. The final stage of the fit convolves a resolution
function with a line shape to describe the data. The speedup
performance is somewhat of an disappointment for the moment
compared to the pipipi0DPFit example.

## Fit options

Fit options can be viewed by using the `-h` option: `./zachFit -h`. The `-p` option makes plots. The other options are:

| Fit Mode (`-m`) | Description       |
| --------------- | ----------------- |
| 0               | Unbinned          |
| 1               | Binned            |
| 2               | Binned Chisquared |


| Dataset (`-d`) | Description |
| --- | --- |
| 0 (simple) | Early testing sample for GooFit before nominal dataset was released. Resolution sample and MCdata for channel `D*+ -> D0 pi+`; `D0 -> K- pi+`. Samples are composed of events that pass the majority of selection criteria, but fail at least one of the stricter tracking cuts. The resulting resolution is worse than in the events of the nominal samples used in the official analysis/publication marked below as data set options "1" and "2".  |
| 1 (kpi)  | Nominal MC resolution sample and data for channel `D*+ -> D0 pi+`; `D0 -> K- pi+`         |
| 2 (k3pi) | Nominal MC resolution sample and data for channel `D*+ -> D0 pi+`; `D0 -> K- pi+ pi- pi+` |

## Data

The data are stored in the GitHub releases mechanism and are available [here](https://github.com/GooFit/GooFit/releases/download/v1.0.0/dataFiles_zachFit.tgz). The data were recorded by the BaBar detector and correspond to approximately 477 inverse femtobarns at and 40 below the Upsilon(4S) resonance.

The `.dat` files are delta M values for reconstructed `D*+ -> D0 pi+` events after selection criteria for the datasets that correspond to their filenames. The `D0` is reconstructed in the `D0 -> K- pi+` and `D0 -> K- pi+ pi- pi+` decay channels. Delta M is defined as the mass difference between the `D*+` and `D0` candidates.

Each dataset has two files: a resolution sample of MC simulated events (using a Geant4 simulation) and the real data events. The figures were created by the PlotOriginalData notebook file (also provided in pdf form) using the `.dat` files and are named according to the corresponding figures in the PRD publication.


## Fit stages

The description of the fit stages paraphrases the information
contained in the 2013 PRD publication of the analysis:

* [Phys. Rev. D 88, 052003 [ERRATUM: Phys. Rev. D 88, 079902]][main-paper]
* [arXiv link]

[main-paper]: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.88.052003
[arXiv link]: https://arxiv.org/abs/1304.5009


### Stage 1

The resolution is determined by fitting MC simulated events generated with an
infinitesimal `D*+` width.

The resolution function is composed of three Gaussian PDFs (independent means
and widths) and a ARGUS background PDF with an extra parameter allowed to float.
Simulation and validation fits showed we needed this additional non-Gaussian
component for an unbiased result. The ARGUS function is used to describe
events where the slow pion track includes a few hits from after it decayed
in flight to a muon. There are relatively few events where this decay in
flight mistake occurs, but validation samples indicate the description
does not have to be terribly precise in order to prevent bias.

### Stage 2

The fit to the data is composed of a signal PDF and a background PDF.
The background PDF is an ARGUS background parameterization.

The signal PDF is composed of the relativistic Breit-Wigner (RBW) line shape
convolved with the Gaussians of the resolution function. The resolution
function parameter values from the Stage 1 fit are held constant.
The RBW line shape (near a kinematic threshold) is described
by KinLimitBWPdf. The non-Gaussian, decay-in-flight, component
found in Stage 1 is omitted from the convolution since it corresponds
to a very small fraction of events. Two additional parameters are
included during the convolution to allow for possible differences
between resolutions in data and simulation. The widths of the Gaussians
found in Stage 1 are scaled by a common factor (1 + epsilon). The means of the Gaussians
found in Stage 1 are allowed a common offset delta.
