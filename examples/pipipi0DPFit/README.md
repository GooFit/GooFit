This example is to perform a Dalitz-plot for $D^0\to \pi^+\pi^-\pi^0$.
You need to make the current project by building GooFit with `-DGOOFIT_EXAMPLES=ON` (default).

Before running a Dalitz plot fit with command "pipipi0DPFit" for the first time, 
a bunch of ascii files containing data events (beefy ones) or pdf descriptions need to be downloaded from the following URL:

https://github.com/GooFit/GooFit/releases/download/v1.0.0/alltxtfiles4pipipi0DPFit.tgz
(Size: ~50 MB, ~144 MB after unpacking)

> Note:
> Previously from https://github.com/GooFit/GooFit/releases/download/v1.0.0/alltxtfiles4pipipi0DPFit.tgz

You then need to unpack the downloaded tar ball file in the current "pipipi0DPFit" directory via:

```
$ tar zxvf alltxtfiles4pipipi0DPFit.tgz
bkg_2_pdf_sideband_6slices.txt
bkg_3_pdf.txt
bkg_4_pdf.txt
signal_sigma_2slices_pdf.txt
signal_sigma_3slices_pdf.txt
signal_sigma_4slices_pdf.txt
signal_sigma_5slices_pdf.txt
signal_sigma_6slices_pdf.txt
signal_sigma_7slices_pdf.txt
signal_sigma_8slices_pdf.txt
dataFiles/
dataFiles/cocktail_pp_0.txt
dataFiles/bkgDalitz_4.txt
dataFiles/sideband1.txt
dataFiles/bkgDalitz_2.txt
dataFiles/efficiency_flat.txt
dataFiles/bkgDalitz_3.txt
dataFiles/sideband2.txt
```

For the executable "pipipi0DPFit", it takes several command-line arguments, and will fail if the arguments are not properly given. Now only for the demonstration of Dalitz plot fits, you just run:

```
$ ./pipipi0DPFit 4 dataFiles/cocktail_pp_0.txt blindSeed=0
```

The first argument is 4, which calls for a "canonical" DP fit over the given "data" stored in "dataFiles/cocktail_pp_0.txt" as the third argument. The last argument is to disable the blinding of mixing parameters x and y which is only necessary for real data. The "data" in "dataFiles/cocktail_pp_0.txt" is completely from the simulation (MC cocktails of signals / different background sources). The included signals events are simulated with the input of x = y = +1%. The statistics of this sample is comparable to that of real data.

Toy MC study is added here for the purpose of performance evaluation. The command to execute toy fits is:

```
$./pipipi0DPFit 0 0 1
```

The first argument is to call for the toy fit function, and the second points to the zeroth toy sample file "dataFiles/toyPipipi0/dalitz_toyMC_000.txt". The third and last argument is the number of times this toy sample is to be loaded. That is to say, the larger this number is, the longer the fitting procedure takes to complete. Roughly, for N times the toy sample file is loaded, the mixing fit would last N x 20 seconds @ Cerberus. 
