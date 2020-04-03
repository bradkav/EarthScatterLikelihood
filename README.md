## Likelihood calculator for Earth-scattering DM

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3725882.svg)](https://doi.org/10.5281/zenodo.3725882) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Last Update:** 02/04/2020

In order to reproduce the results of the paper "*Measuring the local Dark Matter density in the laboratory*" (arXiv:2004.XXXXX), you will also need the tabulated p-values, which were generated with this code: [DOI:10.5281/zenodo.3739341](https://doi.org/10.5281/zenodo.3739341). Simply extract those data files into the [results/](results/) folder.

#### Summary

Code for calculating event rates and likelihoods for Earth-scattering DM, released in association with arXiv:2004.XXXXX. 

The likelihood calculator and event generator currently work in the ranges:

* m_x = [0.1, 0.5] GeV  
* sigma_p^SI = [1e-40, 1e-30] cm^2

For now, we should limit our search to this range - I'm working on extending to lower masses using some kind of scaling relations. Higher cross sections can possibly be accomodated with some approximations.

Everything is written in fortran. It should conform to the fortran95 standard and it should compile straightforwardly using the Makefile: this will make all the modules, as well as `test` (used for testing the likelihood calculator),  `GenerateEvents` and `GenerateEvents_Asimov` which are the command-line event generators (see below). Let me know if there are any problems compiling or interfacing with MultiNest.

Experimental parameters are currently hard-coded in the module `expt.f90`. If you want to change the resolution, mass, exposure, latitude, background rate etc., just edit the appropriate parameter in `expt.f90`. In future, we can do something more sophisticated and read the parameters from file.


**NB: Requires results of the DaMaSCUS simulations.** The location of the DaMaSCUS simulation results can be specified as `data_dir` in the file `modulation.f90`.

#### Use


You can generate a random sample of events by running 

```
./GenerateEvents M_X SIGMA_SI RHO
```
where

* M_X is the WIMP mass in GeV
* SIGMA_SI is the WIMP-nucleon xsec in cm^2
* RHO is the local density in GeV/cm^3

This will generate a sample of events and save it in the file `events.txt` (the first column contains the event times in Julian days, the second column event energies in keV).

You can also generate a binned asimov dataset by running

```
./GenerateEvents_Asimov M_X SIGMA_SI RHO
```

This will output `events_Asimov.txt` (the first column contains the centres of the bins in time, the second the centres of the logarithmic bins in energy, and the third the expected number of events in each bin). You can change the number of bins in t and E by editing the `GenerateEvents_Asimov.f90` file directly, but for now it's 12 bins in each direction, which should be fine).

Once you have an `events_Asimov.txt`, you can run `./test`, which will test the likelihood calculator after loading in `events_Asimov.txt`. It will spit out the file `likes.txt` which contains a list of cross section values (column 1), followed by log-likelihood values (columns 2-4, including energy+time, energy-only and counts-only). There's also `./GridLike` which will calculate the likelihood on a 2-D grid of cross section and density (at fixed mass).

The subroutines for calculating likelihoods are in `like.f90`:

* `loglike` - calculate the unbinned likelihood including event energies and timing
* `loglike_Eonly` - calculate the unbinned likelihood including only event energies
* `loglike_counts` - calculate Poisson likelihood for the observed number of counts

Each subroutine accepts three arguments:

* `Cube` - which is the array of parameter. For now, this should be the cross section and local DM density.
* `slhood` - the output variable, to which the log-likelihood is assigned.
* `binned` - [optional, default=.False.] if `.True.`, use binned Asimov data for calculating the likelihood.

Hopefully these should be easy enough to interface with MultiNest (there's an example wrapper function for `slikelihood` in `like.f90`).

Note that you should always call the `initialise_modulation` subroutine (from `modulation.f90`) before doing any calculations of DD rates, in order to load in the tabulated velocity distributions. Don't worry about calling it multiple times, it will only load in the tables the first time it's called. 

#### Reconstructions

To run the parameter reconstructions, run 
```
./calcContour M_X SIGMA_B DATA FIX_MASS LAT_DET OUTPATH
```

The command line arguments are as follows:
* `M_X` - WIMP mass in GeV  
* `SIGMA_B` - Benchmark WIMP-proton cross section in cm^2  
* `DATA` - flag for which data to use (1 = Energy + time, 2 = time only, 3 = energy only)  
* `FIX_MASS` - flag for whether the WIMP mass should be kept fixed (1 = fix to benchmark value, 0 = profile in range [0.1, 0.5] GeV)  
* `LAT_DET` - detector latitude (in degrees, over the range [-90, 90])  
* `OUTPATH` - output folder to save results to (this will be `./output/OUTPATH/`)


#### List of files:

* `test.f90` - some testing code that I'm using to make sure that all the calculations are correct.
* `like.f90` - the likelihood calculators.
* `utils.f90` - functions for doing interpolation, integration, initialising arrays etc.
* `LabFuncs.f90` - functions for calculation the lab velocity as a function of time (useful for calculating the isodetection angle).
* `modulation.f90` - functions for loading in, storing and interpolating the results of the DaMaSCUS simulations for the local DM density and velocity distribution.
* `DDrate.f90` - functions for calculating direct detection scattering rates, velocity integrals, event numbers etc.
* `expt.f90` - contains (hard-coded) experimental parameters (such as exposure time, detector mass, resolution, latitude, etc.)
* `GenerateEvents.f90` - event generator, samples recoil events for a given mass, cross section and local density and saves the results to file.
* `GenerateEvents_Asimov.f90` - Asimov event generator.
* `rnglib.f90` and `ranlib_poisson.f90` - random number generators needed to generate Poisson random variables. Note that these files have been taken and adapted from https://people.sc.fsu.edu/~jburkardt/f_src/rnglib/rnglib.html. It's possible that some of these could be removed and consolidated if you can come up with a simple Poisson generator. 

#### The inner workings

The basic idea is that calling `initialise_modulation` will read in all the velocity and density tables, ready to interpolate. Then when you call `dRdE(E, t, ...)` (the event rate), it calculates the isodetection angle of the detector at time `t`, calculates the value of `vmin` corresponding to `E` and then performs a linear interpolation to get the correct value of eta, the velocity integral. It's actually a 4-dimensional interpolation, because it also interpolates for the value of the mass and cross section (in this case though, it's a linear interpolation in log10(m) and log10(sigma_SI), which turns out to give much smoother results). 

There's also a similar interpolation which happens for the value of the local density at a given isodetection angle (and cross section).

Functions like `Nevents` then just integrate over time and energy to obtain the number of events per kg of detector. We then just multiply by the detector mass and any correction if we want to change the 'free' value of the local density and this gives the total event rate. There are a number of other functions which calculate, for example, the number of events at a fixed time (`Nevents_fixedt`) or the recoil spectrum integrated over time (`dRdE_res_tint`).

The event rate calculations include a fixed flat background rate (`BG_rate` in `expt.f90`).

#### To-Do

* Extend to lower values of m_x by using scaling relations
