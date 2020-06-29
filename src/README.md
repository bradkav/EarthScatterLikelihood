## Source code

Everything is written in fortran. It should conform to the fortran95 standard and it should compile straightforwardly using the Makefile. The default target is to make `calcContour` which calculates p-value tables for projected reconstructions. You can also use the Makefile to compile some other things, such as event generators and generic likelihood calculators (for example, if you're interested in interfacing with something like MultiNest). Just run `make all` to get everything, or look at the specific targets.

#### Contour calculator

To run parameter reconstructions, you can then do:
```
./calcContour M_X SIGMA_B DATA FIX_MASS LAT_DET OUTPATH
```

The command line arguments are as follows:
* `M_X` - WIMP mass in GeV  
* `SIGMA_B` - Benchmark WIMP-proton cross section in cm^2  
* `DATA` - flag for which data to use (1 = Energy + time, 2 = time only, 3 = energy only)  
* `FIX_MASS` - flag for whether the WIMP mass should be kept fixed (1 = fix to benchmark value, 0 = profile in range [0.058, 0.5] GeV)  
* `LAT_DET` - detector latitude (in degrees, over the range [-90, 90])  
* `OUTPATH` - output folder to save results to (this will be `./results/OUTPATH/`)


#### Event generators

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


#### Generic likelihood calculators

The subroutines for calculating likelihoods are in `like.f90`:

* `loglike` - calculate the unbinned likelihood including event energies and timing
* `loglike_Eonly` - calculate the unbinned likelihood including only event energies
* `loglike_counts` - calculate Poisson likelihood for the observed number of counts

Each subroutine accepts three arguments:

* `Cube` - which is the array of parameter. For now, this should be the cross section and local DM density.
* `slhood` - the output variable, to which the log-likelihood is assigned.
* `binned` - [optional, default=.False.] if `.True.`, use binned Asimov data for calculating the likelihood.

Hopefully these should be easy enough to interface with MultiNest (there's an example wrapper function for `slikelihood` in `like.f90`).


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
* (Incomplete...)

#### The inner workings

The basic idea is that calling `initialise_modulation` will read in all the velocity and density tables (from the `DaMaSCUS_results` folder), ready to interpolate. Then when you call `dRdE(E, t, ...)` (the event rate), it calculates the isodetection angle of the detector at time `t`, calculates the value of `vmin` corresponding to `E` and then performs a linear interpolation to get the correct value of eta, the velocity integral. It's actually a 4-dimensional interpolation, because it also interpolates for the value of the mass and cross section (in this case though, it's a linear interpolation in log10(m) and log10(sigma_SI), which turns out to give much smoother results). 

Note that you should always call the `initialise_modulation` subroutine (from `modulation.f90`) before doing any calculations of DD rates, in order to load in the tabulated velocity distributions. Don't worry about calling it multiple times, it will only load in the tables the first time it's called. 

There's also a similar interpolation which happens for the value of the local density at a given isodetection angle (and cross section).

Functions like `Nevents` then just integrate over time and energy to obtain the number of events per kg of detector. We then just multiply by the detector mass and any correction if we want to change the 'free' value of the local density and this gives the total event rate. There are a number of other functions which calculate, for example, the number of events at a fixed time (`Nevents_fixedt`) or the recoil spectrum integrated over time (`dRdE_res_tint`).

Experimental parameters are currently hard-coded in the module `expt.f90`. If you want to change the resolution, mass, exposure, latitude, background rate etc., just edit the appropriate parameter in `expt.f90`. In future, we can do something more sophisticated and read the parameters from file.

#### To-Do

* Extend to lower values of m_x by using scaling relations
