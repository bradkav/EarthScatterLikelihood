## Likelihood calculator for Earth-scattering DM

#### Summary

Currently, the likelihood calculator and event generator use **a fixed DM mass (0.5 GeV)**, because we only have DaMaSCUS simulations for this mass. The free parameters are therefore the DM cross section and local DM density. Because we only have a finite number of simulations, the code only works for cross sections in the range [1e-45, 3e-34] cm^2.

Everything is written in fortran. It should conform to the fortran95 standard and it should compile straightforwardly using the Makefile: this will make all the modules, as well as `test` (used for testing the likelihood calculator) and `GenerateEvents` which is the command-line event generator (see below). Let me know if there are any problems compiling or interfacing with MultiNest.

Experimental parameters are currently hard-coded in the module `expt.f90`. If you want to change the resolution, mass, exposure, latitude, etc., just edit the appropriate parameter in `expt.f90`. In future, we can do something more sophisticated and read the parameters from file.

The event generator should be fast - it shouldn't take more than about 10-20 seconds no matter how many events you need to generate. The likelihood calculators should also be reasonably fast. For the full likelihood, you should be able to do about 100 evaluations/second with O(1000) events. If you want to use the energy-only likelihood, it's a little slower (because there's an extra time integral), but this can probably be sped up using some kind of pre-tabulation.

#### Use


You can generate events by running 

```
./GenerateEvents SIGMA_SI RHO
```
where

* SIGMA_SI is the WIMP-nucleon xsec in cm^2
* RHO is the local density in GeV/cm^3

This will generate a sample of events and save it in the file `events.txt` (the first column contains the event times in Julian days, the second column event energies in keV).

Once you have an `events.txt` file, you can run `./test`, which will test the likelihood calculator after loading in `events.txt`. It will spit out the file `likes.txt` which contains a list of cross section values (column 1), followed by log-likelihood values (columns 2-4, including energy+time, energy-only and counts-only).

The subroutines for calculating likelihoods are in `like.f90`:

* `slikelihood` - calculate the unbinned likelihood including event energies and timing
* `slikelihood_Eonly` - calculate the unbinned likelihood including only event energies
* `slikelihood_counts` - calculate Poisson likelihood for the observed number of counts

Each subroutine accepts two arguments:

* `Cube` - which is the array of parameter. For now, this should be the cross section and local DM density.
* `slhood` - the output variable, to which the log-likelihood is assigned.

Hopefully these should be easy enough to interface with MultiNest.

Note that you should always call the `initialise_modulation` subroutine (from `modulation.f90`) before doing any calculations of DD rates, in order to load in the tabulated velocity distributions. Don't worry about calling it multiple times, it will only load in the tables the first time it's called. 

#### List of files:

* `test.f90` - some testing code that I'm using to make sure that all the calculations are correct.
* `like.f90` - the likelihood calculators.
* `utils.f90` - functions for doing interpolation, integration, initialising arrays etc.
* `LabFuncs.f90` - functions for calculation the lab velocity as a function of time (useful for calculating the isodetection angle).
* `modulation.f90` - functions for loading in, storing and interpolating the results of the DaMaSCUS simulations for the local DM density and velocity distribution.
* `DDrate.f90` - functions for calculating direct detection scattering rates, velocity integrals, event numbers etc.
* `expt.f90` - contains (hard-coded) experimental parameters (such as exposure time, detector mass, resolution, latitude, etc.)
* `GenerateEvents.f90` - event generator, samples recoil events for a given cross section and local density and saves the results to file.
* `rnglib.f90` and `ranlib_poisson.f90` - random number generators needed to generate Poisson random variables. Note that these files have been taken and adapted from https://people.sc.fsu.edu/~jburkardt/f_src/rnglib/rnglib.html. It's possible that some of these could be removed and consolidated if you can come up with a simple Poisson generator. 

#### The inner workings

The basic idea is that calling `initialise_modulation` will read in all the velocity and density tables, ready to interpolate. Then when you call `dRdE(E, t, ...)` (the event rate), it calculates the isodetection angle of the detector at time `t`, calculates the value of `vmin` corresponding to `E` and then performs a linear interpolation to get the correct value of eta, the velocity integral. It's actually a 3-dimensional interpolation, because it also interpolates for the value of the cross section (in this case though, it's a linear interpolation in log10(sigma_SI), which turns out to give much smoother results). 

There's also a similar interpolation which happens for the value of the local density at a given isodetection angle (and cross section).

Functions like `Nevents` then just integrate over time and energy to obtain the number of events per kg of detector. We then just multiply by the detector mass and any correction if we want to change the 'free' value of the local density and this gives the total event rate. There are a number of other functions which calculate, for example, the number of events at a fixed time (`Nevents_fixedt`) or the recoil spectrum integrated over time (`dRdE_res_tint`).

#### To-Do

* Add header information to events.txt files
* Add backgrounds and background likelihood
* Improve bounds-checking in `interp_eta_scalar` (in `modulation.f90`)