## Plotting

This folder contains some plotting routines in python for loading in the p-value data and generating plots from the paper. Note that in order to reproduce the results of the paper "*Measuring the local Dark Matter density in the laboratory*" (arXiv:2004.XXXXX), you will also need the tabulated p-values, which were generated with the EarthScatterLikelihood code: [DOI:10.5281/zenodo.3739341](https://doi.org/10.5281/zenodo.3739341). Simply extract those data files into the [`results/`](results/) folder and you should be able to use the plotting routines.

To generate figures 2 and 3 from the paper, run:

```
python PlotContours_all.py -m_x 0.2 -hemisphere N -runID Final
python PlotContours_all.py -m_x 0.2 -hemisphere S -runID Final

python PlotContours_zoom.py -m_x 0.2 -hemisphere N -runID Final
python PlotContours_zoom.py -m_x 0.2 -hemisphere S -runID Final

python PlotLikelihoodExamples.py -hemisphere N
python PlotLikelihoodExamples.py -hemisphere S
```
