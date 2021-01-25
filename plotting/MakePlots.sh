#python3 PlotContours_zoom.py -m_x 0.1 -hemisphere N -runID $1 
#python3 PlotContours_zoom.py -m_x 0.1 -hemisphere S -runID $1

python3 PlotContours_all.py -hemisphere N -m_x 0.100 -runID $1 -plottext '$E_\mathrm{th} = 60\,\mathrm{eV}$; $\sigma_E = 18\,\mathrm{eV}$'
python3 PlotContours_all.py -hemisphere S -m_x 0.100 -runID $1 -plottext '$E_\mathrm{th} = 60\,\mathrm{eV}$; $\sigma_E = 18\,\mathrm{eV}$'

#python PlotLikelihoodExamples_mx100.py -m_x 0.1 -hemisphere S -sigtext "sig-34.00"
