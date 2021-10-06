



#TEXT="\$E_\mathrm{th} = 60\,\mathrm{eV}\$; \$\sigma_E = 18\,\mathrm{eV}\$ (Germanium, 1 kg day)"
#ID=Final_60eV_res18eV_EDE_clip2_refine
#MASS=0.400

TEXT="\$E_\mathrm{th} = 10\,\mathrm{eV}\$; \$\sigma_E = 3\,\mathrm{eV}\$ (Sapphire, 1 kg day)"    
ID="Final_10eV_res3eV_Sapphire_clip2_refine"
MASS=0.100

#python3 PlotContours_all.py -hemisphere N -m_x $MASS -runID $ID -plottext "$TEXT"
#python3 PlotContours_all.py -hemisphere S -m_x $MASS -runID $ID -plottext "$TEXT"

python3 PlotContours_zoom.py -m_x $MASS -hemisphere N -runID $ID 
python3 PlotContours_zoom.py -m_x $MASS -hemisphere S -runID $ID

#python3 PlotLikelihoodExamples_mx100.py -m_x 0.1 -hemisphere S -sigtext "sig-32.50"
