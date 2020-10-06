import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors

from scipy.integrate import cumtrapz, quad
from scipy.interpolate import interp1d
from scipy.stats import chi2

import PlottingTools as PT

import argparse
import os

#---------------
# MATPLOTLIB settings

mpl.rcParams.update({'font.size': 18,'font.family':'serif'})

mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['ytick.minor.width'] = 1
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rc('text', usetex=True)

mpl.rcParams['legend.edgecolor'] = 'inherit'

#---------------

parser = argparse.ArgumentParser(description='...')
parser.add_argument('-runID', '--runID', help='Text ID for the results to be plotted', type=str, default="Final")
parser.add_argument('-m_x', '--m_x', help='DM mass in GeV', type=float, default = 0.2)
parser.add_argument('-hemisphere','--hemisphere', help='Hemisphere of the experiment (N or S)', type=str, default="N")
args = parser.parse_args()

runID = args.runID
hemisphere = args.hemisphere
m_x = args.m_x  #DM mass in GeV
m_str = str(int(m_x*1000)) #String of DM mass in MeV



fig = plt.figure(figsize = (7,5))
ax = plt.gca()
ax.set_xscale("log")
#ax.set_yscale("log")

plt.xlabel(r"$\sigma_p^\mathrm{SI}$ [cm$^2$]")
plt.ylabel(r"$\rho_\chi$ [GeV/cm$^3$]")


if (hemisphere == "N"):
    lat_text = "Northern Hemisphere ($46^\circ$N)"
elif (hemisphere == "S"):
    lat_text = "Southern Hemisphere ($37^\circ$S)"


#45.5 N
#37.1 S

plt.title(r"$m_\chi' = " + m_str + " \,\mathrm{MeV}$ (fixed); "+ lat_text)# + data_str)

#List of cross section to plot for
sig_list = np.logspace(-35.5, -30.5, 11)

print("> Plotting results for runID:", runID, " (Hemisphere:",hemisphere, ")")

for i, sig in enumerate(sig_list):
    print("   Cross section:", sig)
    IDstr = runID + "3_" + hemisphere
    PT.plotContour_nomodulation(m_x, sig, IDstr, col="C4", ls='dashed', fixed_mass=True)
    
    IDstr = runID + "1_" + hemisphere
    PT.plotContour_modulation(m_x, sig, IDstr, col="C0", ls='solid', fixed_mass=True)

    plt.plot(sig, 0.4, color='k', marker="+", mew=2, markersize=8)
    
proxy_w = plt.Rectangle((-1,-1),1,1,fc = 'C0', alpha=0.8, edgecolor='C0', linewidth=1.5, linestyle='-')
proxy_wo = plt.Rectangle((-1,-1),1,1,fc = 'C4', alpha=0.8, edgecolor='C4', linewidth=1.5, linestyle='--')
              
plt.legend([proxy_wo, proxy_w], ["Energy-only", "Energy+timing"], loc='upper right',framealpha=0.9)
           
#plt.axvline(1e-34, linestyle='--', color='k')
#plt.axvline(1e-33, linestyle='--', color='k')
plt.xlim(1e-36, 2e-30)
#plt.xlim(3e-36, 3e-34)
plt.ylim(0, 1e0)
#plt.xticks(np.geomspace(1e-37, 1e-30, 8))

plt.savefig("../plots/contour_" + runID + "_" + m_str + "_" + hemisphere + "_all_fixedmass.pdf", bbox_inches='tight')
plt.show()

