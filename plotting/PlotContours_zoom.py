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
parser.add_argument('-m_x', '--m_x', help='DM mass in GeV', type = float, default = 0.2)
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

#hemisphere = "S"
data = "1"

if (hemisphere == "N"):
    lat_text = "Northern Hemisphere ($46^\circ$N)"
elif (hemisphere == "S"):
    lat_text = "Southern Hemisphere ($37^\circ$S)"
    
if (data == "1"):
    data_str = "Energy+Time"
elif (data == "3"):
    data_str = "Energy-only"

#45.5 N
#37.1 S

plt.title(r'$m_\chi^\prime = ' + m_str + ' \,\mathrm{MeV}$; ' + lat_text)#; ' + data_str)
#print(np.geomspace(0.1, 0.5, 10))

#N_bench = 5
#sig_list = np.geomspace(3.16e-37, 1e-30, 14)
#sig_list = np.geomspace(1e-37, 1e-30, 8)

sig_list = [10**-34.0,]

for i, sig in enumerate(sig_list):
    IDstr = runID + data + "_" + hemisphere
    PT.plotContour_single(m_x, sig, IDstr, col="C0", overlay_mass=True, fill_contour=False)
    plt.plot(sig, 0.4, color='k', marker="+", mew=2, markersize=12)
    
plt.axvline(5e-35, linestyle=':', color='k')
plt.text(3.8e-35, 0.90, "A", color='k')

plt.axvline(8.5e-35, linestyle=':', color='k')
plt.text(6.6e-35, 0.90, "B", color='k')

plt.axvline(1e-34, linestyle=':', color='k')
plt.text(1.1e-34, 0.90, "C", color='k')

plt.axvline(2e-34, linestyle=':', color='k')
plt.text(2.18e-34, 0.90, "D", color='k')

plt.xlim(1e-35, 1e-33)
#plt.xlim(3e-36, 3e-34)
plt.ylim(0.01, 1e0)

plt.savefig("../plots/contour_" + runID + "_" + m_str + "_" + hemisphere + "_zoom.pdf", bbox_inches='tight')
plt.show()

