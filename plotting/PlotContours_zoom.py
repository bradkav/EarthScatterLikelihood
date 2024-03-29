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

plt.title(r'$m_\chi^\prime = ' + m_str + ' \,\mathrm{MeV}$; ' + lat_text, fontsize=14)#; ' + data_str)
#print(np.geomspace(0.1, 0.5, 10))

#N_bench = 5
#sig_list = np.geomspace(3.16e-37, 1e-30, 14)
#sig_list = np.geomspace(1e-37, 1e-30, 8)

sig_list = [10**-32.5,]
#sig_list = [10**-33.5,]

ADD_SLICES = True

for i, sig in enumerate(sig_list):
    IDstr = runID + data + "_" + hemisphere
    PT.plotContour_single(m_x, sig, IDstr, col="k", overlay_mass=True, fill_contour=False)
    plt.plot(sig, 0.4, color='k', marker="+", mew=2, markersize=12)
    
if (ADD_SLICES):
#if ((sig_list[0]/1e34 - 1)**2 < 1e-3):
    #plt.axvline(5e-35, linestyle=':', color='k')
    #plt.text(3.8e-35, 0.90, "A", color='k')
    
    plt.axvline(1e-33, linestyle=':', color='k')
    plt.text(8e-34, 0.05, "A", color='k')
    
    plt.axvline(3.16e-33, linestyle=':', color='k')
    plt.text(2.5e-33, 0.05, "B", color='k')
    
    plt.axvline(1e-32, linestyle=':', color='k')
    plt.text(1.12e-32, 0.05, "C", color='k')

plt.xlim(sig_list[0]/10, sig_list[0]*10)
#plt.xlim(3e-36, 3e-34)
plt.ylim(0.01, 1e0)

plt.savefig("../plots/contour_" + runID + "_" + m_str + "_sig%.1f_"%(np.log10(sig_list[0]),) + hemisphere + "_zoom.pdf", bbox_inches='tight')
plt.show()

