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
parser.add_argument('-m_x', '--m_x', help='DM mass in GeV', default = 0.2)
parser.add_argument('-hemisphere','--hemisphere', help='Hemisphere of the experiment (N or S)', type=str, default="N")
args = parser.parse_args()

runID = args.runID
hemisphere = args.hemisphere
m_x = args.m_x  #DM mass in GeV
m_str = str(int(m_x*1000)) #String of DM mass in MeV

if (hemisphere == "N"):
    lat_text = "Northern Hemisphere ($46^\circ$N)"
elif (hemisphere == "S"):
    lat_text = "Southern Hemisphere ($37^\circ$S)"


fig =  plt.figure(figsize=(7,5))

ax = fig.add_subplot(111) 
axes = []
axes.append(fig.add_subplot(141))
axes.append(fig.add_subplot(142))
axes.append(fig.add_subplot(143))
axes.append(fig.add_subplot(144))

plt.subplots_adjust(wspace=0.1)

ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

file_labels = ["A", "B", "C", "D"]

sig_labels = [r"5 \times 10^{-35}", r"8.5 \times 10^{-35}", r"1 \times 10^{-34}", r"2 \times 10^{-34}"]

ratio_labels = [r"0.5 ", r"0.85 ", r"", r"2 "]

m_list = np.linspace(0.1, 0.5, 100)
rho_list = np.linspace(0.01, 1.0, 90)

for i in range(4):
    
    m, rho, L = np.loadtxt("../results/example_" + hemisphere + "_" + file_labels[i] + ".txt", unpack=True, usecols=(0, 2,3))
    L_grid = L.reshape(90, 100)
    
    cont = axes[i].contourf(m_list, rho_list, np.clip(L_grid, -20, 0.1), levels=np.linspace(-20, 0, 11))
    irho,im = np.unravel_index(np.argmax(L_grid), (90, 100))
    
    axes[i].plot([0.2, 0.2], [0.01, 0.8], linestyle='--',color='w')
    
    axes[i].plot(m_list[im], rho_list[irho],"k^", mew=2, ms=4)

    if (i > 0):
        axes[i].get_yaxis().set_ticklabels([])
    axes[i].text(0.3, 0.92, file_labels[i], color='w', ha='center', va='center', fontsize=18)
    axes[i].set_xticks([0.1, 0.2, 0.3, 0.4, 0.5])
    axes[i].set_xticklabels([" ", "0.2", " ", "0.4", " "])
    axes[i].tick_params(axis='x', colors='white')
    axes[i].tick_params(axis='y', colors='white')
    axes[i].tick_params(labelcolor='k')
    axes[i].spines['top'].set_color('white')
    axes[i].spines['bottom'].set_color('white')
    axes[i].spines['left'].set_color('white')
    axes[i].spines['right'].set_color('white')

    
    #axes[i].text(0.3, 0.07, r'$' + sig_labels[i] + '\, \mathrm{cm}^2$', color='w', ha='center', va='center', fontsize=12)
    #0.3, 0.07
    axes[i].text(0.3, 0.85, r"$\sigma_{p}^{\mathrm{SI}} = " + ratio_labels[i] + " \,\sigma_{p}^{\mathrm{SI}}{}'$",color='w', ha='center', va='center', fontsize=13)
    
axes[1].text(0.2, 1.02,lat_text,fontsize=14)
    
cb_ax = fig.add_axes([0.94, 0.09, 0.02, 0.8])
cbar = fig.colorbar(cont, cax=cb_ax, label=r'$\Delta \log \mathcal{L}$')


ax.set_xlabel(r'$m_\chi$ [GeV]')
axes[0].set_ylabel(r'$\rho_\chi$ [GeV/cm$^3$]')

    
plt.savefig("../plots/Likelihood_examples_" + hemisphere + ".pdf", bbox_inches='tight')
plt.show()

