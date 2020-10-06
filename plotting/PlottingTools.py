import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors

from scipy.integrate import cumtrapz, quad
from scipy.interpolate import interp1d
from scipy.stats import chi2

from scipy.ndimage import uniform_filter

import os

rootdir = "../results/"

def plotContour_single(m_x, sigma, outpath, col='C0', ls='solid', overlay_mass=False, add_fixed_mass=False, fill_contour=True):
    #print(np.geomspace(0.05, 0.5, 9))
    ID_str = "m%d_lsig%6.2f_mprof.txt"%(1000*m_x, np.log10(sigma))
    
    if not (os.path.isfile(rootdir + outpath + "/stat_sigma_" + ID_str)):
        ID_str = "m%d_lsig%6.2f.txt"%(1000*m_x, np.log10(sigma))
    
    sig_list = np.loadtxt(rootdir + outpath + "/stat_sigma_" + ID_str)
    rho_list = np.loadtxt(rootdir + outpath + "/stat_rho_" + ID_str)
    p_grid = np.loadtxt(rootdir + outpath + "/stat_p_" + ID_str)

    #print(np.min(p_grid))
    sig_grid, rho_grid = np.meshgrid(sig_list, rho_list)

    #print(0.4 - np.min(rho_grid[p_grid > 0.05]), np.max(rho_grid[p_grid > 0.05]) - 0.4)
    #print(sigma - np.min(sig_grid[p_grid > 0.05]), np.max(sig_grid[p_grid > 0.05]) - sigma)    
    #print(p_grid.shape)
    
    if (overlay_mass):
        m_grid = np.loadtxt(rootdir + outpath + "/stat_mbf_" + ID_str)
        plt.contourf(sig_grid, rho_grid, np.log10(m_grid), alpha=0.5, levels=np.log10(np.geomspace(0.5e-1, 5e-1, 17)),
                #ticks=np.arange(50, 549, 25),
                cmap=plt.get_cmap('seismic'), norm=colors.DivergingNorm(vmin=np.log10(0.05), vcenter=np.log10(0.1), vmax=np.log10(0.5)))
                 #norm=colors.LogNorm(vmin=0.1, vmax=0.5))
        cb = plt.colorbar(label=r'$\hat{m}_\chi$ [MeV]')
        #cb.set_ticks(np.log10(np.arange(50, 549, 25)))        
        #cb.set_ticklabels(["0.1","", "0.2","", "0.3","", "0.4","", "0.5"])
        cb.set_ticklabels(["50","", "90","", "160","", "280","", "500"])
        
        
    if (add_fixed_mass):
        ID_str = "m%d_lsig%6.2f_mfix.txt"%(1000*m_x, np.log10(sigma))
        sig_list2 = np.loadtxt(rootdir + outpath + "/stat_sigma_" + ID_str)
        rho_list2 = np.loadtxt(rootdir + outpath + "/stat_rho_" + ID_str)
        p_grid2 = np.loadtxt(rootdir + outpath + "/stat_p_" + ID_str)
        #m_grid = np.loadtxt("EarthScatterLikelihood/output/" + outpath + "/stat_mbf_" + ID_str)
        sig_grid2, rho_grid2 = np.meshgrid(sig_list2, rho_list2)
        
        plt.contour(sig_grid2, rho_grid2, p_grid2, levels = (0.05,), colors=col, linewidths=1.5, linestyles='--')
        
    #p_smooth = uniform_filter(p_grid, size=5)
        
    plt.contour(sig_grid, rho_grid, p_grid, levels = (0.05,), colors=col, linewidths=1.5, linestyles=ls)
    if (fill_contour):
        plt.contourf(sig_grid, rho_grid, p_grid, levels = (0.05,1), colors=col, alpha=0.5)


def plotContour_modulation(m_x, sigma, outpath, col='C0', ls='solid', fixed_mass=False):
    #print(np.geomspace(0.05, 0.5, 9))
    if (fixed_mass):
        ID_str = "m%d_lsig%6.2f_mfix.txt"%(1000*m_x, np.log10(sigma))
    else:
        ID_str = "m%d_lsig%6.2f_mprof.txt"%(1000*m_x, np.log10(sigma))
    
    if not (os.path.isfile(rootdir + outpath + "/stat_sigma_" + ID_str)):
        ID_str = "m%d_lsig%6.2f.txt"%(1000*m_x, np.log10(sigma))
    
    sig_list = np.loadtxt(rootdir + outpath + "/stat_sigma_" + ID_str)
    rho_list = np.loadtxt(rootdir + outpath + "/stat_rho_" + ID_str)
    p_grid = np.loadtxt(rootdir + outpath + "/stat_p_" + ID_str)

    #print(np.min(p_grid))
    sig_grid, rho_grid = np.meshgrid(sig_list, rho_list)

    #print(0.4 - np.min(rho_grid[p_grid > 0.05]), np.max(rho_grid[p_grid > 0.05]) - 0.4)
    #print(sigma - np.min(sig_grid[p_grid > 0.05]), np.max(sig_grid[p_grid > 0.05]) - sigma)    
    #print(p_grid.shape)
        
    plt.contour(sig_grid, rho_grid, p_grid, levels = (0.05,), colors=col, linewidths=1.5, linestyles=ls)
    plt.contourf(sig_grid, rho_grid, p_grid, levels = (0.05,1), colors=col, alpha=0.5)

def plotContour_nomodulation(m_x, sigma, outpath, col='C0', ls='dashed', fixed_mass = False):
    if (fixed_mass):
        ID_str = "m%d_lsig%6.2f_mfix.txt"%(1000*m_x, np.log10(sigma))
    else:
        ID_str = "m%d_lsig%6.2f_mprof.txt"%(1000*m_x, np.log10(sigma))
    
    if not (os.path.isfile(rootdir + outpath + "/stat_sigma_" + ID_str)):
        ID_str = "m%d_lsig%6.2f.txt"%(1000*m_x, np.log10(sigma))
    
    sig_list = np.loadtxt(rootdir + outpath + "/stat_sigma_" + ID_str)
    rho_list = np.loadtxt(rootdir + outpath + "/stat_rho_" + ID_str)
    p_grid = np.loadtxt(rootdir + outpath + "/stat_p_" + ID_str)

    #print(np.min(p_grid))
    #sig_grid, rho_grid = np.meshgrid(sig_list, rho_list)
    
    rho_min = np.zeros(len(sig_list))
    rho_max = np.zeros(len(sig_list))
    
    for i in range(len(sig_list)):
        inds = np.where(p_grid[:,i] > 0.05)[0]
        if (len(inds) == 0):
            rho_min[i] = np.nan
            rho_max[i] = np.nan
        else:
            rho_min[i] = rho_list[inds[0]]
            rho_max[i] = rho_list[inds[-1]]
        
    plt.plot(sig_list, rho_min, color=col, linewidth=1.5, linestyle=ls)
    plt.plot(sig_list, rho_max, color=col, linewidth=1.5, linestyle=ls)
    plt.fill_between(sig_list, rho_min, rho_max, color=col, alpha=0.1)
        

