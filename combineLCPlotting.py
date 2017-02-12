import numpy as np
import math as math
import cmath as cmath
import psutil as psutil
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from matplotlib import gridspec as gridspec
import argparse as argparse
import operator as operator
import warnings as warnings
import copy as copy
import time as time
import pdb
import os as os
import random

import kali.k2
import kali.s82
import kali.carma
import kali.util.mcmcviz as mcmcviz
from kali.util.mpl_settings import set_plot_params
import kali.util.triangle as triangle

import matplotlib.patches as patches

# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height

fig = plt.figure()
ax = fig.add_axes([0,0,1,1])

# axes coordinates are 0,0 is bottom left and 1,1 is upper right
p = patches.Rectangle(
    (left, bottom), width, height,
    fill=False, transform=ax.transAxes, clip_on=False
    )



def doubleMADsfromMedian(y,thresh=3.):
    # warning: this function does not check for NAs
    # nor does it address issues when 
    # more than 50% of your data have identical values
    m = np.median(y)
    abs_dev = np.abs(y - m)
    left_mad = np.median(abs_dev[y <= m])
    right_mad = np.median(abs_dev[y >= m])
    y_mad = left_mad * np.ones(len(y))
    y_mad[y > m] = right_mad
    modified_z_score = 0.6745 * abs_dev / y_mad
    modified_z_score[y == m] = 0
    return np.where(modified_z_score < thresh)[0]
    
# make new kepler flux
def kepFlux(sdsslc_r,sdsslc_g):
    sdss_g = sdsslc_g.y
    sdss_r = sdsslc_r.y
    sdss_t = sdsslc_g.t
    sdss_gerr = sdsslc_g.yerr
    
    if len(sdss_r) != len(sdss_g):
        print sdsslc_r.t[-1], sdsslc_g.t[-1]
        #catch missing g band measurements in the middle of the array 
        for i in range(0, len(sdss_t)):
            tolerance = 0.0002
            missed = np.isclose(sdsslc_r.t, sdss_t[i],tolerance)
            m = np.where(missed == True)[0]
            #print m, i
            if m != i:
                print m, i, len(sdss_t)
                sdss_g = np.insert(sdss_g, m ,0.)
                sdss_gerr = np.insert(sdss_gerr, m ,0.)
                sdss_t = np.insert(sdss_t, m ,sdsslc_r.t[m])
                print sdss_g[m], sdss_t[i], len(sdss_t),len(sdsslc_r.t)
            if len(sdss_r) == len(sdss_g):
                break
    c = sdss_r*0.8 + sdss_g*0.2
    fullr = np.where(sdss_g == 0.)[0]
    c[fullr] = sdss_r[fullr]
    print(c[fullr],sdss_r[fullr-1])
    c_err = np.zeros(len(sdss_t))
    for i in range (len(sdss_t)):
    	#c_err[i] = np.maximum(sdsslc_r.yerr[i], sdss_gerr[i])
    	
    	#compute error by adding in quadrature
    	thing1 = 0.8
    	thing2 = 0.2
    	c_err = np.sqrt( (thing1**2. * sdsslc_r.yerr**2.) + (thing2**2. * sdss_gerr**2.))

    c_t = sdss_t
    #return 0,0,0                   
    return c, c_err, c_t



#Loading SDSS g and r band flux lightcurves: object #1
#id = '220191502'
#sdss_id = '005523.82-001941.9'
#id = '220212788'
#sdss_id = '012147.73+002718.7'
#id = '220176624'
#sdss_id = '011417.18-005518.8'
#id  = '220224214'
#sdss_id = "005232.43+005123.0"
#---------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-id', '--id', type = str, default = '220173631', help = r'EPIC ID')
parser.add_argument('-sdssid', '--sdssid', type = str, default = '005102.42-010244.4', help = r'EPIC ID')
parser.set_defaults(viewer = True)
args = parser.parse_args()

id = args.id
sdss_id = args.sdssid


sdsslc_r = kali.s82.sdssLC(name = sdss_id, band = 'r')
sdsslc_g = kali.s82.sdssLC(name = sdss_id, band = 'g')
lastT = sdsslc_r.t[-1] +100.

import seaborn as sns

sns.palplot(sns.color_palette("GnBu_d"))


c, c_err, c_t = kepFlux(sdsslc_r,sdsslc_g)
#plt.scatter(c_t, c)
#plt.savefig('TestCombineLC.png')
#plt.errorbar(sdsslc_r.t, c, yerr = c_err)

plt.clf()
#sanity check: Is the lightcurve normalized? Answer : no
c_med = np.median(c)
c_norm = c/c_med
print("the average flux is %d",c_med)
#plt.show()
#plt.scatter(sdsslc_r.t, c_norm)
#plt.savefig('normedComboLC.png')
#plt.clf()

k2lc = kali.k2.k2LC(name = id, band = 'Kep', processing = 'k2sff', campaign = 'c08')


sdss_std = np.std(c_norm)
#write function to do all of this 
full_lcy = np.concatenate([c_norm, np.array(k2lc.y)])
full_lcyerr = np.concatenate([c_err, np.array(k2lc.yerr)])
full_lct = np.concatenate([c_t, np.array(k2lc.t)])
w = np.where(full_lcy > 0)[0]
#compare time attibutes in sdss and k2 objects from kali
print(len(sdsslc_g.t), len(k2lc.t),len(full_lct))
k2t = k2lc.t + 3310.
full_lct = np.concatenate([c_t, k2t])
w = np.where(full_lcy > 0)[0]
#plt.plot(full_lct[w], full_lcy[w])
#plt.savefig('combined.png')
#plt.clf()
#match mask index
#mmatch errors and deal with outliers
y = full_lcy[w]
yt = full_lct[w]
z = doubleMADsfromMedian(y )
#plt.plot(yt[z], y[z])
#plt.savefig('combined_noOutliers.png')


#reset K2LC
k2lc.y = y

k2lc.yerr = full_lcyerr[w]

k2lc.t = yt

k2lc.mask = np.zeros(len(k2lc.y))
k2lc.mask[z] = 1
k2lc.cadence = np.arange(0,len(k2lc.y))
plt.scatter(k2lc.t , k2lc.y*k2lc.mask)
plt.savefig('combined_noOutliers_k2lc.png')

#correct other timescale params 
#sigma = maxSigma*var(lc)
k2lc.startT = 0.
k2lc._dt = 0.5 ## Increment between epochs.
k2lc._mindt = 0.02
k2lc._maxdt = 3010.
k2lc._T =k2lc.t[-1] - k2lc.t[0] ## Total duration of the light curve.
k2lc._numCadences = len(k2lc.y)
#self.cadence = 0.02

from matplotlib import gridspec
#sns.set_style("whitegrid")
#f, (ax, ax2) = plt.subplots(1, 2, sharey=True)
f = plt.figure()
gs = gridspec.GridSpec(1, 3)

ax = f.add_subplot(gs[0,:2])
ax.scatter(k2lc.t[z],k2lc.y[z], marker = '+',color = '#6e016b', label=str('SDSS '+sdss_id))
sdsserr = np.zeros(np.shape(c_norm))
sdsserr[:len(c_err)]= c_err/c_med
ax.fill_between(c_t, c_norm-sdsserr, c_norm+sdsserr, color='b', alpha=0.2)

ax.legend(loc="upper left", prop={'size':15})
ax.set_xlim(c_t[0], c_t[-1])  
ax.set_ylim(0.6, 1.3)  
#-------------------------------------------
ax2 = f.add_subplot(gs[0,2])
ax2.fill_between(k2lc.t[z],k2lc.y[z]-k2lc.yerr[z], k2lc.y[z]+k2lc.yerr[z], color='c', alpha=0.2)
ax2.scatter(k2lc.t[z],k2lc.y[z],  marker = "+",color =  '#ff7f00', label=str("K2 "+id))
ax2.set_xlim(k2t[0], k2t[-1])  
ax2.set_ylim(0.6, 1.3)  

ax.spines['right'].set_visible(False)
#ax2.spines['left'].set_visible(False)
#
ax2.tick_params(labelleft='off')  # don't put tick labels at the top

x = [k2t[0], k2t[0]+20, k2t[0]+40, k2t[0]+60,k2t[0]+80]
labels = [str(int(k2t[0])), str(int(k2t[0]+20)), str(int(k2t[0]+40)), str(int(k2t[0]+60)), str(int(k2t[0]+80))]

ax2.set_xticks(x)
ax2.set_xticklabels(labels,fontname = "Times New Roman")


#ax2.xaxis.tick_bottom()
ax2.legend(loc="upper left", prop={'size':15})

d = .020  # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)


ax.plot((1 - d, 1 + d), (-d, +d), **kwargs) # top-right diagona
ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes

ax2.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal l
#f.subplots_adjust(hspace=0.1)
#gs.update(wspace=0.5, hspace=0.05)
# What's cool about this is that now if we vary the distance between
# ax and ax2 via f.subplots_adjust(hspace=...) or plt.subplot_tool(),
f.subplots_adjust(hspace=0.01)
f.tight_layout()
f.savefig(id+'ExtendedLC.png', dpi=400)




