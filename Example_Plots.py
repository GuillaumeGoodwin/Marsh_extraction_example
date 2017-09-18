"""
Example_Plots.py

This script produces nice plots to visualise your results

"""


#------------------------------------------------------------------
#0. Set up display environment in putty if you are working on a terminal with no graphical interface.
import matplotlib
matplotlib.use('Agg')

#----------------------------------------------------------------
#1. Load useful Python packages
import os
import sys

import numpy as np
import functools
import math as mt
import cmath
import scipy as sp
import scipy.stats as stats
from datetime import datetime
import cPickle
from pylab import *
import functools
import itertools as itt
from osgeo import gdal, osr
from osgeo import gdal, gdalconst
from osgeo.gdalconst import *
from copy import copy
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import *
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap, cm
from matplotlib.patches import Rectangle
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

#------------------------------------------------------------------
# Import the marsh-finding functions
from Example_functions import ENVI_raster_binary_to_2d_array
from Example_functions import ENVI_raster_binary_from_2d_array
from Example_functions import add_subplot_axes
from Example_functions import Distribution
from Example_functions import define_search_space


# This is a function straight out of the Python cookbook: https://stackoverflow.com/questions/26911898/matplotlib-curve-with-arrow-ticks
def add_arrow_to_line2D(
    axes, line, arrow_locs=[0.2, 0.4, 0.6, 0.8],
    arrowstyle='-|>', arrowsize=1, transform=None):
    """
    Add arrows to a matplotlib.lines.Line2D at selected locations.

    Parameters:
    -----------
    axes: 
    line: Line2D object as returned by plot command
    arrow_locs: list of locations where to insert arrows, % of total length
    arrowstyle: style of the arrow
    arrowsize: size of the arrow
    transform: a matplotlib transform instance, default to data coordinates

    Returns:
    --------
    arrows: list of arrows
    """
    if not isinstance(line, mlines.Line2D):
        raise ValueError("expected a matplotlib.lines.Line2D object")
    x, y = line.get_xdata(), line.get_ydata()

    arrow_kw = {
        "arrowstyle": arrowstyle,
        "mutation_scale": 10 * arrowsize,
    }

    color = line.get_color()
    use_multicolor_lines = isinstance(color, np.ndarray)
    if use_multicolor_lines:
        raise NotImplementedError("multicolor lines not supported")
    else:
        arrow_kw['color'] = color

    linewidth = line.get_linewidth()
    if isinstance(linewidth, np.ndarray):
        raise NotImplementedError("multiwidth lines not supported")
    else:
        arrow_kw['linewidth'] = linewidth

    if transform is None:
        transform = axes.transData

    arrows = []
    for loc in arrow_locs:
        s = np.cumsum(np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2))
        n = np.searchsorted(s, s[-1] * loc)
        arrow_tail = (x[n], y[n])
        arrow_head = (np.mean(x[n:n + 2]), np.mean(y[n:n + 2]))
        p = mpatches.FancyArrowPatch(
            arrow_tail, arrow_head, transform=transform,
            **arrow_kw)
        axes.add_patch(p)
        arrows.append(p)
    return arrows



# This is a function that makes an outline out of a raster where each object has a single given value

def Outline (Raster, Outline_value):

    P1 = np.where(Raster[:,1:] != Raster[:,:-1]) 
    Raster[P1] = Outline_value           

    P2 = np.where(Raster[1:,:] != Raster[:-1,:])
    Raster[P2] = Outline_value
    
    return Raster






#------------------------------------------------------------------
#2. Set up the important variables

#Select site names. Simply add a site in the list to analyse multiple sites simultaneously.
Sites = ["FEL"]

# Set the value for empty DEM cells
Nodata_value = -9999




#------------------------------------------------------------------
#3. Start Plotting

#Plot 1: Draw the platform on a DEM, superimposed on a hillshade
for site in Sites:
    fig=plt.figure(1, facecolor='White',figsize=[10,10])
    ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1, axisbg='white')

    # Set up the fonts and stuff
    matplotlib.rc('xtick', labelsize=9) 
    matplotlib.rc('ytick', labelsize=9)

    # Name the axes
    ax1.set_xlabel('x (m)', fontsize = 12)
    ax1.set_ylabel('y (m)', fontsize = 12)

    # Load the relevant data
    HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array ("Input/%s_hs_clip.bil" % (site), site)
    DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array ("Input/%s_DEM_clip.bil" % (site), site)
    Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array ("Output/%s_Marsh.bil" % (site), site)

    
    # Make a map!
    Platform_mask = np.ma.masked_where(Platform <=0, Platform)
    Platform_mask[Platform_mask>0] = DEM[Platform_mask>0]

    Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
    Map_DEM = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(DEM[DEM!=Nodata_value]), vmax = np.amax(DEM), alpha = 0.5)
    Map_Marsh = ax1.imshow(Platform_mask, interpolation='None', cmap=plt.cm.gist_earth, vmin=np.amin(DEM[DEM!=Nodata_value]), vmax=np.amax(DEM), alpha = 0.5)
    


plt.savefig('Output/Figure_1.png')




#Plot 2: Draw the marsh outline on a DEM, superimposed on a hillshade
for site in Sites:
    fig=plt.figure(2, facecolor='White',figsize=[10,10])
    ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1, axisbg='white')

    # Set up the fonts and stuff
    matplotlib.rc('xtick', labelsize=9) 
    matplotlib.rc('ytick', labelsize=9)

    # Name the axes
    ax1.set_xlabel('x (m)', fontsize = 12)
    ax1.set_ylabel('y (m)', fontsize = 12)

    # Load the relevant data
    HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array ("Input/%s_hs_clip.bil" % (site), site)
    DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array ("Input/%s_DEM_clip.bil" % (site), site)
    Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array ("Output/%s_Marsh.bil" % (site), site)

    # Outline the marsh
    Platform[Platform > 0] = 1
    Marsh_outline = Outline (Platform,2)

    
    # Make a map!
    Outline_mask = np.ma.masked_where(Marsh_outline <=1, Marsh_outline)


    Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
    Map_DEM = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(DEM[DEM!=Nodata_value]), vmax = np.amax(DEM), alpha = 0.5)
    
    Map_Marsh = ax1.imshow(Outline_mask, interpolation='None', cmap=plt.cm.Reds, vmin = 0, vmax = 2, alpha = 1)
    

plt.savefig('Output/Figure_2.png')



#Plot 3: Draw the confusion map, superimposed on a hillshade
for site in Sites:
    fig=plt.figure(1, facecolor='White',figsize=[10,10])
    ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1, axisbg='white')

    # Set up the fonts and stuff
    matplotlib.rc('xtick', labelsize=9) 
    matplotlib.rc('ytick', labelsize=9)

    # Name the axes
    ax1.set_xlabel('x (m)', fontsize = 12)
    ax1.set_ylabel('y (m)', fontsize = 12)

    # Load the relevant data
    HS, post_HS, envidata_HS = ENVI_raster_binary_to_2d_array ("Input/%s_hs_clip.bil" % (site), site)
    DEM, post_DEM, envidata_DEM = ENVI_raster_binary_to_2d_array ("Input/%s_DEM_clip.bil" % (site), site)
    Platform, post_Platform, envidata_Platform = ENVI_raster_binary_to_2d_array ("Output/%s_Marsh.bil" % (site), site)
    Reference, post_Reference, envidata_Reference = ENVI_raster_binary_to_2d_array ("Input/%s_ref_DEM_clip.bil" % (site), site)

    # Outline the reference
    Reference[Reference > 0] = 1
    Ref_outline = Outline (Reference,2)

    
    # Make a map!
    Outline_mask = np.ma.masked_where(Ref_outline <=1, Ref_outline)
    
    
    # Make a map!
    Platform_mask = np.ma.masked_where(Platform <=0, Platform)
    Platform_mask[Platform_mask>0] = DEM[Platform_mask>0]

    Map_HS = ax1.imshow(HS, interpolation='None', cmap=plt.cm.gist_gray, vmin = 100, vmax = 200)
    Map_DEM = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.gist_gray, vmin = np.amin(DEM[DEM!=Nodata_value]), vmax = np.amax(DEM), alpha = 0.5)
    Map_Marsh = ax1.imshow(Platform_mask, interpolation='None', cmap=plt.cm.gist_earth, vmin=np.amin(DEM[DEM!=Nodata_value]), vmax=np.amax(DEM), alpha = 0.5)
    
    Map_Marsh = ax1.imshow(Outline_mask, interpolation='None', cmap=plt.cm.Reds, vmin = 0, vmax = 2, alpha = 1)

plt.savefig('Output/Figure_3.png')

