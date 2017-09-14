"""
This is a Python script. It does quite a lot of useful things:

1. If you give it a DEM, it will prepare it for analysis

2. It will then find the marsh platforms and scarps inside your DEM

"""

#------------------------------------------------------------------
#0. Set up display environment in putty if you are working on a terminal with no graphical interface.
import matplotlib
matplotlib.use('Agg')

#------------------------------------------------------------------
#1. Load useful Python packages
import os
import sys
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
import functools
import math as mt
import cmath
import scipy as sp
import scipy.stats as stats
from datetime import datetime
import cPickle
from matplotlib import cm
from pylab import *
import functools
import matplotlib.ticker as tk
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import *
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import itertools as itt
from osgeo import gdal, osr
import matplotlib.ticker as tk
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import *
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from osgeo import gdal, gdalconst
from osgeo.gdalconst import *
import fileinput, sys

#------------------------------------------------------------------
# Import the costum-made marsh-finding functions
from Example_functions import Open_tide_stats
from Example_functions import MARSH_ID
from Example_functions import Confusion
from Example_functions import METRIC1_DEM
from Example_functions import METRIC2_DEM
from Example_functions import METRIC3_DEM
from Example_functions import ENVI_raster_binary_to_2d_array
from Example_functions import ENVI_raster_binary_from_2d_array



#------------------------------------------------------------------
#2. Set up the important variables

#Select site names. Simply add a site in the list to analyse multiple sites simultaneously.
Sites = ["FEL"]

# Set the value for empty DEM cells
Nodata_value = -9999

# Set the value for optimised detection parameters
opt1 = -2.0
opt2 = 1
opt3 = 1


#------------------------------------------------------------------
#3. Run the Marsh identification script

for site in Sites:
    print "Loading input data"
    print " Loading DEM"
    DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Input/%s_DEM_clip.bil" % (site), site)
    print " Loading Slopes"
    Slope, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array ("Input/%s_slope_clip.bil" % (site), site)
    print " Loading Curvature"
    Curvature, post_Curvature, envidata_Curvature =  ENVI_raster_binary_to_2d_array ("Input/%s_curvature_clip.bil" % (site), site)

    
    print "Identifying the platform and scarps"
    DEM_work = np.copy(DEM)
    Search_space, Scarps, Platform = MARSH_ID(DEM, Slope, Curvature, Nodata_value, opt1, opt2, opt3)
    Platform_work = np.copy(Platform)
    Scarps[Scarps == 0] = Nodata_value



    print "Saving marsh features"
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/%s_Search_space.bil" % (site), post_DEM, Search_space)
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/%s_Scarps.bil" % (site), post_DEM, Scarps)
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/%s_Marsh.bil" % (site), post_DEM, Platform)



    print " Loading Detected Marsh"
    Platform_work, post_Platform, envidata_Platform =  ENVI_raster_binary_to_2d_array ("Output/%s_Marsh.bil" % (gauge,gauge,res), gauge)

    
    
    print "Measuring performances"
    Confusion_matrix, Performance, Metrix = Confusion (Platform_work, Reference, Nodata_value)
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_Platform, "Output/%s/%s_%s_Confusion_DEM.bil" % (gauge, gauge,res), post_Platform, Confusion_matrix)

    
    cPickle.dump(Performance,open("Output/%s/%s_%s_Performance.pkl" % (gauge,gauge,res), "wb"))
    cPickle.dump(Metrix,open("Output/%s/%s_%s_Metrix.pkl" % (gauge,gauge,res), "wb"))


