#----------------------------------------------------------------
#1. Load useful Python packages
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal, osr, gdalconst
from osgeo.gdalconst import *
import cPickle




#---------------------------------------------------------------
# This function generates data distributions out of a 2D raster
def Distribution(Data2D, Nodata_value):
    

    Data1D = Data2D.ravel()

    Max_distribution = max(Data1D)    
    if len(Data1D[Data1D>Nodata_value]) == 0:
        Min_distribution = -1
    else:
        Min_distribution = min(Data1D[Data1D>Nodata_value])
    
    bin_size = (Max_distribution - Min_distribution) / 100
    
    X_values = np.arange(Min_distribution, Max_distribution, bin_size)
    
    
    hist, bins = np.histogram (Data1D, X_values, density=True)
    hist=hist/sum(hist)
    bins=bins[:-1]
   
    
    return bins,hist




#---------------------------------------------------------------
# This function opens tidal stats data
def Open_tide_stats (file_location, gauge):
    print 'Loading tidal statistics of' + file_location + ' %s' % (gauge)


    with open (file_location + "Metric1.pkl", 'rb') as input_file:
        Metric1 = cPickle.load(input_file)
    with open (file_location + "Metric2.pkl", 'rb') as input_file:
        Metric2 = cPickle.load(input_file)
    with open (file_location + "Metric3.pkl", 'rb') as input_file:
        Metric3 = cPickle.load(input_file)
    with open (file_location + "Subsample.pkl", 'rb') as input_file:
        Subsample = cPickle.load(input_file)

    return Metric1, Metric2, Metric3, Subsample




#---------------------------------------------------------------
# This function calculates the pdf and cdf of the DEM

def METRIC1_DEM (DEM, Zstep):
    print 'Calculating pdf and cdf'
    Z_data = DEM.ravel()
    Z_value = np.arange(min(Z_data[Z_data>-100]), max(Z_data), Zstep)
    METRIC1 = np.zeros((3, len(Z_value)-1),dtype=np.float)

    # This is the PDF
    hist,bins = np.histogram (Z_data,Z_value,density=True); hist=hist/sum(hist); bins=bins[:-1]

    # And this is the CDF
    cdf=[0]
    for i in range(1,len(hist)-1):
        if hist[i]>=5*hist[i-1] and hist[i]>=5*hist[i+1]:
            hist[i]=0
        cdf.append(hist[i]+cdf[i-1])
    cdf.append(hist[-1]+cdf[-2])
    METRIC1[0]=bins; METRIC1[1]=hist; METRIC1[2]=cdf

    return METRIC1



#---------------------------------------------------------------
# This function replaces the values in the DEM with the calculated values of Tidal percentiles
# It then calculates the statistics of area, perimeter, total channel length, mean (and std) of elevation/slope, curvature

# Call this function with DEM_new, DEM_MarshChannels, DEM_TFChannels, Average_Z, Stdev_Z, Area, Perimeter, Channel_length = METRIC2_DEM (DEM, Metric2, Channels)
def METRIC2_DEM (DEM, Metric2, Channels):
    DEM_new = np.copy(DEM)

    Average_Z=np.zeros(len(Metric2), dtype=np.float)
    Stdev_Z=np.zeros(len(Metric2), dtype=np.float)
    Area=np.zeros(len(Metric2), dtype=np.float)
    Perimeter=np.zeros(len(Metric2), dtype=np.float)
    Channel_length=np.zeros(2, dtype=np.float)

    for i in range(len(Metric2)):
        Threshold = np.mean(Metric2[i])

        Average_Z[i] = np.mean (DEM_new [DEM_new<=Threshold]) # Calculate the average elevation of elements below threshold
        Stdev_Z[i] = np.std (DEM_new [DEM_new<=Threshold]) # Calculate the average elevation of elements below threshold

        DEM_new [DEM_new<=Threshold] = (i+1)*100 # Replace elevation values with category value (100: Subtidal ; 200: [MLWS -> MLWN]; 300: [MLWN -> MHWN]; 400: [MHWN -> MHWS])

        Area[i] = DEM_new[DEM_new == (i+1)*100].size
        Perimeter[i] = np.sum(DEM_new[:,1:] != DEM_new[:,:-1]) + np.sum(DEM_new[1:,:] != DEM_new[:-1,:])


    Channels[Channels<200] = 0 # filtering

    DEM_MarshChannels= DEM_new
    DEM_TFChannels= DEM_new

    DEM_MarshChannels [DEM_MarshChannels != 400] = 0
    DEM_TFChannels [DEM_TFChannels != 400] = 0

    DEM_MarshChannels = DEM_MarshChannels * Channels /400
    DEM_TFChannels = DEM_TFChannels * Channels /300

    print max(sum(np.isnan (DEM_MarshChannels))),', ', max(sum(np.isnan (DEM_TFChannels)))

    Channel_length[0] = sum(DEM_TFChannels [DEM_TFChannels > 0] > 0) # Tidal flats
    Channel_length[1] = sum(DEM_MarshChannels [DEM_MarshChannels > 0] > 0) # Marsh


    return DEM_new, DEM_MarshChannels, DEM_TFChannels, Average_Z, Stdev_Z, Area, Perimeter, Channel_length


#---------------------------------------------------------------
# This function replaces the values in the DEM with the calculated values of Flooding parameters

# Call this function with DEM_end = METRIC3_DEM (DEM, Metric3)

def METRIC3_DEM (DEM, Metric3):
    DEM_work = np.copy(DEM)

    Index1 = np.where (Metric3[1]>12)
    Threshold1 = min(Metric3[0, min(Index1)])

    Index2 = np.where (Metric3[1]>24)
    Threshold2 = min(Metric3[0, min(Index2)])

    Index3 = np.where (Metric3[1]>72)
    Threshold3 = min(Metric3[0, min(Index3)])

    print Threshold1, Threshold2, Threshold3


    DEM_work [DEM_work <= Threshold1] = 100
    DEM_work [DEM_work <= Threshold2] = 200
    DEM_work [DEM_work <= Threshold3] = 300


    return DEM_work





#-----------------------------------------------------------------------------------------------------
# This functions calculates local slope using the maximum slope method
# http://www.onlinegeographer.com/slope/Dunn_hickey_slope.pdf, equation 7

# It takes as input:
# 1/ The DEM array

# It returns:
# 1/ A Slope array


def maximum_slope (DEM):
    Height = len(DEM_work); Width = len(DEM_work[0,:])
    Slope_max = np.zeros((Height,Width), dtype=np.float)

    for i in range(1,len(DEM_work)-1): # rows
        for j in range(1,len(DEM_work[0,:])-1):  # cols
            # Here are the cells we consider
            Cells = np.zeros(9, dtype=np.float)
            Slopes = np.zeros(9, dtype=np.float)
            Cells[0]=DEM_work[i,j]; Cells[1]=DEM_work[i-1,j]; Cells[2]=DEM_work[i-1,j+1]; Cells[3]=DEM_work[i,j+1]; Cells[4]=DEM_work[i+1,j+1]; Cells[5]=DEM_work[i+1,j]; Cells[6]=DEM_work[i+1,j-1]; Cells[7]=DEM_work[i,j-1]; Cells[8]=DEM_work[i-1,j-1]

            # Calculate Slopes
            for a in range(1,len(Cells)):
                Slopes[a]= (Cells[a]-Cells[0])/1

            # Largest Slope is the slope value for i,j
            Slope_max[i,j]=max(Slopes)


    return Slope_max





#-----------------------------------------------------------------------------------------------------
# This functions defines a search space within an input raster

# It takes as input:
# 1/ The data array
# 2/ The slope array

# It returns:
# 1/ A search space within the data array

def define_search_space (DEM, Slope, Nodata_value,opt):
    print 'Choosing a holiday destination ...'
    Height = len(DEM); Width = len(DEM[0,:])
    Search_space = np.zeros((Height,Width), dtype=np.float)

    # We calculate the relative relief of the DEM to have values of elevation between 0 and 1
    Relief = DEM-np.amin(DEM[DEM > Nodata_value])
    Rel_relief = Relief/np.amax(Relief)
    Rel_relief[DEM == Nodata_value] = Nodata_value

    # We then do the same thing for slope
    Rel_slope = Slope/np.amax(Slope)
    Rel_slope[Slope == Nodata_value] = Nodata_value

    # We then multiply these new relative relief and slope arrays and biologically name them "Crossover"
    Crossover = Rel_relief * Rel_slope
    Crossover[DEM == Nodata_value] = Nodata_value

    # We make a curve of the frequency of values in this Crossover
    # That curve should look like a decreasing exponential function
    data = Crossover.ravel(); data = data[data>0]
    step = (max(data) - min(data)) / 100
    value = np.arange(min(data), max(data), step)
    hist, bins = np.histogram (data, value, density=True)
    hist=hist/sum(hist); bins=bins[:-1]

    # We now find the slope of that curve
    hist_der = np.zeros(len(hist), dtype = np.float)
    for j in range(1, len(hist), 1):
        hist_der[j] = (hist[j]-hist[j-1])/step

    # If the slope gets above the -1 threshold, now that we have hit the closest point to the origin.
    # We call it the inflexion point even though it's not really an inflexion point.
    for j in range(1, len(hist)-1, 1):
        if hist_der[j] < opt and hist_der[j+1] >= opt:
            Inflexion_point = bins[j]

    # Points within the search space should have a Crossover value above the inflexion point
    Search = np.where(Crossover > Inflexion_point)
    Search_space[Search] = 1

    # We get rid of the borders of the DEM because otherwise it will be difficult to work with the smaller slope array
    Search_space[0,:] = 0; Search_space[Height-1,:] = 0; Search_space[:,0] = 0; Search_space[:,Width-1] = 0

    # And update the search locations for the shaved edges
    Search = np.where(Search_space == 1)

    # If this happens, your landscape is weird
    if np.amax(Search_space) == 0:
        print
        print " ... Your search space is empty! Are you sure there's a marsh platform here?"
        print
        STOP
        
    
    
    
  
    #fig=plt.figure(3, facecolor='White',figsize=[4.7,4])
    #ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1, axisbg='white')
    #ax1.plot(bins, hist)
    #plt.savefig('Output/Paper/0_Fig3.png')
    #STOP

    return Search_space, Crossover, bins, hist, Inflexion_point


#-----------------------------------------------------------------------------------------------------
# This functions makes a kernel in an input raster

# It takes as input:
# 1/ The data array
# 2/ The desired kernel size (the total length or height of it!)
# 3/ The coordinates of the kernel centre

# It returns:
# 1/ A kernel, which is a lower ranking officer than Hollywood would have you believe

def kernel (array, kernel_size, x_centre, y_centre):
    if (-1)**kernel_size < 0:
        X_to_0 = x_centre
        X_to_End = len(array)-x_centre
        Y_to_0 = y_centre
        Y_to_End = len(array[0,:])-y_centre

        Lim_left = x_centre - min(np.floor(kernel_size/2), X_to_0)
        Lim_right = x_centre + min(np.floor(kernel_size/2)+1, X_to_End)
        Lim_top = y_centre - min(np.floor(kernel_size/2), Y_to_0)
        Lim_bottom = y_centre + min(np.floor(kernel_size/2)+1, Y_to_End)

        kernel = array [Lim_left:Lim_right, Lim_top:Lim_bottom]

    else:
        print
        print " ... WARNING: you need to choose an odd kernel size, buddy"
        print
        pass

    return kernel


#-----------------------------------------------------------------------------------------------------
# This functions detects and flags local peaks of a data array

# It takes as input:
# 1/ The data array
# 2/ The desired search space (0 = don't search; 1 = search)

# It returns:
# 1/ An array where local peaks have a value of 1 (all other values are 0)
# 2/ A copy of the data array where the locations of the peaks are 0.

def peak_flag (Slope, Search_space, Order):
    print 'Preparing the Kilimanjaro expedition ...'
    Slope_copy = np.copy(Slope) # the copy of the initial data array
    Search = np.where(Search_space == 1) # the searched locations
    Peaks = np.zeros((len(Slope),len(Slope[0,:])),dtype = np.float)

    for i in range(len(Search[0])):
        x=Search[0][i]; y=Search[1][i] # coordinates of the kernel's centre
        Kernel_slope = kernel (Slope, 3, x, y)
        Kernel_search = kernel(Search_space, 3, x, y)

        # if the centre of the kernel is its maximum and is not an isolated point
        if Kernel_slope[1,1] == np.amax(Kernel_slope) and np.amax(Kernel_search[Kernel_search<=Kernel_search[1,1]] > 0):
            Peaks[x,y] = Order # The kernel centre becomes a local peak
            Slope_copy[x,y] = 0 # The slope of the modified data array drops to 0

    return Peaks, Slope_copy



#-----------------------------------------------------------------------------------------------------
# This functions detects the 2 or more points from which ridges are initiated. We call them ridge starters

# It takes as input:
# 1/ The data array
# 2/ The desired search space (0 = don't search; 1 = search)
# 3/ The array containing the location of local peaks

# It returns:
# 1/ An array where local peaks are 1 and ridge starters are 2 (all other values are 0)
# 2/ A copy of the data array where the locations of the peaks and ridge starters are 0.


def initiate_ridge (Slope, Search_space, Peaks, Order):
    print ' ... Rolling off Mt. Kilimanjaro ...'
    Slope_copy = np.copy(Slope) # the copy of the initial data array
    Search = np.where(Search_space == 1) # the searched locations
    Search_peaks = np.where(Peaks == Order-1) # the searched locations where the peaks are
    Ridges = np.copy(Peaks)

    # Define Kernels
    for i in range(len(Search_peaks[0])):
        x=Search_peaks[0][i]; y=Search_peaks[1][i] # coordinates of the kernel's centre
        Kernel_slope = kernel (Slope, 3, x, y)
        Kernel_slope_copy = kernel (Slope_copy, 3, x, y)
        Kernel_ridges = kernel (Ridges, 3, x, y)
        Kernel_search = kernel (Search_space, 3, x, y)

        # 1/ If there are no other peaks, we have two ridge starters
        if np.count_nonzero(Kernel_ridges) == 1:
            Ridge_starter1 = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
            X1=Ridge_starter1[0][0]; Y1=Ridge_starter1[1][0]

            # if it is within the initial search space
            if Search_space[x+X1-1, y+Y1-1] != 0:
                Ridges[x+X1-1, y+Y1-1] = Order
                Slope_copy[x+X1-1, y+Y1-1] = 0

                # Look for a second ridge starter
                Ridge_starter2 = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
                X2=Ridge_starter2[0][0]; Y2=Ridge_starter2[1][0]
                Distance = np.sqrt((X2-X1)**2+(Y2-Y1)**2)

                # if it is within the initial search space AND not next to the first ridge starter
                if Search_space[x+X2-1, y+Y2-1] != 0 and Distance > np.sqrt(2):
                    Ridges[x+X2-1, y+Y2-1] = Order
                    Slope_copy[x+X2-1, y+Y2-1] = 0

                # Otherwise, look for second ridge starter elsewhere in the kernel
                elif Search_space[x+X2-1, y+Y2-1] != 0 and Distance <= np.sqrt(2):
                    for j in np.arange(0,9,1):
                        Kernel_slope_copy[X2, Y2] = 0

                        Ridge_starter2 = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
                        X2=Ridge_starter2[0][0]; Y2=Ridge_starter2[1][0]
                        Distance = np.sqrt((X2-X1)**2+(Y2-Y1)**2)

                        if Search_space[x+X2-1, y+Y2-1] != 0 and Distance > np.sqrt(2):
                            Ridges[x+X2-1, y+Y2-1] = Order
                            Slope_copy[x+X2-1, y+Y2-1] = 0
                            break


        # 2/ If there are two peaks, we have one ridge starter
        elif np.count_nonzero(Kernel_ridges) == 2:
            Ridge_starter1 = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
            X1=Ridge_starter1[0][0]; Y1=Ridge_starter1[1][0]

            # if it is within the initial search space
            if Search_space[x+X1-1, y+Y1-1] != 0:
                Ridges[x+X1-1, y+Y1-1] = Order
                Slope_copy[x+X1-1, y+Y1-1] = 0

    return Ridges, Slope_copy






#-----------------------------------------------------------------------------------------------------
# This functions continues the ridges from the ridge starters. It incorporates tidal range to filter out noise

# It takes as input:
# 1/ The data array
# 2/ The desired search space (0 = don't search; 1 = search)
# 3/ The array containing the location of local peaks
# 4/ The Tidal ranges metric
# 5/ The ridge Order

# It returns:
# 1/ An array where local peaks are 1 and ridge starters are 2 (all other values are 0)
# 2/ A copy of the data array where the locations of the peaks and ridge starters are 0.


def Continue_ridge (DEM, Slope, Search_space, Peaks, Order, Tidal_ranges):
    #print ' ... Rolling down ...'

    DEM_copy = np.copy(DEM) # the copy of the initial DEM array
    Slope_copy = np.copy(Slope) # the copy of the initial slope array
    Search = np.where(Search_space == 1) # the searched locations
    Search_peaks = np.where(Peaks == Order-1) # the searched locations where the peaks are
    Ridges = np.copy(Peaks)

    # Define Kernels
    for i in range(len(Search_peaks[0])):
        x=Search_peaks[0][i]; y=Search_peaks[1][i] # coordinates of the kernel's centre

        Kernel_slope = kernel (Slope, 3, x, y)
        Kernel_slope_copy = kernel (Slope_copy, 3, x, y)
        Kernel_ridges = kernel (Ridges, 3, x, y)
        Kernel_search = kernel (Search_space, 3, x, y)

        # Count the number of nonzero points in the kernel of the ridge array
        Ridge_count = np.count_nonzero(Kernel_ridges)

        # If there are only the 2 previous ridge points, draw a third point that is far enough from the previous point
        if Ridge_count == 2:
            New_point = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
            X=New_point[0][0]; Y=New_point[1][0]
            Grandad_point = np.where (Kernel_ridges == Order-2)
            Xgd=Grandad_point[0][0]; Ygd=Grandad_point[1][0]
            Distance = np.sqrt((X-Xgd)**2+(Y-Ygd)**2)

            if Search_space[x+X-1, y+Y-1] != 0 and Distance > np.sqrt(2):
                Ridges[x+X-1, y+Y-1] = Order
                Slope_copy[x+X-1, y+Y-1] = 0

            elif Search_space[x+X-1, y+Y-1] != 0 and Distance <= np.sqrt(2):
                for j in np.arange(0,9,1):
                    Kernel_slope_copy[X, Y] = 0

                    New_point = np.where (Kernel_slope_copy == np.amax (Kernel_slope_copy))
                    X=New_point[0][0]; Y=New_point[1][0]
                    Distance = np.sqrt((X-Xgd)**2+(Y-Ygd)**2)

                    if Search_space[x+X-1, y+Y-1] != 0 and Distance > np.sqrt(2):
                        Ridges[x+X-1, y+Y-1] = Order
                        Slope_copy[x+X-1, y+Y-1] = 0
                        break

    return Ridges, Slope_copy




#-----------------------------------------------------------------------------------------------------
# This functions cleans up the ridges

# It takes as input:
# 1/ The ridge array
# 2/ The elevation array
# 3/ The slope array
# $/ the tidal range

# It returns:
# 1/ A ridge array cleansed of unnecessary ridges


def Clean_ridges (Peaks, DEM, Slope, Tidal_metric, Nodata_value,opt):
    print "Cleaning up: I want to break free ..."
    DEM_copy = np.copy(DEM)
    DEM_copy[DEM_copy==Nodata_value] = 0

    TR = np.mean(Tidal_metric[3])-np.mean(Tidal_metric[0])
    Search_ridge = np.where (Peaks != 0)


    """# We calculate the relative relief of the DEM to have values of elevation between 0 and 1
    Relief = DEM-np.amin(DEM[DEM > Nodata_value])
    Rel_relief = Relief/np.amax(Relief)
    Rel_relief[DEM == Nodata_value] = Nodata_value

    # We then do the same thing for slope
    Rel_slope = Slope/np.amax(Slope)
    Rel_slope[Slope == Nodata_value] = Nodata_value

    # We then multiply these new relative relief and slope arrays and biologically name them "Crossover"
    # We apply this only to ridges
    Crossover = Rel_relief * Rel_slope
    Crossover[DEM == Nodata_value] = Nodata_value
    Crossover[Peaks <= 0] = Nodata_value

    print 'High relief', np.amax(Rel_relief)
    print 'Low relief', np.amin(Rel_relief)
    print 'High slope', np.amax(Rel_slope)
    print 'Low slope', np.amin(Rel_slope)
    print 'High X', np.amax(Crossover)
    print 'Low X', np.amin(Crossover)"""



    """# We now find the slope of that curve
    hist_der = np.zeros(len(hist), dtype = np.float)
    for j in range(1, len(hist), 1):
        hist_der[j] = (hist[j]-hist[j-1])/step

    # If the slope gets above the -1 threshold, now that we have hit the closest point to the origin.
    # We call it the inflexion point even though it's not really an inflexion point.
    for j in range(1, len(hist)-1, 1):
        if hist_der[j] < -1 and hist_der[j+1] >= -1:
            Inflexion_point = bins[j]

    # Points within the search space should have a Crossover value above the inflexion point
    Search = np.where(Crossover > Inflexion_point)
    Search_space[Search] = 1"""




    print " ... What a relief ..."

    Cutoff = np.percentile(DEM_copy,75)
    Threshold = np.amax(DEM_copy[DEM_copy<Cutoff])
    DEM_copy[DEM_copy>Threshold]=Threshold

    for i in range(len(Search_ridge[0])):
        x=Search_ridge[0][i]; y=Search_ridge[1][i] # coordinates of the kernel's centre
        Kernel_DEM = kernel (DEM_copy, 9, x, y)
        Kernel_DEM[Kernel_DEM==Nodata_value]=0
        Kernel_relief = Kernel_DEM - np.amin(Kernel_DEM)
        Kernel_slope = kernel(Slope, 9, x, y)

        #KEEP THIS
        if np.amax(Kernel_DEM)/Threshold < 0.6:
            Peaks[x,y] = 0



    print " ... Shave the stubble ..."
    Search_ridge = np.where (Peaks != 0)
    for i in range(len(Search_ridge[0])):
        x=Search_ridge[0][i]; y=Search_ridge[1][i] # coordinates of the kernel's centre
        Kernel_ridges = kernel (Peaks, 9, x, y)
        # If there aren't at least 8 ridge points in the neighbourhood of 10 by 10
        if np.count_nonzero(Kernel_ridges) < 8:
            Peaks[x,y] = 0


    return Peaks






#-----------------------------------------------------------------------------------------------------
# This functions fills areas above the steep bits up the raster

# It takes as input:
# 1/ The ridge array
# 2/ The DEM
# 3/ Tidal properties

# It returns:
# 1/ An array with the marsh bits

def Fill_high_ground (DEM, Peaks, Tidal_metric, Nodata_value,opt):
    print "Paint me a platform ..."
    TR = np.mean(Tidal_metric[3])-np.mean(Tidal_metric[0])
    DEM_copy = np.copy(DEM)
    Marsh = np.zeros((len(DEM), len(DEM[0,:])), dtype = np.float)

    print " ... Start close to your sketch lines ..."
    Counter = 1
    Search_ridges = np.where (Peaks > 0)
    for i in range(len(Search_ridges[0])):
        x=Search_ridges[0][i]; y=Search_ridges[1][i]
        Kernel_ridges = kernel (Peaks, 3, x, y)
        Kernel_DEM = kernel (DEM, 3, x, y)

        Marsh_point = np.where (np.logical_and (Kernel_DEM >= Kernel_DEM[1,1], Kernel_ridges == 0))
        for j in range(len(Marsh_point[0])):
            X=Marsh_point[0][j]; Y=Marsh_point[1][j]
            Marsh[x+X-1, y+Y-1] = Counter

    print " ... Erase when you've overstepped the line ..."
    Search_marsh_start = np.where (Marsh == 1)
    for i in range(len(Search_marsh_start[0])):
        x=Search_marsh_start[0][i]; y=Search_marsh_start[1][i]
        Kernel_marsh = kernel (Marsh, 3, x, y)
        Kernel_ridges = kernel (Peaks, 3, x, y)
        if np.count_nonzero(Kernel_marsh) <=2:
            Marsh[x,y] = 0


    while Counter < 100:
        Counter = Counter+1
        #print ' ... Filling ... ', Counter
        Search_marsh = np.where (Marsh == Counter-1)
        for i in range(len(Search_marsh[0])):
            x = Search_marsh[0][i]; y = Search_marsh[1][i]
            Kernel_DEM = kernel (DEM, 3, x, y)
            Kernel_DEM_copy = kernel (DEM_copy, 3, x, y)
            Kernel_ridges = kernel (Peaks, 3, x, y)
            Kernel_marsh = kernel (Marsh, 3, x, y)
            Big_Kernel_DEM = kernel (DEM, 11, x, y)
            Big_Kernel_DEM_copy = kernel (DEM_copy, 11, x, y)
            

            Conditions = np.zeros((len(Kernel_DEM), len(Kernel_DEM[0,:])), dtype = np.float)
            # 1: free space
            Condition_1 = np.where (np.logical_and(Kernel_ridges == 0, Kernel_marsh == 0)); Conditions[Condition_1] = 1
            # 2: not topped
            Condition_2 = np.where (np.logical_and(Kernel_DEM_copy > np.amax(Big_Kernel_DEM_copy)-0.2, Conditions == 1)); Conditions[Condition_2] = 2
            
        
            
            #This is a distance thing to make sure you don't cross the ridges agin
            Here_be_ridges = np.where (Kernel_ridges != 0)
            Here_be_parents = np.where (Kernel_marsh == Counter-1)

            for j in range(len(Condition_2[0])):
                X=Condition_2[0][j]; Y=Condition_2[1][j]
                Distance_to_ridges = []
                Distance_to_parents = []

                for k in range(len(Here_be_ridges[0])):
                    Xr=Here_be_ridges[0][k]; Yr=Here_be_ridges[1][k]
                    Distance = np.sqrt((X-Xr)**2+(Y-Yr)**2)
                    Distance_to_ridges.append(Distance)

                for k in range(len(Here_be_parents[0])):
                    Xp=Here_be_parents[0][k]; Yp=Here_be_parents[1][k]
                    Distance = np.sqrt((X-Xp)**2+(Y-Yp)**2)
                    Distance_to_parents.append(Distance)

                if len(Distance_to_ridges)>0:
                    if min(Distance_to_ridges) > min(Distance_to_parents):
                        Marsh[x+X-1, y+Y-1] = Counter
                else:
                    Marsh[x+X-1, y+Y-1] = Counter
                    DEM_copy[x+X-1, y+Y-1] = 0

            """These conditions work! They generate a distribution where you can claerly see if there is some colouring on the tidal flat because you will see a peak at the lowest elevations. All you need to do now is identify that peak and cut it off!"""
            

    
    #This is where you define the cutoff spot!
    Platform = np.copy(Marsh)
    Platform[Platform > 0] = DEM [Platform > 0]
    Platform_bins, Platform_hist = Distribution(Platform,0)

    #1. Find the highest and biggest local maximum of frequency distribution
    
    # Initialize Index
    Index = len(Platform_hist)-1
    # Initiate Cutoff_Z value
    Cutoff_Z = 0
    
    for j in range(1,len(Platform_hist)-1):
        if Platform_hist[j]>0.9*max(Platform_hist) and Platform_hist[j]>Platform_hist[j-1] and Platform_hist[j]>Platform_hist[j+1]:
            Index  = j

    #2. Now run a loop from there toward lower elevations.
    Counter = 0
    for j in range(Index,0,-1):
        # See if you cross the mean value of frequency. Count for how many indices you are under.
        if Platform_hist[j] < np.mean(Platform_hist):
            Counter = Counter + 1
        # Reset the counter value if you go above average again
        else:
            Counter = 0 
    
        #If you stay long enough under (10 is arbitrary for now), initiate cutoff and stop the search
        if Counter > 10:
            Cutoff = j
            Cutoff_Z = Platform_bins[Cutoff]
            break
        
    # If you stay under for more than 5, set a Cutoff_Z value but keep searching    
    if Counter > 5:
        Cutoff = j
        Cutoff_Z = Platform_bins[Cutoff]
        
            
    
    Marsh[Platform<Cutoff_Z] = 0

   
    
    print " ... Fill high gaps ..."
    Search_marsh_condition = np.zeros((len(DEM), len(DEM[0,:])), dtype = np.float)
    Search_marsh = np.where (DEM >= Platform_bins[Index])
    Search_marsh_condition [Search_marsh] = 1
    Search_marsh_2 = np.where (np.logical_and(Marsh == 0, Search_marsh_condition == 1))
    Marsh[Search_marsh_2] = 3



        
    for Iteration in np.arange(0,10,1):
        Counter = 100
        while Counter > 2:
            Counter = Counter-1
            #print " ... Reverse filling ... ", Counter
            Search_marsh = np.where (Marsh == Counter+1)
            Non_filled = 0
            for i in range(len(Search_marsh[0])):
                x = Search_marsh[0][i]; y = Search_marsh[1][i]
                Kernel_DEM = kernel (DEM, 3, x, y)
                Kernel_ridges = kernel (Peaks, 3, x, y)
                Kernel_marsh = kernel (Marsh, 3, x, y)

                if Non_filled <len(Search_marsh[0]):
                    if np.count_nonzero(Kernel_marsh) > 6:
                        Condition = np.where (np.logical_and(Kernel_marsh == 0, Kernel_ridges == 0))
                        for j in range(len(Condition[0])):
                            X=Condition[0][j]; Y=Condition[1][j]
                            Marsh[x+X-1, y+Y-1] = Counter
                    else:
                        Non_filled = Non_filled + 1
                        
    # Reapply the cutoff because the straight line thing is ugly
    Platform = np.copy(Marsh)
    Platform[Platform > 0] = DEM [Platform > 0]
    Marsh[Platform<Cutoff_Z] = 0

 

    # We fill in the wee holes
    print " ... Filling ISOs ... "
    Search_marsh = np.where (np.logical_and(Marsh == 0, Peaks == 0))
    for i in range(len(Search_marsh[0])):
        x = Search_marsh[0][i]; y = Search_marsh[1][i]
        Kernel_marsh = kernel (Marsh, 3, x, y)
        if np.count_nonzero(Kernel_marsh) == 8:
            Marsh[x,y] = 105



    # We get rid of scarps that do not have a marsh next to them
    print " ... Eliminating false scarps ..."
    Search_false_scarp = np.where (Peaks > 0)
    for i in range(len(Search_false_scarp[0])):
        x = Search_false_scarp[0][i]; y = Search_false_scarp[1][i]
        Kernel_marsh = kernel (Marsh, 3, x, y)
        if np.count_nonzero (Kernel_marsh) == 0:
            Peaks[x, y] = 0


    # We get rid of the sticky-outy bits
    print " ... Shave the stubble ..."
    Search_ridge = np.where (Peaks > 0)
    for i in range(len(Search_ridge[0])):
        x=Search_ridge[0][i]; y=Search_ridge[1][i]
        Kernel_ridges = kernel (Peaks, 9, x, y)
        if np.count_nonzero(Kernel_ridges) < 8:
            Peaks[x,y] = 0


    
    
    # We put the scarps in the platform
    print " ... Filling ridges ..."
    Search_side = np.where (Peaks > 0)
    Marsh[Search_side] = 110
    
    
    
    
    
    
    # Some of our platforms are patchy. Try filling them now that we have added the scarps
    
    
    print " ... Fill high gaps ..."
    Search_marsh_condition = np.zeros((len(DEM), len(DEM[0,:])), dtype = np.float)
    Search_marsh = np.where (DEM >= Platform_bins[Index])
    Search_marsh_condition [Search_marsh] = 1
    Search_marsh_2 = np.where (np.logical_and(Marsh == 0, Search_marsh_condition == 1))
    Marsh[Search_marsh_2] = 3
    
    for Iteration in np.arange(0,10,1):
        Counter = 110
        while Counter > 2:
            Counter = Counter-1
            #print " ... Reverse filling ... ", Counter
            Search_marsh = np.where (Marsh == Counter+1)
            Non_filled = 0
            for i in range(len(Search_marsh[0])):
                x = Search_marsh[0][i]; y = Search_marsh[1][i]
                Kernel_DEM = kernel (DEM, 3, x, y)
                Kernel_ridges = kernel (Peaks, 3, x, y)
                Kernel_marsh = kernel (Marsh, 3, x, y)

                if Non_filled <len(Search_marsh[0]):
                    if np.count_nonzero(Kernel_marsh) > 6:
                        Condition = np.where (np.logical_and(Kernel_marsh == 0, Kernel_ridges == 0))
                        for j in range(len(Condition[0])):
                            X=Condition[0][j]; Y=Condition[1][j]
                            Marsh[x+X-1, y+Y-1] = Counter
                    else:
                        Non_filled = Non_filled + 1
                        
    # Reapply the cutoff because the straight line thing is ugly
    Platform = np.copy(Marsh)
    Platform[Platform > 0] = DEM [Platform > 0]
    Marsh[Platform<Cutoff_Z] = 0
    

    
    
    
    Marsh[DEM == Nodata_value] = Nodata_value

    return Marsh

 
    #-------------------------------------------
    # This is the cleaning process

    """print " ... Find the cutoff points ..."
    Marsh_topo = np.copy(Marsh); Marsh_topo[Marsh>0] = DEM[Marsh>0]
    Z_data = Marsh_topo.ravel(); Z_data = Z_data[Z_data>0]
    Z_value = np.arange(min(Z_data), max(Z_data), 0.05)
    hist,bins = np.histogram (Z_data, Z_value, density=True); hist=hist/sum(hist); bins=bins[:-1]

    Highest_platform_freq = hist[np.where(hist == max(hist))]
    Highest_platform_Z = min(bins[np.where(hist == max(hist))])

    High_Cutoff_fraction = 0.08
    Low_Cutoff_fraction = 0.02
    Low_cutoff = []
    High_cutoff = []

    for i in range(1,len(hist)-1):
        if hist[i] <= High_Cutoff_fraction*Highest_platform_freq and bins[i] < Highest_platform_Z:
            Low_cutoff.append(bins[i])
        if hist[i] <= Low_Cutoff_fraction*Highest_platform_freq and bins[i] > Highest_platform_Z:
            High_cutoff.append(bins[i])

    if len(Low_cutoff) > 0:
        Cutoff_Z_min = max(Low_cutoff)
    else:
        Cutoff_Z_min = 0

    if len(High_cutoff) > 0:
        Cutoff_Z_max = min(High_cutoff)
    else:
        Cutoff_Z_max = 0




    print " ... Salting isolated lowlands ..."
    Search_marsh = np.where (np.logical_and(Marsh > 0, DEM < Cutoff_Z_min))
    for i in range(len(Search_marsh[0])):
        x = Search_marsh[0][i]; y = Search_marsh[1][i]
        Kernel_DEM = kernel (DEM, 3, x, y)
        Kernel_ridges = kernel (Peaks, 3, x, y)

        if np.count_nonzero(Kernel_ridges) == 0:
            Marsh[x,y] = 0


            
    print " ... Eliminating false platforms ..."
    Search_false_platform = np.where (np.logical_and(Marsh > 0, Marsh < 2))
    for i in range(len(Search_false_platform[0])):
        x = Search_false_platform[0][i]; y = Search_false_platform[1][i]
        Kernel_marsh = kernel (Marsh, 5, x, y)
        if np.amax(Kernel_marsh) < 3:
            Marsh[x, y] = 0"""


    #fig = plt.figure(3, facecolor='White',figsize=[9,9])
    #ax5 = plt.subplot2grid((1,1),(0,0),colspan=1)
    #ax5.plot(bins, hist)
    #ax5.axvline(Cutoff_Z_min, color='green', lw=2, alpha=0.5)
    #ax5.axvline(Cutoff_Z_max, color='green', lw=2, alpha=0.5)
    #ax5.set_ylim(ymin=0)
    #plt.show()
    #plt.savefig('Output/0_%d_Platform_Distribution.png' % (TR))




    """print " ... Fill high gaps ..."
    Search_marsh_condition = np.zeros((len(DEM), len(DEM[0,:])), dtype = np.float)
    Search_marsh = np.where (DEM >= Highest_platform_Z)
    Search_marsh_condition [Search_marsh] = 1
    Search_marsh_2 = np.where (np.logical_and(Marsh == 0, Search_marsh_condition == 1))
    Marsh[Search_marsh_2] = 3"""



        



    """for A in np.arange(0,10,1):
        Counter = 100
        while Counter > 2:
            Counter = Counter-1
            print " ... Reverse filling ... ", Counter
            Search_marsh = np.where (Marsh == Counter+1)
            Non_filled = 0
            for i in range(len(Search_marsh[0])):
                x = Search_marsh[0][i]; y = Search_marsh[1][i]
                Kernel_DEM = kernel (DEM, 3, x, y)
                Kernel_ridges = kernel (Peaks, 3, x, y)
                Kernel_marsh = kernel (Marsh, 3, x, y)

                if Non_filled <len(Search_marsh[0]):
                    if np.count_nonzero(Kernel_marsh) > 6:
                        Condition = np.where (np.logical_and(Kernel_marsh == 0, Kernel_ridges == 0))
                        for j in range(len(Condition[0])):
                            X=Condition[0][j]; Y=Condition[1][j]
                            Marsh[x+X-1, y+Y-1] = Counter
                    else:
                        Non_filled = Non_filled + 1"""








    """print " ... Find the cutoff points ... again ..."
    Marsh_topo = np.copy(Marsh); Marsh_topo[Marsh>0] = DEM[Marsh>0]
    Z_data = Marsh_topo.ravel(); Z_data = Z_data[Z_data>0]
    Z_value = np.arange(min(Z_data), max(Z_data), 0.05)
    hist,bins = np.histogram (Z_data, Z_value, density=True); hist=hist/sum(hist); bins=bins[:-1]

    Highest_platform_freq = hist[np.where(hist == max(hist))]
    Highest_platform_Z = min(bins[np.where(hist == max(hist))])

    High_Cutoff_fraction = 0.5
    Low_cutoff = []

    for i in range(1,len(hist)-1):
        if hist[i] <= High_Cutoff_fraction*Highest_platform_freq and bins[i] < Highest_platform_Z:
            Low_cutoff.append(bins[i])
        else:
            Low_cutoff.append(0)

    Cutoff_Z_min = max(Low_cutoff)

    fig = plt.figure(3, facecolor='White',figsize=[9,9])
    ax5 = plt.subplot2grid((1,1),(0,0),colspan=1)
    ax5.plot(bins, hist)
    ax5.axvline(Cutoff_Z_min, color='green', lw=2, alpha=0.5)
    ax5.set_ylim(ymin=0)
    plt.show()
    plt.savefig('Output/0_%d_Platform_Cutoff.png' % (TR))"""




    #print " ... Filling ridges ..."
    #Search_side = np.where (Peaks > 0)
    #Marsh[Search_side] = 20

    #for i in range(len(Search_side[0])):
        #x = Search_side[0][i]; y = Search_side[1][i]
        #Kernel_marsh = kernel (Marsh, 7, x, y)

        #if np.sum(Kernel_marsh) > 0 and  np.amin(Kernel_marsh[Kernel_marsh>0]) != Counter + 5:
        #if np.amin(Kernel_marsh) != 105:
            #Marsh[x, y] = 105





    """print " ... Filling ISOs ... "

    Search_marsh = np.where (np.logical_and(Marsh == 0, Peaks == 0))

    for i in range(len(Search_marsh[0])):
        x = Search_marsh[0][i]; y = Search_marsh[1][i]
        Kernel_marsh = kernel (Marsh, 3, x, y)
        if np.count_nonzero(Kernel_marsh) == 8:
            Marsh[x,y] = 200"""









    """Here_be_ridges = np.where (Kernel_ridges != 0)
    Here_be_parents = np.where (Kernel_marsh == Counter-1)


    Distance_to_ridges = []
    Distance_to_parents = []

    for k in range(len(Here_be_ridges[0])):
        Xr=Here_be_ridges[0][k]; Yr=Here_be_ridges[1][k]
        Distance = np.sqrt((X-Xr)**2+(Y-Yr)**2)
        Distance_to_ridges.append(Distance)

    for k in range(len(Here_be_parents[0])):
        Xp=Here_be_parents[0][k]; Yp=Here_be_parents[1][k]
        Distance = np.sqrt((X-Xp)**2+(Y-Yp)**2)
        Distance_to_parents.append(Distance)

    if len(Distance_to_ridges)>0:
        if min(Distance_to_ridges) > min(Distance_to_parents):
            Marsh[x+X-1, y+Y-1] = Counter
    else:
        Marsh[x+X-1, y+Y-1] = Counter
        DEM_copy[x+X-1, y+Y-1] = 0"""



    """Holes = np.zeros((len(DEM), len(DEM[0,:])), dtype = np.float)
    Search_holes = np.where (np.logical_and(Marsh == 0, Peaks == 0))
    Holes[Search_holes] = -1
    Search_holes = np.where (Holes == -1)

    for i in range(len(Search_holes[0])):
        x = Search_holes[0][i]; y = Search_holes[1][i]
        Kernel_holes = kernel (Holes, 21, x, y)


        Distance_to_marshes = []
        Here_be_marshes = np.where(Kernel_holes == 0)


        for j in range(len(Here_be_marshes[0])):
            Xm = Here_be_marshes[0][j]; Ym = Here_be_marshes[1][j]
            Distance = np.sqrt((10-Xm)**2+(10-Ym)**2)
            Distance_to_marshes.append(Distance)

            #print 'D', len(len(Here_be_marshes[0]))
            #print 'D', Distance
            #print 'Dlist', Distance_to_marshes

        if len(Distance_to_marshes)>0:
            Closest_marsh = min(Distance_to_marshes)
            Furthest_marsh= max(Distance_to_marshes)
            Holes[x,y] = Furthest_marsh


            Holes[Holes == 0] = 16"""




    """#Search_marsh = np.where (np.logical_and(Marsh > 2, Peaks == 0, DEM > Cutoff_Z_min))
    for i in range(len(Search_marsh[0])):
        x = Search_marsh[0][i]; y = Search_marsh[1][i]
        Kernel_marsh = kernel (Marsh, 9, x, y)
        Kernel_ridges = kernel (Peaks, 9, x, y)
        Kernel_DEM = kernel (DEM, 9, x, y)

        Here_be_holes = np.where(Kernel_marsh == 0)
        Here_be_ridges = np.where (Kernel_ridges > 0)

        Distance_hole_to_ridges = np.zeros((len(Kernel_marsh), len(Kernel_marsh[0,:])), dtype = np.float)


        for k in range(len(Here_be_holes[0])):
            Xh=Here_be_holes[0][k]; Yh=Here_be_holes[1][k]
            for l in range(len(Here_be_ridges[0])):
                Xr=Here_be_ridges[0][l]; Yr=Here_be_ridges[1][l]

                Distance = np.sqrt((Xh-Xr)**2+(Yh-Yr)**2)
                Distance_hole_to_ridges.append(Distance)

            Closest_ridge = np.where()


        if len(Distance_to_ridges)>0:
            if min(Distance_to_ridges) > min(Distance_to_holes):
                Marsh[x+X-1, y+Y-1] = 3"""





    """print " ... Eliminating false platforms ..."
    Search_false_platform = np.where (Marsh <= 3)
    for i in range(len(Search_false_platform[0])):
        x = Search_false_platform[0][i]; y = Search_false_platform[1][i]
        Kernel_marsh = kernel (Marsh, 5, x, y)
        Kernel_ridges = kernel (Peaks, 5, x, y)
        if np.count_nonzero (Kernel_marsh) < 10: #np.amax (Kernel_ridges) < 3 and
            Marsh[x, y] = 0"""



    """print " ... Find the cutoff points ..."
    Marsh_topo = np.copy(Marsh); Marsh_topo[Marsh>0] = DEM[Marsh>0]
    Z_data = Marsh_topo.ravel(); Z_data = Z_data[Z_data>0]
    Z_value = np.arange(min(Z_data), max(Z_data), 0.05)
    hist,bins = np.histogram (Z_data, Z_value, density=True); hist=hist/sum(hist); bins=bins[:-1]

    Highest_platform_freq = hist[np.where(hist == max(hist))]
    Highest_platform_Z = min(bins[np.where(hist == max(hist))])

    High_Cutoff_fraction = 0.05
    Low_Cutoff_fraction = 0.02
    Low_cutoff = []
    High_cutoff = []

    for i in range(1,len(hist)-1):
        if hist[i] <= High_Cutoff_fraction*Highest_platform_freq and bins[i] < Highest_platform_Z:
            Low_cutoff.append(bins[i])
        if hist[i] <= Low_Cutoff_fraction*Highest_platform_freq and bins[i] > Highest_platform_Z:
            High_cutoff.append(bins[i])

    Cutoff_Z_min = max(Low_cutoff)
    Cutoff_Z_max = min(High_cutoff)


    fig = plt.figure(3, facecolor='White',figsize=[9,9])
    ax5 = plt.subplot2grid((1,1),(0,0),colspan=1)
    ax5.plot(bins, hist)
    ax5.axvline(Cutoff_Z_min, color='green', lw=2, alpha=0.5)
    ax5.axvline(Cutoff_Z_max, color='green', lw=2, alpha=0.5)
    ax5.set_ylim(ymin=0)
    plt.show()
    #plt.savefig('Output/Poster/0_%d_Platform_Distribution.png' % (TR))"""





    """print " ... Shaving hairy ridges ..."
    Search_marsh = np.where (Marsh > Counter)
    #Search_marsh = np.where (np.logical_and(Marsh > 0, DEM < Cutoff_Z_min))
    for i in range(len(Search_marsh[0])):
        x = Search_marsh[0][i]; y = Search_marsh[1][i]
        #Kernel_DEM = kernel (DEM, 9, x, y)
        #Kernel_ridges = kernel (Peaks, 9, x, y)
        Kernel_marsh = kernel (Marsh, 3, x, y)

        if np.amax(Kernel_marsh[Kernel_marsh < Counter]) == 0:
            Marsh[x,y] = 0"""




            #Here_be_sons = np.where (Kernel_marsh > 2)
            #Here_be_parents = np.where (Kernel_marsh == 1)
            #Distance_to_sons = []
            #Distance_to_parents = []

            #for k in range(len(Distance_to_sons[0])):
                #Xs=Distance_to_sons[0][k]; Ys=Distance_to_sons[1][k]
                #Distance = np.sqrt((1-Xr)**2+(1-Yr)**2)
                #Distance_to_sons.append(Distance)

            #for k in range(len(Here_be_parents[0])):
                #Xp=Here_be_parents[0][k]; Yp=Here_be_parents[1][k]
                #Distance = np.sqrt((1-Xp)**2+(1-Yp)**2)
                #Distance_to_parents.append(Distance)

            #if len(Distance_to_sons)>0 and len(Distance_to_parents)>0:
                #if min(Distance_to_sons) > min(Distance_to_sons):
                    #Marsh[x+X-1, y+Y-1] = Counter + 10



"""Search_holes = np.where (np.logical_and(Marsh == 0, Peaks == 0, DEM>=Cutoff_Z_min))
for i in range(len(Search_holes[0])):
    x = Search_holes[0][i]; y = Search_holes[1][i]
    Kernel_marsh = kernel (Marsh, 3, x, y)
    Kernel_ridges = kernel (Peaks, 3, x, y)

    Here_be_ridges = np.where (Kernel_ridges > 0)
    Here_be_parents = np.where (Kernel_marsh > 0)
    Distance_to_ridges = []
    Distance_to_parents = []

    for k in range(len(Here_be_ridges[0])):
        Xr=Here_be_ridges[0][k]; Yr=Here_be_ridges[1][k]
        Distance = np.sqrt((1-Xr)**2+(1-Yr)**2)
        Distance_to_ridges.append(Distance)

    for k in range(len(Here_be_parents[0])):
        Xp=Here_be_parents[0][k]; Yp=Here_be_parents[1][k]
        Distance = np.sqrt((1-Xp)**2+(1-Yp)**2)
        Distance_to_parents.append(Distance)

    if len(Distance_to_ridges)>0 and len(Distance_to_parents)>0:
        if min(Distance_to_ridges) > min(Distance_to_parents):
            Marsh[x+X-1, y+Y-1] = Counter + 10"""


# Then fill holes
# Then get rid of ridge/counter1 combos
# Then clean the ridgelets
# Then get rid of low platform that are lower than the closest ridge


    #Marsh[DEM<Cutoff_Z_min] = 0
    #Marsh[DEM>Cutoff_Z_max] = 0




"""print " Don't colour over the lines"
    while Counter < 100:
        Counter = Counter+1
        print "Stalling ... ", Counter
        Search_marsh = np.where (Marsh == Counter-1)
        for i in range(len(Search_marsh[0])):
            x = Search_marsh[0][i]; y = Search_marsh[1][i]
            Kernel_DEM = kernel (DEM, 3, x, y)
            Kernel_relief = Kernel_DEM-np.amin(Kernel_DEM)
            Kernel_ridges = kernel (Peaks, 3, x, y)
            Kernel_marsh = kernel (Marsh, 3, x, y)

            Conditions = np.zeros((len(Kernel_DEM), len(Kernel_DEM[0,:])), dtype = np.float)
            # 1: free, flat space
            #Condition_1 = np.where (np.logical_and(Kernel_ridges == 0, Kernel_marsh == 0))
            Condition_1 = np.where (np.logical_and(Kernel_ridges == 0, Kernel_marsh == 0, Kernel_relief < 0.2))
            Conditions[Condition_1] = 1
            # 2: Inside the limits
            Condition_2 = np.where (np.logical_and(Conditions == 1, Kernel_DEM > 0.8*Cutoff_Z_min, Kernel_DEM < 1.2*Cutoff_Z_max))
            Conditions[Condition_2] = 2

            for j in range(len(Condition_1[0])):
                X=Condition_1[0][j]; Y=Condition_1[1][j]
                Marsh[x+X-1, y+Y-1] = Counter


                print Marsh[x+X-1, y+Y-1]

    print np.amax(Marsh)"""




"""print " Eliminating false platforms"
    Search_false = np.where (Marsh == 1)
    for i in range(len(Search_false[0])):
        x = Search_false[0][i]; y = Search_false[1][i]
        Kernel_marsh = kernel (Marsh, 3, x, y)
        Kernel_marsh2 = kernel (Marsh, 9, x, y)
        if np.amax(Kernel_marsh) == 1 and np.amax(Kernel_marsh2) <= 2:
            False_point = np.where (Kernel_marsh == 1)
            for j in range(len(False_point[0])):
                X=False_point[0][j]; Y=False_point[1][j]
                if x+X-1 < len(DEM) and y+Y-1 < len(DEM[0,:]):
                    Marsh[x+X-1, y+Y-1] = 0


    print " Filling internal ridges"
    Search_interior = np.where (Peaks > 0)
    for i in range(len(Search_interior[0])):
        x = Search_interior[0][i]; y = Search_interior[1][i]
        Kernel_peaks = kernel (Peaks, 3, x, y)
        Kernel_marsh = kernel (Marsh, 3, x, y)
        Kernel_sum = abs(Kernel_peaks)+abs(Kernel_marsh)

        if np.amin(Kernel_sum) != 0 :
            Marsh[x, y] = 100"""







"""print " Filling holes"
    for number in np.arange(0,9,1):
        print number

        Search_holes = np.where (np.logical_and(Marsh == 0, Peaks == 0, DEM>=Cutoff_Z_min))
        for i in range(len(Search_holes[0])):
            x = Search_holes[0][i]; y = Search_holes[1][i]
            Kernel_marsh = kernel (Marsh, 3, x, y)

            Here_be_ridges = np.where (Kernel_marsh == 100)
            Here_be_parents = np.where (Kernel_marsh > 0)
            Distance_to_ridges = []
            Distance_to_parents = []

            for k in range(len(Here_be_ridges[0])):
                Xr=Here_be_ridges[0][k]; Yr=Here_be_ridges[1][k]
                Distance = np.sqrt((1-Xr)**2+(1-Yr)**2)
                Distance_to_ridges.append(Distance)

            for k in range(len(Here_be_parents[0])):
                Xp=Here_be_parents[0][k]; Yp=Here_be_parents[1][k]
                Distance = np.sqrt((1-Xp)**2+(1-Yp)**2)
                Distance_to_parents.append(Distance)

            if len(Distance_to_ridges)>0 and len(Distance_to_parents)>0:
                if min(Distance_to_ridges) > min(Distance_to_parents):
                    Marsh[x+X-1, y+Y-1] = 200

    for number in np.arange(0,9,1):
        print number
        Search_holes = np.where (np.logical_and(Marsh == 0, Peaks == 0, DEM>=Cutoff_Z_min))
        for i in range(len(Search_holes[0])):
            x = Search_holes[0][i]; y = Search_holes[1][i]
            Kernel_marsh = kernel (Marsh, 3, x, y)

        if np.count_nonzero(Kernel_marsh)>6:
            Marsh[x,y] = 300"""



"""print " Sinking islands"
    Search_islands = np.where (Marsh != 0)
    for i in range(len(Search_islands[0])):
            x = Search_islands[0][i]; y = Search_islands[1][i]
            Kernel_marsh = kernel (Marsh, 3, x, y)
            if np.amin(Kernel_marsh) == 0:
                if np.count_nonzero(Kernel_marsh) < 3:
                    Marsh[x, y] = 0"""




"""print " Eliminating false scarps"
    Search_false_scarp = np.where (Peaks > 0)
    for i in range(len(Search_false_scarp[0])):
        x = Search_false_scarp[0][i]; y = Search_false_scarp[1][i]
        Kernel_marsh = kernel (Marsh, 3, x, y)
        if np.amax (Kernel_marsh) == 0:
            Peaks[x, y] = 0"""




"""# But also we WHAT???
for i in range(len(Search_side[0])):
    x = Search_side[0][i]; y = Search_side[1][i]
    Kernel_marsh = kernel (Marsh, 7, x, y)

    if np.sum(Kernel_marsh) > 0 and  np.amin(Kernel_marsh[Kernel_marsh>0]) != Counter + 5:
        if np.amin(Kernel_marsh) != 105:
            Marsh[x, y] = 110"""




    



#---------------------------------------------------------------
# This function is the MarshFinder: it finds marsh platforms, scarps and pioneer zones
# Here's how the MarshFinder works:
# STEP 0: Identify a sesarch space where elevation is above the Neap Low Tide and slope is in the higher 50%
# STEP 1: Identify scarps by:
#                          identifying local slope maxima (peaks)
#                          starting ridge lines along the highest slope values stemming from these peaks
#                          continuing the ridges along the shortest path of high slopes
# STEP 2: Clean the ridges if:
#                           they are too short
#                           the local relief is too small
#                           they are lower than the most frequent elevation of ridges minus a constant depending on tidal range
# STEP 3: Fill ground above ridges by:
#                                   filling ground directly in contact with ridges and above the ridge
#
# STEP 4:
# STEP 5:
# STEP 6:
# STEP 7:
# STEP 8:
# It takes as input:
# 1/ the DEM
# 2/ the Slopes
# 3/ the Curvature
# 4/ the Channels
# 5/ the Tidalstatistix
# It returns:
# 1/ the marsh platforms
# 2/ the marsh scarps
# 3/ the marsh channels
# 4/ the pioneer zones



#def MARSH_ID (DEM, Slope, Curvature, Metric2, Nodata_value):
#Added this line for optimisation purposes
def MARSH_ID (DEM, Slope, Curvature, Nodata_value, opt1, opt2, opt3):
    DEM_work = np.copy(DEM); Slope_work = np.copy(Slope); Curvature_work = np.copy(Curvature) #; Channels_work = np.copy(Channels)

    Platform = np.copy(DEM_work)
    Ridge = np.copy(DEM_work)
    Marsh = np.copy(DEM_work)

    Platform[Platform != Nodata_value] = 0
    Summit = np.where (Platform==np.amax(Platform))
    Platform[Summit] = 1

    #---------------------------------------
    TR_spring = np.mean (Metric2[3])-np.mean (Metric2[0])

    Search_space, Crossover, bins, hist, Inflexion_point = define_search_space (DEM_work, Slope_work, Nodata_value,opt1)

    Order = 1
    Ridge, Slope_temp = peak_flag (Slope_work, Search_space, Order)

    Order = Order+1
    Ridge, Slope_temp = initiate_ridge (Slope_temp, Search_space, Ridge, Order)

    while Order < 50:
        Order = Order+1
        Ridge, Slope_temp = Continue_ridge (DEM, Slope_temp, Search_space, Ridge, Order, Metric2)

    Ridge = Clean_ridges (Ridge, DEM_work, Slope_work, Metric2, Nodata_value, opt2)

    Marsh = Fill_high_ground (DEM_work, Ridge, Metric2, Nodata_value, opt3)


    print "My hovercraft is full of eels!"
    print


    return Search_space, Ridge, Marsh




#-----------------------------------------------------------------------------------------------------
# This functions compares the marsh identified automatically to a reference marsh, usually digitised from a picture
# It then builds a confusion matrix to see if the MarshFinder has performed as expected.

# It takes as input:
# 1/ The subject marsh array
# 2/ The reference marsh array


# It returns:
# 1/ An array of the confusion matrix and its associated metrix (in a raster and a figure)

def Confusion (Subject, Reference, Nodata_value):
    Height = len(Subject[:,0]); Width = len(Subject[0,:])
    Height_R = len(Reference[:,0]); Width_R = len(Reference[0,:])
    
    
    
    print Height, Width
    print Height_R, Width_R
    
    
    H = min (Height, Height_R)
    W = min (Width, Width_R)
    
    

    Confusion_matrix = Nodata_value*np.ones((Height, Width), dtype = np.float)


    Subject_marsh = np.where (np.logical_and(Subject != 0, Subject != Nodata_value))
    Reference_marsh = np.where (np.logical_and(Reference != 0, Reference != Nodata_value))

    Subject[Subject_marsh] = 1.
    Reference[Reference_marsh] = 1.

    for i in range (H):
        for j in range (W):
            if Subject[i,j] == 1 and Reference[i,j] == 1: # TRUE POSITIVE
                Confusion_matrix[i,j] = 1
            elif Subject[i,j] == 0 and Reference[i,j] == 0: # TRUE NEGATIVE
                Confusion_matrix[i,j] = 2
            elif Subject[i,j] == 1 and Reference[i,j] == 0: # FALSE POSITIVE
                Confusion_matrix[i,j] = -1
            elif Subject[i,j] == 0 and Reference[i,j] == 1: # FALSE NEGATIVE
                Confusion_matrix[i,j] = -2

    True_positive = np.sum(Confusion_matrix[Confusion_matrix == 1])
    True_negative = np.sum(Confusion_matrix[Confusion_matrix == 2])/2
    False_positive = -np.sum(Confusion_matrix[Confusion_matrix == -1])
    False_negative = -np.sum(Confusion_matrix[Confusion_matrix == -2])/2

    Reliability = True_positive / (True_positive+False_positive)
    Sensitivity = True_positive / (True_positive+False_negative)
    Accuracy = (True_positive+True_negative) / (True_positive+True_negative+False_positive+False_negative)
    F1 = 2*True_positive/(2*True_positive+False_positive+False_negative)

    Performance = np.array([True_positive,True_negative,False_positive,False_negative])
    Metrix = np.array([Accuracy, Reliability, Sensitivity, F1])



    return Confusion_matrix, Performance, Metrix





"""x = np.arange(Width)
    y = np.arange(Height)

    fig = plt.figure(1, facecolor='White',figsize=[9,9])
    ax5 = plt.subplot2grid((2,2),(1,0),colspan=1)
    ax6 = plt.subplot2grid((2,2),(0,1),colspan=1)
    ax7 = plt.subplot2grid((2,2),(1,1),colspan=1)
    ax8 = plt.subplot2grid((2,2),(0,0),colspan=1)

    p5 = ax5.imshow(Ridge, interpolation='none')
    #p6 = ax6.imshow(Slope_work, interpolation='none')
    p6 = ax6.imshow(Search_space*Slope_work, interpolation='none')
    p7 = ax7.imshow(Marsh, interpolation='none')
    p8 = ax8.imshow(DEM_work, interpolation='none')

    ax5.set_title("Scarps (order)")
    ax6.set_title("Search space (slopy x high)")
    ax7.set_title("Detected Marsh platform")
    ax8.set_title("Regional DEM")

    fig.colorbar(p5, ax=ax5)
    fig.colorbar(p6, ax=ax6)
    fig.colorbar(p7, ax=ax7)
    fig.colorbar(p8, ax=ax8)

    plt.savefig('Find_marsh.png')
    plt.show()"""







"""#--------------------------------------------------------------
# This is a tryout function that calculates local relief at every point

def Relief (DEM):

    print "TESTING RELIEF-BASED SHITE"

    Height = len(DEM)
    Width = len(DEM[0,:])

    DEM_copy = np.copy(DEM)
    local_relief = np.zeros((Height, Width), dtype=np.float)


    for i in range(Height):
        for j in range(Width):
            Kernel_DEM = kernel (DEM_copy, 21, i, j)
            Kernel_relief = Kernel_DEM-np.amin(Kernel_DEM)

            if np.amax(Kernel_relief) < 0.2:
                local_relief[i,j] = np.amax(Kernel_relief)


    return local_relief"""












"""# STEP 4: Search for connexions
print 'Connecting ridges'

Order = Order+1

#Search in the initial search space
for i in range(len(Search[0])):
    x=Search[0][i]; y=Search[1][i]
    Kernel_slope = Slope_work[x-1:x+2, y-1:y+2] # The kernel in the slope array
    Kernel_slope_temp = Slope_temp[x-1:x+2, y-1:y+2] # The kernel in the modified slope array
    Kernel_ridge = Ridge[x-1:x+2, y-1:y+2] # The kernel in the ridge array
    Kernel_search = Search_space[x-1:x+2, y-1:y+2] # The kernel in the search space
    Kernel_search_mob = Search_space_mob[x-1:x+2, y-1:y+2] # The kernel in the mobile search space

    # Count the number of nonzero points in the kernel of the ridge array
    Ridge_points = np.count_nonzero(Kernel_ridge)

    # But not within the summits already identified and only those with 2 neighbours
    if Kernel_ridge[1,1] == 0 and Ridge_points == 2:
        Ridge[x,y]=Order
        Slope_temp[x,y] = 0"""




"""#---------------------------------------------------------------
# Fictional DEM and Slope to study. This is a simplified version of a salt marsh
DEM_work = np.zeros((30,30), dtype=np.float)
Height = len(DEM_work)
Width = len(DEM_work[0,:])

Slope_work = np.zeros((Height,Width), dtype=np.float)


# Reference elevation for the high marsh, tidal flat, scarp, pioneer zone and low marsh
Zma = 1.2; Ztf = 0.4; Zpz = Ztf + 0.3 * (Zma-Ztf); Zloma = Zma - 0.2 * (Zma-Ztf); Zsc = (Zma+Ztf)/2
Zsc_pz = (Zma+Zpz)/2
Zsc_loma = (Zloma+Ztf)/2

#Slopes of the marsh, tidal flat and channel
Sma = 0.004; Stf = 0.005; Sch = 0.05

# Rugosity of the marsh, tidal flat, scarp and pioneer zone
Rma = 0; Rtf = 0; Rsc = 0; Rpz = 0


for i in range(len(DEM_work)): # rows
    for j in range(len(DEM_work[0,:])): # cols

        # Rugosity of the marsh, tidal flat, scarp and pioneer zone
        #Rma = np.random.rand()/4; Rtf = np.random.rand()/5; Rsc = np.random.rand()/4; Rpz = np.random.rand()/4


        if i==j+1: # the high scarp.
            DEM_work[i,j] = 1.1*Zsc + Rsc
            if j>22: # transition to low marsh
                DEM_work[i,j] = 1.1*Zsc_loma + Rsc
            if j<8: # transition to pioneer zone
                DEM_work[i,j] = 1.1*Zsc_pz + Rsc

        elif i==j: # the mid scarp.
            DEM_work[i,j] = Zsc + Rsc
            if j>22: # transition to low marsh
                DEM_work[i,j] = Zsc_loma + Rsc
            if j<8: # transition to pioneer zone
                DEM_work[i,j] = Zsc_pz + Rsc


        elif i==j-1: # the low scarp.
            DEM_work[i,j] = 0.9*Zsc + Rsc
            if j>22: # transition to low marsh
                DEM_work[i,j] = 0.9*Zsc_loma + Rsc
            if j<8: # transition to pioneer zone
                DEM_work[i,j] = 0.9*Zsc_pz + Rsc

        elif i<j: # the tidal flat. It dips NE
            DEM_work[i,j] = Ztf - Stf * (Height-i+j) + Rtf

            if j<8: # a pioneer zone
                DEM_work[i,j] = Zpz - Stf * (Height-i+j) + Rpz

        elif i>j: # the high marsh. It slightly dips NE
            DEM_work[i,j] = Zma - Sma * (Height-i+j) + Rma

            if j>22: # the low marsh
                DEM_work[i,j] = Zloma - Sma * (Height-i+j) + Rma


# Now make a channel
for i in range(len(DEM_work)): # rows
    for j in range(len(DEM_work[0,:])): # cols
        if (i==Height-j or i==Height-j+1 or i==Height-j-1 or i==Height-j+2 or i==Height-j-2) and i<=26:
            DEM_work[i,j] = DEM_work[26,6] - Sch * (Height-i+j)**0.85


# Now make a slope DEM: http://www.onlinegeographer.com/slope/Dunn_hickey_slope.pdf, equation 7 with the maximum slope method for simplicity
for i in range(1,len(DEM_work)-1): # rows
    for j in range(1,len(DEM_work[0,:])-1):  # cols
        # Here are the cells we consider
        Cells = np.zeros(9, dtype=np.float)
        Slopes = np.zeros(9, dtype=np.float)
        Cells[0]=DEM_work[i,j]; Cells[1]=DEM_work[i-1,j]; Cells[2]=DEM_work[i-1,j+1]; Cells[3]=DEM_work[i,j+1]; Cells[4]=DEM_work[i+1,j+1]; Cells[5]=DEM_work[i+1,j]; Cells[6]=DEM_work[i+1,j-1]; Cells[7]=DEM_work[i,j-1]; Cells[8]=DEM_work[i-1,j-1]

        # Calculate Slopes
        for a in range(1,len(Cells)):
            Slopes[a]= (Cells[a]-Cells[0])/1

        # Largest Slope is the slope value for i,j
        Slope_work[i,j]=max(Slopes)"""


"""print " elevation*slope cleaning"
fig=plt.figure(2, facecolor='White',figsize=[10,10])
ax = plt.subplot2grid((2,1),(0,0),colspan=1, rowspan=1)

Ridge_ZxS = np.copy(Peaks)
Ridge_ZxS[Peaks!=0] = DEM[Peaks!=0]*Slope[Peaks!=0]
ZxS_data = Ridge_ZxS.ravel(); ZxS_data = ZxS_data[ZxS_data!=0]
ZxS_value = np.arange(min(ZxS_data), max(ZxS_data),(max(ZxS_data)-min(ZxS_data))/20)
hist,bins = np.histogram (ZxS_data, ZxS_value, density=True); hist=hist/sum(hist); bins=bins[:-1]
ax.plot(bins,hist)

ax = plt.subplot2grid((2,1),(1,0),colspan=1, rowspan=1)
Ridge_ZxS = np.copy(Peaks)
Ridge_ZxS[Peaks!=0] = DEM[Peaks!=0]*Slope[Peaks!=0] / TR
ZxS_data = Ridge_ZxS.ravel(); ZxS_data = ZxS_data[ZxS_data!=0]
ZxS_value = np.arange(min(ZxS_data), max(ZxS_data), (max(ZxS_data)-min(ZxS_data))/20)
hist,bins = np.histogram (ZxS_data, ZxS_value, density=True); hist=hist/sum(hist); bins=bins[:-1]
ax.plot(bins,hist)

plt.savefig('Output/%d_Distribution.png' % (TR))

# This is the most frequent ridge elevation
#High_ridge_proba = np.where(hist == max(hist)); High_ridge_Z = bins[High_ridge_proba]

# If ridges are well below this elevation (empirical cut-off based on TR), it's unlikely they are marsh ridges

Tidal_factor = np.sqrt(TR)/1000
print Tidal_factor
ZxS_condition = Tidal_factor*np.amax(Ridge_ZxS)
#ZxS_condition = np.percentile(Ridge_ZxS[Ridge_ZxS!=0], (TR-1)*4)
low_ridges = np.where (Ridge_ZxS < ZxS_condition)

Peaks[low_ridges] = 0"""







#---------------------------------------------------------------
def ENVI_raster_binary_to_2d_array(file_name, gauge):
    print 'Opening %s' % (gauge)
    #Source : http://chris35wills.github.io/python-gdal-raster-io/
    '''
    Converts a binary file of ENVI type to a numpy array.
    Lack of an ENVI .hdr file will cause this to crash.
    '''

    driver = gdal.GetDriverByName('ENVI')

    driver.Register()

    inDs = gdal.Open(file_name, GA_ReadOnly)

    if inDs is None:
        print "Couldn't open this file: " + file_name
        print "Perhaps you need an ENVI .hdr file? "
        sys.exit("Try again!")
    else:
        print "%s opened successfully" %file_name

        #print '~~~~~~~~~~~~~~'
        #print 'Get image size'
        #print '~~~~~~~~~~~~~~'
        cols = inDs.RasterXSize
        rows = inDs.RasterYSize
        bands = inDs.RasterCount

        #print "columns: %i" %cols
        #print "rows: %i" %rows
        #print "bands: %i" %bands

        #print '~~~~~~~~~~~~~~'
        #print 'Get georeference information'
        #print '~~~~~~~~~~~~~~'
        geotransform = inDs.GetGeoTransform()
        originX = geotransform[0]
        originY = geotransform[3]
        pixelWidth = geotransform[1]
        pixelHeight = geotransform[5]

        #print "origin x: %i" %originX
        #print "origin y: %i" %originY
        #print "width: %2.2f" %pixelWidth
        #print "height: %2.2f" %pixelHeight

        # Set pixel offset.....
        print '~~~~~~~~~~~~~~'
        print 'Convert image to 2D array'
        print '~~~~~~~~~~~~~~'
        band = inDs.GetRasterBand(1)
        image_array = band.ReadAsArray(0, 0, cols, rows)
        image_array_name = file_name
        print type(image_array)
        print image_array.shape

        return image_array, pixelWidth, (geotransform, inDs)




#---------------------------------------------------------------------
def ENVI_raster_binary_from_2d_array(envidata, file_out, post, image_array):
    #Source : http://chris35wills.github.io/python-gdal-raster-io/
    #util.check_output_dir(file_out)

    driver = gdal.GetDriverByName('ENVI')

    original_geotransform, inDs = envidata

    rows, cols = image_array.shape
    bands = 1

    # Creates a new raster data source
    outDs = driver.Create(file_out, cols, rows, bands, gdal.GDT_Float32)

    # Write metadata
    originX = original_geotransform[0]
    originY = original_geotransform[3]

    outDs.SetGeoTransform([originX, post, 0.0, originY, 0.0, -post])
    outDs.SetProjection(inDs.GetProjection())

    #Write raster datasets
    outBand = outDs.GetRasterBand(1)
    outBand.WriteArray(image_array)

    new_geotransform = outDs.GetGeoTransform()
    new_projection = outDs.GetProjection()

    print "Output binary saved: ", file_out

    return new_geotransform,new_projection,file_out


#-------------------------------------------------------------------
#7.
# Some subaxes
#http://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib
def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax
