# LSD_TopoTools_MarshExtraction #
This repository contains several Python scripts you can use to detect saltmarsh platforms and outlines on a DEM.
In this README you will find:

1) Instructions to use this programme successfully.

2) A description of the input data files provided.

3) Some comments.


## User Manual ##



### Example_Input_prep.py ###
This file prepares your DEM (or DEMs) for analysis.
Instructions and process breakdown are included in the script.
Please note that this script requires:
1. GDAL
2. a Wiener filter, also used in LSDTopoTools for channel extractions. The makefile "Wiener_filter.make" can be found at: https://github.com/LSDtopotools/LSDTopoTools_ChannelExtraction/tree/master/driver_functions_ChannelExtraction
3. a basic LSD TopoTools topographic analysis package and the associated driver file. You can find the analysis package here: https://github.com/LSDtopotools/LSDTopoTools_AnalysisDriver. A template for the driver file is located in the Input/ folder under the name "FEL_slope_curv.LSDTT_driver".

We recommend that you run this file independently from the rest of the analysis as it may need to be run from several locations, depending on your file organisation.

#### Do I really need this file? #### 
If your DEM is already associated to a slope raster, you may run the marsh identification without using this file. However, we do not guarantee results if the slope was calculated via a method different from that used in our preprocessing file.



### Example_Motherscript.py ###
This file runs the identification script "Example_Marsh_ID.py" and the plotting script "Example_Plots.py".
You can also set it to run the preparation script "Example_Input_prep.py" if you have already computed the slope, curvature and hillshade rasters. This is particularly useful if you want to change the domain of your analysis or modify your reference.


### Example_Marsh_ID.py ###
This is the file which runs the identification process. In this file you need to specify the site you are working with and the parameter values. In exchange, you will get three .bil rasters: one for the initial search space used to narrow down the search for scarps, one for the scarps, and one for the marsh platforms.
Note that the scarps file does not necessarily correspond to the edges of the marsh platform. It is merely an intermediate step.

If you have a reference marsh outline, you will also get a confusion map as a .bil raster to evaluate the quality of the detection process.



### Example_Plots.py ###
This file produces plots so that you may visualise your results. Please read the description of each plot in the file to know what you are looking at.


### Example_functions.py ###
This file contains all the useful functions used to identify marsh platforms and plot your results. Please do not touch this file unless you know what you are doing.




## Use and Abuse ##
Preferably use this script on saltmarshes with distinct scarp features separating them from the tidal flat.
Sparsely vegetated pioneer zones that have no impact on topography will lead to a strong decrease in accuracy.
