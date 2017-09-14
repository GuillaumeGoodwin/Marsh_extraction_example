# LSD_TopoTools_MarshExtraction #
This repository contains Python scripts that detect saltmarsh platforms and outlines on a DEM.
In this README you will find indications to use this programme and a few comments on the use and abuse of it.


## How to use ##


### Motherscript.py ###
This file runs the entire analysis.
If the other files are correctly set to read your input, it will produce .bil files of saltmarsh platforms and outlines as well as a variety of plots.
However, unless you use this with the original dataset, these plots will be difficult to interpret or may be misleading.
Please edit the plotting file if you want plots that mean something for you.


### Prepare_Input_DEM.py ###
This file prepares your DEM (or DEMs) for analysis.
Instructions and process breakdown are included in the script.
Please note that this script requires:
1. GDAL
2. a Wiener filter, also used in LSDTopoTools for channel extractions. The makefile "Wiener_filter.make" can be found at: https://github.com/LSDtopotools/LSDTopoTools_ChannelExtraction/tree/master/driver_functions_ChannelExtraction
3. a basic LSD TopoTools topographic analysis package and the associated driver file. You can find the analysis package here: https://github.com/LSDtopotools/LSDTopoTools_AnalysisDriver. A template for the driver file is located in this repository under the name "Example_driver.LSDTT_driver".




## Use and Abuse ##
Preferably use this script on saltmarshes with distinct scarp features separating them from the tidal flat.
Sparsely vegetated pioneer zones that have no impact on topography will lead to a strong decrease in accuracy.
