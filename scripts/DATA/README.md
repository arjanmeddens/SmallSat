# DATA

This sub-directory contains the relevant files that are `INPUT` and `OUTPUT` for the classificaion pipeline. The sub-folders in this directory include:

* `cloud_mask`: contains the cloud mask raster used in the classification pipeline.
* `highres_images`: contains the very high resolution (VHR) imagery used to train and predict the classification models on.
* `output`: contains any relevant `OUTPUT` data from the classification pipeline. The `OUTPUT` data include: classification raster, confusion matrices, metrics, and run parameters.
* `shapefiles`: contains the vector file (.shp) used to define the extent of the processing for the classification workflow. This folder also holds the .csv file containing the run parameters for the classification workflow.
