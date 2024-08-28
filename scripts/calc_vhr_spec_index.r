################################################################################################
# Purpose: Function to calculate Spectral indices from WV, QB, or GeoEye imagery for 
#   detecting insect disturbance using very high-res satellite data
################################################################################################
# Author: Arjan Meddens, Washington State University (arjan.meddens@wsu.edu)
# Funding: NASA SmallSat grant: 80NSSC21K1155
# Date: July 25, 2024
# Publication: Meddens et al. In review: Detecting multiple insect disturbance types across the 
#   western US using very high-resolution satellite data, Remote Sensing of Environment.  
#- - - - - - - - - - - - - - - 
# R version: 4.2.2 (2022-10-31 ucrt) -- "Innocent and Trusting" 
#- - - - - - - - - - - - - - - 
################################################################################################
# Input format:  .tif (VHR satelite spectral data) 
# Output format: .tif (index + spectral data cube) 
################################################################################################

calc_wv_index = function(wv_data,num_bands=num_bands) {
#' Calculate Indices for Very-High Resolution Imagery 
#'
#' @param wv_data `RasterStack`; input imagery
#' @param num_bands `numeric` (4 or 8); number of bands of the input imagery. Currently supports calulcating indices for `num_bands = 4` and `num_bands = 8`.
#'
#' @return `RasterStack`; Bands returned (in order) for:
#'   * 4-band input: `c('BLU','GRE','RED','NR1','NDVI1','EVI','RGI')`
#'   * 8-band input: `c('OCE','BLU','GRE','YEL','RED','EDG','NR1','NR2','NDVI1','NDVI2','NDRE1','NDRE2','EVI','RGI')`
#' @export
#'
#' @examples
  if (class(wv_data) == "RasterBrick") {
    #print(“wv_data correct format”)
  } else {
    message(paste('wv_data should be a RasterBrick'))
  }

  if (length(num_bands) == 0) {
    num_bands == dims(wv_data)[3]
    message(paste('num_bands taken from RasterBrick...'))
  } else {
    message(paste('num_bands is ',num_bands))
  }

  #--- Define constants ---#
  #-- See Chen et al 2004, RSE IKONOS LAI study --#
  c_g  = c(2.5)
  c_c1 = c(6.0)
  c_c2 = c(7.5)
  c_l  = c(1)

   #-- FOR 8 Band WV data --#
   if (num_bands == 8) {
     ptm <- Sys.time()
     message(paste(Sys.time()))
     list_return = c('OCE','BLU','GRE','YEL','RED','EDG','NR1','NR2','NDVI1',
                     'NDVI2','NDRE1','NDRE2','EVI','RGI')
     OCE   =  wv_data[[1]]
     BLU   =  wv_data[[2]]
     GRE   =  wv_data[[3]]
     YEL   =  wv_data[[4]]
     RED   =  wv_data[[5]]
     EDG   =  wv_data[[6]]
     NR1   =  wv_data[[7]]
     NR2   =  wv_data[[8]]
     NDVI1 = (wv_data[[7]]-wv_data[[5]])/(wv_data[[7]]+wv_data[[5]])
     NDVI2 = (wv_data[[8]]-wv_data[[5]])/(wv_data[[8]]+wv_data[[5]])
     NDRE1 = (wv_data[[7]]-wv_data[[6]])/(wv_data[[7]]+wv_data[[6]])
     NDRE2 = (wv_data[[8]]-wv_data[[6]])/(wv_data[[8]]+wv_data[[6]])
     EVI   = c_g*(wv_data[[7]]-wv_data[[5]])/(wv_data[[7]]+c_c1*wv_data[[5]]+c_c2*wv_data[[2]]+c_l)
     RGI   = wv_data[[5]]/wv_data[[3]]
     message(paste('Calculated:'))
     message(paste(list_return,'...'))
     message(paste('Processing time:'))
     message(paste(Sys.time() - ptm))
     
     #browser() # for problem solving
     
     image = stack(OCE,BLU,GRE,YEL,RED,EDG,NR1,NR2,NDVI1,NDVI2,NDRE1,NDRE2,EVI,RGI)
     names(image) = list_return
     return(image) 
  }

  #-- For 4 Band WV DATA --#
  if (num_bands == 4) {
     message(paste(Sys.time()))
     ptm <- Sys.time()
     list_return = c('BLU','GRE','RED','NR1','NDVI1','EVI','RGI')
     BLU   =  wv_data[[1]]
     GRE   =  wv_data[[2]]
     RED   =  wv_data[[3]]
     NR1   =  wv_data[[4]]
     NDVI1 = (wv_data[[4]]-wv_data[[3]])/(wv_data[[4]]+wv_data[[3]])
     EVI = c_g*(wv_data[[4]]-wv_data[[3]])/(wv_data[[4]]+c_c1*wv_data[[3]]+c_c2*wv_data[[1]]+c_l)
     RGI   = wv_data[[3]]/wv_data[[2]]
     message(paste('Calculated:'))
     message(paste(list_return,'...'))
     message(paste('Processing time:'))
     message(paste(Sys.time() - ptm))

     return(stack(BLU,GRE,RED,NR1,NDVI1,EVI,RGI))
  }
}


# this is added for extracting data later in main code
get_value <- function(mykey, mylookupvector){
#' Match Extracted Classification Values with Key-Value Pairs (Utility Function)
#'
#' @param mykey `array` of extracted raster value
#' @param mylookupvector `list` of classes used as a dictionary for key-value pairs
#'
#' @return `array` of converted class values 
#' @export
#'
#' @examples
  myvalue <- mylookupvector[mykey]
  myvalue <- unname(myvalue)
  return(myvalue)
}

#### END ####
