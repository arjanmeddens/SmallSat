#-------------------------------------------------------------------#
#- Function to calculate different indices from:
#    NAIP, WV02, and WV03 data -#
#-------------------------------------------------------------------#
# Meddens - WSU/UIdaho - 8 Feb 2018 
#-------------------------------------------------------------------#
# Updated to include 4 band 2021 NAIP images -- 19 May 2022
# Updated to include 5 band MX RedEdge Sensor (UAS) -- 06 Jan 2023
# Updated to add `roxygen` skeleton to function -- 19 Aug 2024
#-------------------------------------------------------------------#

calc_wv_index = function(input_img,num_bands=num_bands) {
    
    #' Calculate index for WorldView2, NAIP, MX Rededge imagery
    #'
    #' @param input_img input imagery
    #' @param num_bands number of bands of the `input_img`. `num_bands = 4` for NAIP imagery, `num_bands = 5` for MX RedEdge Sensor (UAS), and `num_bands = 8` for WorldView 2 imagery
    #'
    #' @return input_imagery with calculated indices and renamed bands.
    #'   Bands returned (in order) for:
    #'   * NAIP imagery input: `c('BLU','GRE','RED','NR1','NDVI1','EVI','RGI')`
    #'   * MX RedEdge input: `c('BLU','GRE','RED', 'REDEDGE', 'NIR', 'NDVI', 'NDRE', 'SR', 'GLI', 'RGI')`
    #'   * WorldView imagery input: `c('OCE','BLU','GRE','YEL','RED','EDG','NR1','NR2','NDVI1','NDVI2','NDRE1','NDRE2','EVI','RGI')`
    #' @export
    #'
    #' @examples
    
    
  if (class(input_img) == "RasterBrick") {
    #print(“input_img not specified”)
  } else {
    message(paste('input_img should be a RasterBrick'))
  }

  if (length(num_bands) == 0) {
    num_bands == dims(input_img)[3]
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
     
     #browser()
  
     ptm <- Sys.time()
     message(paste(Sys.time()))
     list_return = c('OCE','BLU','GRE','YEL','RED','EDG','NR1','NR2','NDVI1',
                     'NDVI2','NDRE1','NDRE2','EVI','RGI')
     OCE   =  input_img[[1]]
     BLU   =  input_img[[2]]
     YEL   =  input_img[[3]]
     GRE   =  input_img[[4]]
     RED   =  input_img[[5]]
     EDG   =  input_img[[6]]
     NR1   =  input_img[[7]]
     NR2   =  input_img[[8]]
     NDVI1 = (input_img[[7]]-input_img[[5]])/(input_img[[7]]+input_img[[5]])
     NDVI2 = (input_img[[8]]-input_img[[5]])/(input_img[[8]]+input_img[[5]])
     NDRE1 = (input_img[[7]]-input_img[[6]])/(input_img[[7]]+input_img[[6]])
     NDRE2 = (input_img[[8]]-input_img[[6]])/(input_img[[8]]+input_img[[6]])
     EVI   = c_g*(input_img[[7]]-input_img[[5]])/(input_img[[7]]+c_c1*input_img[[5]]+c_c2*input_img[[2]]+c_l)
     RGI   = input_img[[5]]/input_img[[4]]
     message(paste('Calculated:'))
     message(paste(list_return,'...'))
     message(paste('Processing time:'))
     message(paste(Sys.time() - ptm))
     
     image = stack(OCE,BLU,YEL,GRE,RED,EDG,NR1,NR2,NDVI1,NDVI2,NDRE1,NDRE2,EVI,RGI)
     names(image) = list_return
     return(image) 
     #img = brick(OCE,BLU,YEL,GRE,RED,EDG,NR1,NR2,NDVI1,NDVI2,NDRE1,NDRE2,EVI,RGI)
     #return(brick(list(OCE=OCE,BLU=BLU,YEL=YEL,GRE=GRE,RED=RED,EDG=EDG,NR1=NR1,NR2=NR2,
		 #NDVI1=NDVI1,NDVI2=NDVI2,NDRE1=NDRE1,NDRE2=NDRE2,EVI=EVI,RGI=RGI)))
  }

  #-- For 4 Band WV/NAIP DATA --#
  if (num_bands == 4) {
     message(paste(Sys.time()))
     ptm <- Sys.time()
     list_return = c('BLU','GRE','RED','NR1','NDVI1','EVI','RGI')
     BLU   =  input_img[[1]]
     GRE   =  input_img[[2]]
     RED   =  input_img[[3]]
     NR1   =  input_img[[4]]
     NDVI1 = (input_img[[4]]-input_img[[3]])/(input_img[[4]]+input_img[[3]])
     EVI = c_g*(input_img[[4]]-input_img[[3]])/(input_img[[4]]+c_c1*input_img[[3]]+c_c2*input_img[[1]]+c_l)
     RGI   = input_img[[3]]/input_img[[2]]
     message(paste('Calculated:'))
     message(paste(list_return,'...'))
     message(paste('Processing time:'))
     message(paste(Sys.time() - ptm))
     return(raster::brick(list(BLU=BLU,GRE=GRE,RED=RED,NR1=NR1,NDVI1=NDVI1,EVI=EVI,RGI=RGI)))
  }

  #-- For 5 Band RedEdge MX (UAV MS sensor) DATA --#
  if (num_bands == 5) {

    img = input_img
    tmp_img <- (img/10000) 
    
    ptm <- Sys.time()
    list_return = c('BLU','GRE','RED', 'REDEDGE', 'NIR', 'NDVI', 'NDRE', 'SR', 'GLI', 'RGI') # bands 
    BLU     =  tmp_img[[1]]
    GRE     =  tmp_img[[2]]
    RED     =  tmp_img[[3]]
    REDEDGE =  tmp_img[[4]]
    NIR     =  tmp_img[[5]]
    NDVI    = (tmp_img[[5]]-tmp_img[[3]]) / (tmp_img[[5]]+tmp_img[[3]]) 
    NDRE    = (tmp_img[[5]]-tmp_img[[4]]) / (tmp_img[[5]]+tmp_img[[4]]) 
    SR      = (tmp_img[[5]]/tmp_img[[3]])
    GLI     = ((tmp_img[[2]]-tmp_img[[3]]) + (tmp_img[[2]]-tmp_img[[1]])) / ((2*tmp_img[[2]])+tmp_img[[3]]+tmp_img[[1]])   
    RGI     = (tmp_img[[3]]/tmp_img[[2]])
    message(paste('Calculated:'))
    message(paste(list_return,'...'))
    message(paste('Processing time:'))
    message(paste(Sys.time() - ptm))
    if (class(tmp_img) == "RasterBrick") {
      return(raster::brick(list(BLU=BLU,GRE=GRE,RED=RED,REDEDGE=REDEDGE,NIR=NIR, NDVI=NDVI, NDRE=NDRE, SR=SR, GLI=GLI, RGI=RGI)))
    } else {
      return(raster::brick(list(BLU=BLU,GRE=GRE,RED=RED,REDEDGE=REDEDGE,NIR=NIR, NDVI=NDVI, NDRE=NDRE, SR=SR, GLI=GLI, RGI=RGI)))
    }
  }




}

#### END ####
