################################################################################################
# Purpose: Workflow to classify high-res SmallSAT data
################################################################################################
# Author: Arjan Meddens 
# Date: Dec, 2023
# R version & packages: Raster, sp, etc (see p_load) 
# - Note I moved back to Raster after Terra gave too many bugs.  
################################################################################################
# Input format: .tif (WV spectral data) & .shp (identified classes, point data), & 
#               .csv (run params file)
# Output format: .tif (classified imagery), figures (pdf), output data (csv, txt) 
################################################################################################
# clean environment, start code timer
rm(list = ls()); gc()
ptm <- Sys.time()
message(paste(Sys.time()))

#--- Set up directory structure -------
dir = getwd()

dataPATH = file.path(dir, "DATA") #- input data path
exportsPATH = file.path(dir, "EXPORTS") #- output data path

# note: use of `file.path()` because it creates folder paths for:
# - windows-based systems with "\" 
# - unix (Mac) or linux-based systems with "/"  

# Input paths (ensure files are in DATA folder)
shp_dir = dataPATH
img_dir = dataPATH
img_cloud_dir = dataPATH

# Output paths
fig_dir = exportsPATH 
outdir  = exportsPATH



#------------------------------------------------------------------------
# Load packages & Functions
#------------------------------------------------------------------------
library(pacman)             #- To use p_load for faster package loading
p_load(raster,GMCM,rgl,rgdal,randomForest,rfUtilities,maptools,sp,spatial,RColorBrewer,
       ggplot2,caret,e1071,RStoolbox,plotrix,mapaccuracy,doParallel,graph4lg,utils,kernlab)
registerDoParallel(cores=1)
source("calc_wv_index.r") # make sure this is in working dir
get_value <- function(mykey, mylookupvector){
  myvalue <- mylookupvector[mykey]
  myvalue <- unname(myvalue)
  return(myvalue)
}
#------------------------------------------------------------------------
#- Run options 
#------------------------------------------------------------------------
plot = 1                  # Flag to plot imagery etc. 0:no/1:yes
flag_varsel = 1           # Flag to run variable selection
check_maj = 0             # Flag to run majority filter (turned off for now)
run.list = c(1:1)        # Enter list to be run, taken from input parameter file:
                          #  "smallsat_run_params_231022.csv", Enter "c(1:1)" for running single images     
message(paste("Entered argument:....",run.list[1]))

#--- Read input param file ------------ 
input_arr = suppressWarnings(read.csv(file.path(shp_dir,"smallsat_run_params_231221.csv")))
img_no_arr        = input_arr[,1]
img_id_arr        = input_arr[,2]
shp_id_arr        = input_arr[,3]   # Example: 'GE01_20210824_M1BS_105001002694C800_eval6_UTM_p50.shp' 
ext_shp_file_arr  = input_arr[,4]   # Example: extent_mt_img1_utm.shp"
fnf_mask_file_arr = input_arr[,5]   # Example: c("/Users/arjanmeddens/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/MeddensLab/3_NASA_SmallSat/11_FNF_mask/FNF_site2_mt_clip.tif")
cld_shadow_arr    = input_arr[,6]   # Example: Masked at band: 3 
cld_trhd_arr      = input_arr[,7]   # Example: Masked at band: 2
no_bands_arr      = input_arr[,8]   # Example: 
seed_arr          = input_arr[,9]   # Example: 
img_cloud_arr     = input_arr[,15]  # Example: 
varsel_arr        = input_arr[,32]  # Example: 
cor.var.arr       = input_arr[,33]  # Example: 0.925
num_run_img = length(input_arr[,1])

################################################################################################
# Starting for loop
################################################################################################
for (x in 1:length(run.list)) { # here we input with the Args()
  imgno = run.list[x]
  img_id              = img_id_arr[imgno] # here we can loop over CATID...
  shp_id              = shp_id_arr[imgno]
  ext_shp_file        = ext_shp_file_arr[imgno]
  fnf_mask_file       = fnf_mask_file_arr[imgno]
  cloud_mask_file     = img_cloud_arr[imgno]
  cld_shadow_trhd     = cld_shadow_arr[imgno]
  cld_trhd            = cld_trhd_arr[imgno]
  no.bands.class      = no_bands_arr[imgno]
  seed                = seed_arr[imgno]
  varsel.treshold     = varsel_arr[imgno]
  cor.var.threshold   =cor.var.arr[imgno]
  run_id = paste("_s",seed,sep='') # name of batch
  
  message(paste("#----------------------------------------------------------"))
  message(paste("Running.......",img_id,"....Image number:",imgno))
  message(Sys.time())
  message(paste("Running.......",shp_id))
  message(paste("#----------------------------------------------------------"))
  
  #---------------------------------------------------------------------------------------------  
  #- Load Eval/Training shape file, imagery, and calculate indices
  #---------------------------------------------------------------------------------------------
  shp = raster::shapefile(file.path(shp_dir,shp_id))
  names(shp)[2]  = c("class")
  levels.class = levels(as.factor(shp$class))
  if (length(which(shp$class == "yt")) > 1){
    n.class = length(unique(shp$class))
    paste("Gray tree class, number of classes:",n.class)
  } else {
    n.class = length(unique(shp$class))
    paste("No grey tree class, number of classes:",n.class)
  }
  #- Load image with indices or create indices if image w/indices does not exist
  file.name = (file.path(img_dir,paste(img_id,"_index.tif",sep='')))
  if (file.exists(file.name)==1) {
    message(paste("File with Indices exists...skipping:",file.name)) 
    img = raster::brick(paste(file.name,sep=''))
  } else {
    message(paste("File with Indices does not exists...creating:",file.name)) 
    tmp_img = raster::brick(file.path(img_dir, paste(img_id,".tif",sep='')))
    #-- Calculate indices
    num_bands = length(names(tmp_img))
    img2 = calc_wv_index(tmp_img,num_bands=num_bands)
    img = img2
    # img = raster::brick(img2) # WHY brick twice?
    raster::writeRaster(img, filename=file.name,overwrite=TRUE)
  }

  # Setting the number of spectral bands to plot 
  n.bands=length(names(img))
  if (n.bands == 14) {
    n.spec = 8 
    band1 = 5 # plot which bands? 
    band2 = 9 # plot which bands? 
    nir_band = 7
  } else {
    n.spec = 4
    band1 = 3 # plot which bands? 
    band2 = 5 # plot which bands? 
    nir_band = 4
  }

  #---------------------------------------------------------------------------------------------
  # Select extent
  #---------------------------------------------------------------------------------------------
  if (ext_shp_file != "") {
    shp_ext_id = file.path(shp_dir,ext_shp_file)
    ext_shp = raster::extent(shapefile(shp_ext_id))
  } else {
    ext_shp = raster::extent(img)  # Running entire extent of input image
  }

  #---------------------------------------------------------------------------------------------
  # Create training and evaluation data (randomly select 50%, 25 pixels per set)
  # for each class (Green tree, red tree, gray tree (if present), herbaceous (green and brown), 
  # bare ground, shadow.
  #---------------------------------------------------------------------------------------------
  factor(shp$class)
  shp.gt = shp[c(which((shp$class == "gt") == TRUE)), ]
  shp.rt = shp[c(which((shp$class == "rt") == TRUE)), ]
  shp.hb = shp[c(which((shp$class == "hb") == TRUE)), ]
  shp.bg = shp[c(which((shp$class == "bg") == TRUE)), ]
  shp.sh = shp[c(which((shp$class == "sh") == TRUE)), ]
  sampleSize = floor(length(shp.bg[,1])/2) # NOTE changed from rt
  sam.train.gt <- shp.gt[sample(1:length(shp.gt),sampleSize),]
  sam.eval.gt  <- shp.gt[-sample(1:length(shp.gt),sampleSize),]
  sam.train.rt <- shp.rt[sample(1:length(shp.rt),sampleSize),]
  sam.eval.rt  <- shp.rt[-sample(1:length(shp.rt),sampleSize),]
  sam.train.hb <- shp.hb[sample(1:length(shp.hb),sampleSize),]
  sam.eval.hb  <- shp.hb[-sample(1:length(shp.hb),sampleSize),]
  sam.train.bg <- shp.bg[sample(1:length(shp.bg),sampleSize),]
  sam.eval.bg  <- shp.bg[-sample(1:length(shp.rt),sampleSize),]
  sam.train.sh <- shp.sh[sample(1:length(shp.sh),sampleSize),]
  sam.eval.sh  <- shp.sh[-sample(1:length(shp.sh),sampleSize),]

  # If there is a grey tree class #
  if (n.class == 6) {
    shp.yt = shp[c(which((shp$class == "yt") == TRUE)), ]
    sam.train.yt <- shp.yt[sample(1:length(shp.yt),sampleSize),] 
    sam.eval.yt  <- shp.yt[-sample(1:length(shp.yt),sampleSize),]
    n.train = c(length(sam.train.gt),length(sam.train.rt),length(sam.train.hb),length(sam.train.bg),length(sam.train.sh),length(sam.train.yt))
    n.eval = c(length(sam.eval.gt),length(sam.eval.rt),length(sam.eval.hb),length(sam.eval.bg),length(sam.eval.sh),length(sam.eval.yt))
    trainData <- rbind(sam.train.gt, sam.train.rt,sam.train.hb,sam.train.bg,sam.train.sh,sam.train.yt)
    evalData  <- rbind(sam.eval.gt,  sam.eval.rt, sam.eval.hb, sam.eval.bg, sam.eval.sh, sam.eval.yt)
  } else {
    n.train = c(length(sam.train.gt),length(sam.train.rt),length(sam.train.hb),length(sam.train.bg),length(sam.train.sh))
    n.eval = c(length(sam.eval.gt),length(sam.eval.rt),length(sam.eval.hb),length(sam.eval.bg),length(sam.eval.sh))
    trainData <- rbind(sam.train.gt, sam.train.rt,sam.train.hb,sam.train.bg,sam.train.sh)
    evalData  <- rbind(sam.eval.gt,  sam.eval.rt, sam.eval.hb, sam.eval.bg, sam.eval.sh)
  }

  #--- Extract values for each of the training/evaluation data points
  n.bands = length(names(img))
  if (n.bands == 7) {band_names_arr <- c("BLU","GRE","RED","NR1","NDVI1","EVI","RGI")}
  if (n.bands == 14) {band_names_arr <- c("OCE","BLU","YEL","GRE","RED","EDG","NR1","NR2","NDVI1","NDVI2","NDRE1","NDRE2","EVI","RGI")}
  class_arr = unique(shp$class)
  n.bands = length(names(img))
  class_band_arr = matrix(data=NA,nrow=10000,ncol=100)
  train.arr = matrix(data=NA,nrow=sum(n.train),ncol=n.bands+1)
  eval.arr  = matrix(data=NA,nrow=sum(n.eval),ncol=n.bands+1)
  for (b in 1:n.bands) {
    print(paste("Band:",b))
    if (n.class == 6) {
      train.arr[,b] =  c(extract(img[[b]],sam.train.gt),extract(img[[b]],sam.train.rt),extract(img[[b]],sam.train.hb),
                         extract(img[[b]],sam.train.bg),extract(img[[b]],sam.train.sh),extract(img[[b]],sam.train.yt))
      eval.arr[,b] =  c(extract(img[[b]],sam.eval.gt),extract(img[[b]],sam.eval.rt),extract(img[[b]],sam.eval.hb),
                        extract(img[[b]],sam.eval.bg),extract(img[[b]],sam.eval.sh),extract(img[[b]],sam.eval.yt))
    } else {
      train.arr[,b] =  c(extract(img[[b]],sam.train.gt),extract(img[[b]],sam.train.rt),extract(img[[b]],sam.train.hb),
                         extract(img[[b]],sam.train.bg),extract(img[[b]],sam.train.sh))
      eval.arr[,b] =  c(extract(img[[b]],sam.eval.gt),extract(img[[b]],sam.eval.rt),extract(img[[b]],sam.eval.hb),
                        extract(img[[b]],sam.eval.bg),extract(img[[b]],sam.eval.sh))
    }
  }

  if (n.class == 6) {
    train.arr[1                   :cumsum(n.train)[1],n.bands+1] = 1       # class 1: gt
    train.arr[c(cumsum(n.train)[1]+1):c(cumsum(n.train)[2]),n.bands+1] = 2 # class 2: rt
    train.arr[c(cumsum(n.train)[2]+1):c(cumsum(n.train)[3]),n.bands+1] = 3 # class 3: bg
    train.arr[c(cumsum(n.train)[3]+1):c(cumsum(n.train)[4]),n.bands+1] = 4 # class 1: hb
    train.arr[c(cumsum(n.train)[4]+1):c(cumsum(n.train)[5]),n.bands+1] = 5 # class 1: sh
    train.arr[c(cumsum(n.train)[5]+1):c(cumsum(n.train)[6]),n.bands+1] = 6 # class 1: yt
    eval.arr[1                   :cumsum(n.eval)[1],n.bands+1] = 1      # class 1: gt
    eval.arr[c(cumsum(n.eval)[1]+1):c(cumsum(n.eval)[2]),n.bands+1] = 2 # class 2: rt
    eval.arr[c(cumsum(n.eval)[2]+1):c(cumsum(n.eval)[3]),n.bands+1] = 3 # class 3: bg
    eval.arr[c(cumsum(n.eval)[3]+1):c(cumsum(n.eval)[4]),n.bands+1] = 4 # class 4: hb
    eval.arr[c(cumsum(n.eval)[4]+1):c(cumsum(n.eval)[5]),n.bands+1] = 5 # class 5: sh
    eval.arr[c(cumsum(n.eval)[5]+1):c(cumsum(n.eval)[6]),n.bands+1] = 6 # class 6: yt
  } else {
    train.arr[1                   :cumsum(n.train)[1],n.bands+1] = 1       # class 1: gt
    train.arr[c(cumsum(n.train)[1]+1):c(cumsum(n.train)[2]),n.bands+1] = 2 # class 2: rt
    train.arr[c(cumsum(n.train)[2]+1):c(cumsum(n.train)[3]),n.bands+1] = 3 # class 3: bg
    train.arr[c(cumsum(n.train)[3]+1):c(cumsum(n.train)[4]),n.bands+1] = 4 # class 1: hb
    train.arr[c(cumsum(n.train)[4]+1):c(cumsum(n.train)[5]),n.bands+1] = 5 # class 1: sh
    eval.arr[1                   :cumsum(n.eval)[1],n.bands+1] = 1      # class 1: gt
    eval.arr[c(cumsum(n.eval)[1]+1):c(cumsum(n.eval)[2]),n.bands+1] = 2 # class 2: rt
    eval.arr[c(cumsum(n.eval)[2]+1):c(cumsum(n.eval)[3]),n.bands+1] = 3 # class 3: bg
    eval.arr[c(cumsum(n.eval)[3]+1):c(cumsum(n.eval)[4]),n.bands+1] = 4 # class 4: hb
    eval.arr[c(cumsum(n.eval)[4]+1):c(cumsum(n.eval)[5]),n.bands+1] = 5 # class 5: sh
  }

  #------------------------------------------------------------------------
  # (0) Create a random forest model to list important variables
  #------------------------------------------------------------------------
  train.arr = data.frame(train.arr) 
  eval.arr  = data.frame(eval.arr) 
  colnames(train.arr) <-c(band_names_arr,"Class")
  colnames(eval.arr)  <-c(band_names_arr,"Class")
  train.arr$Class <- factor(train.arr$Class)
  eval.arr$Class  <- factor(eval.arr$Class)

  which(is.na(factor(train.arr$Class))==TRUE)
  rf.classifier = randomForest::randomForest(factor(train.arr$Class) ~ ., type = classification,
                                              data = train.arr, 
                                              ntree = 500, 
                                              importance = TRUE,na.action=na.exclude)
  print(rf.classifier)
  pred.arr = randomForest::predict(rf.classifier,eval.arr) 
  table(Observed=eval.arr[,n.bands+1],Predicted=pred.arr)
  class_txt_arr = c("Healthy trees","Red trees","Herbaceous","Bare ground","Shadow","Gray trees")
  class_txt_arr
  accuracy(eval.arr[,n.bands+1],pred.arr)
  #index_list_arr = names(sort(importance(rf.classifier)[,6],decreasing = T))
  index_list_arr = rownames(randomForest::importance(rf.classifier))

  #------------------------------------------------------------------
  # Inset "break" for testing the variable selection
  #------------------------------------------------------------------
  #if (x == 1) {  
  #  print("breaking for loop for testing")
  #  break
  #  } 

  #------------------------------------------------------------------
  # RF variable selection 
  #------------------------------------------------------------------
  if (flag_varsel == 1) {
    iterations<-1
    VarSel<-matrix(ncol=length(index_list_arr), nrow=iterations)
    #bootstrapp the rf.modelSel
    for (i in 1:iterations){
      seed_n<-sample(1:1000, 1, replace=FALSE)  
      fit_msel<- rfUtilities::rf.modelSel(train.arr[1:n.bands],train.arr$Class, imp.scale = "mir",
                            final.model = FALSE, seed = seed_n, ntree=1000, parsimony=0.03)
      varimp<-fit_msel$imp
      vardat<-varimp[,1]
      df<-t(vardat)
      VarSel[i,] <- c(df)
      names(vardat) <- rownames(varimp)
      Sys.sleep(0.1)
      sort.vardat = (sort(vardat,decreasing = T))
    }
    ######################
    # Rules for selection: 
    # 1. r<0.95, choose one: Here Red and Blue are super correlated, so Red goes.
    # 2. rf.modelSel with Varsel include > 0.7
    ######################
    train.check.arr = train.arr[1:n.bands]
    for (v in 1:20){
      if (v > length(sort.vardat)) { 
        print("reached max variables... Breaking out of loop...")
        break; 
      }
      index.var = which(colnames(train.check.arr) == names(sort.vardat[v])) 
      #--- Selecting only variables that are correlated less than input threshold, eg: which(cor.var.arr < 0.9)
      index.var.cor = which(cor((train.check.arr))[,index.var] > cor.var.threshold & cor((train.check.arr))[,index.var] != c(1.0))
      print(paste("Checking",names(sort.vardat[v])))
      if (length(index.var.cor[]) == 0) {
        print("no problems next variable")
      } else {  
        remove.colums.train = which(names(train.check.arr) %in% names(index.var.cor))
        train.check.arr = train.check.arr[,-c(remove.colums.train)]
        remove.colums.vardat = which(names(sort.vardat) %in% names(index.var.cor))
        sort.vardat     =     sort.vardat[-c(remove.colums.vardat)]
        print("problems cull some variables")
        print(paste("Culling ",names(index.var.cor)))
        print(paste("sort.vardat",sort.vardat))
        print(paste("train.check.arr",names(train.check.arr)))
      }  
    }
  }
  cor.train.mtx = cor((train.arr[1:n.bands]))
  #include<-VarSel>varsel.treshold
  include2<-sort.vardat>varsel.treshold
  index.selected = which(include2 == T)
  select.arr = include2[index.selected]
  #-----------------------------------------------
  index_indices = which(index_list_arr %in% names(select.arr))
  no.bands.class = length(index_indices)
  index_list_arr = index_list_arr[index_indices] #summary.rf.selboot[index_indices,1]

  #------------------------------------------------------------------------
  # Create a (1) maximum likelihood classification and 
  #  (2) randomForest classification and (3) NN with selected bands.
  #------------------------------------------------------------------------
  index_arr = array(data=NA,dim=no.bands.class)
  for (i in 1:no.bands.class) {
    index_arr[i] = which(band_names_arr == index_list_arr[i])
  }
  tmp_img = subset(img,index_arr)
  if (ext_shp_file != "") {
    tmp_img = crop(tmp_img,ext_shp) # .... Note all training data needs to be in shp_ext!!!
  }

  SC.mlhc       <- RStoolbox::superClass(tmp_img, trainData = trainData, responseCol = "class", 
                       model = "mlc", verbose=1)
  val.mlhc  <- RStoolbox::validateMap(SC.mlhc$map, valData = evalData, responseCol = "class", 
                       classMapping = SC.mlhc$classMapping)
  SC.rf       <- RStoolbox::superClass(tmp_img, trainData = trainData, responseCol = "class", 
                       model = "rf", verbose=1)
  val.rf  <- RStoolbox::validateMap(SC.rf$map, valData = evalData, responseCol = "class", 
                         classMapping = SC.rf$classMapping)

  #https://urldefense.com/v3/__https://rpubs.com/uky994/593668__;!!JYXjzlvb!hOLHnMhmrMM9FiGDneDIoKNWZSgAdlCDbAFaml3Ql0HdLPd4ooC8gHNGqWrubAdfPL4IR6bNULMH9K6525NzkgkZ8XCNWqcFWcA$ 
  #https://urldefense.com/v3/__https://stats.stackexchange.com/questions/73032/linear-kernel-and-non-linear-kernel-for-support-vector-machine__;!!JYXjzlvb!hOLHnMhmrMM9FiGDneDIoKNWZSgAdlCDbAFaml3Ql0HdLPd4ooC8gHNGqWrubAdfPL4IR6bNULMH9K6525NzkgkZ8XCNYZ3319g$  
  #https://urldefense.com/v3/__https://blog.revolutionanalytics.com/2015/10/the-5th-tribe-support-vector-machines-and-caret.html__;!!JYXjzlvb!hOLHnMhmrMM9FiGDneDIoKNWZSgAdlCDbAFaml3Ql0HdLPd4ooC8gHNGqWrubAdfPL4IR6bNULMH9K6525NzkgkZ8XCN2e6UvA8$  
  SC.svm <- RStoolbox::RStoolbox::superClass(tmp_img, trainData = trainData, responseCol = "class", 
                     model = "svmRadial",            # linear kernel
                     tuneGrid = expand.grid(sigma = c(.01, .015, 0.2),
                     C = seq(0.1, 1.5, length = 15)),	
                     preProc = c("center","scale"))  # Center and scale data
  val.svm  <- RStoolbox::validateMap(SC.svm$map, valData = evalData, responseCol = "class", 
                        classMapping = SC.svm$classMapping)

  #------------------------------------------------------------------------
  #- Develop mask to mask clouds and clouds shadows
  #------------------------------------------------------------------------
  file.name = (paste(img_cloud_dir,cloud_mask_file,sep=''))
  if (file.exists(file.name)==1) {
    message(paste("Cloud images exists...opening:",file.name)) 
    cloud_img = raster::raster(paste(file.name,sep=''))

    if (ext_shp_file != "") {   # Crop Cloud img to same extent!!!
      cloud_img = raster::crop(cloud_img,ext_shp)
    }  # Cloud shadow
    index_mask1 = cloud_img < 0.5   # Water/cloud shadow
    index_mask2 = cloud_img < 0.5   # Water/cloud shadow
  } else {
    # Cloud shadow
    img_crop = raster::crop(img[[3]],ext_shp)
    index_mask1 = img_crop > cld_shadow_trhd   # Water/cloud shadow
    # Clouds/bright objects
    img_crop = raster::crop(img[[2]],ext_shp)
    index_mask2 = img_crop < cld_trhd          # Clouds   
  }

  #########################################################################
  #- ADD IN IF FNF!!! - See older code how to create fnf mask
  #########################################################################
  if (file.exists(fnf_mask_file)==1) {
    img_crop = raster::crop(fnf_mask,ext_shp)
    index_mask3 = img_crop < 3   # FNF radar mask --> 3:Non-for// 4:water
    mask_all = index_mask1 * index_mask2 * index_mask3 #Including FNF mask
  } else {
    mask_all = index_mask1 * index_mask2               #Excluding FNF mask
  }

  #--- Building the entire mask
  index_mask = which(mask_all[] == 0)
  SC.mlhc$map[index_mask] = 5
  SC.rf$map[index_mask]   = 5 
  SC.svm$map[index_mask]  = 5 

  #---------------------------------------------------------------------------------------------
  ## Calculating pixel %'s --> Note after Masking!!! Note all masking is set to shadow 
  #---------------------------------------------------------------------------------------------
  class_list = unlist(levels(factor(evalData$class)))
  n.pixels.mask = length(index_mask)
  n.pixels = length(index_mask1[])
  n.pixels = length(which(is.na(SC.mlhc$map[]) == FALSE))
  n.class.mlhc = array(array(data = NA, dim = n.class+2))
  n.class.rf = array(array(data = NA, dim = n.class+2))
  n.class.svm = array(array(data = NA, dim = n.class+2))
  class.code.arr = c(0,1,2,3,4,5,6,NA)
  # Mask(0) # Bare Ground (1) # Green Trees (2) # Herbaceous (3) # Red trees (4) # Shadow (5) # Gray trees (6)  # EXtra class (7) # No data
  for (i in 1:(n.class+3)) {
    n.class.mlhc[i] =  length(which(values(SC.mlhc$map) %in% class.code.arr[i])) 
    n.class.rf[i]   =  length(which(values(SC.rf$map) %in% class.code.arr[i])) # Mask # Bare Ground # Green Trees # Herbaceous # Red trees # Shadow # Gray trees
    n.class.svm[i]  =  length(which(values(SC.svm$map) %in% class.code.arr[i])) # Mask # Bare Ground # Green Trees # Herbaceous # Red trees # Shadow # Gray trees
  } 
  p.class.mlhc = (100*n.class.mlhc)/n.pixels
  p.class.rf   = (100*n.class.rf)  /n.pixels
  p.class.svm   = (100*n.class.svm)  /n.pixels
  area.img     = (n.pixels * 2 * 2)/1e6        # km2
  area.na      = (n.class.mlhc[8]* 2 * 2)/1e6  # Nodata land: is equal for each image
  area.mask   = (n.class.mlhc[1]* 2 * 2)/1e6  # Mask: is equal for each image
  area.data    =  area.img -(area.na+area.mask)
  p.area       = 100*c(area.img,area.na,area.mask,area.data)/area.img

  if (n.class == 6) {n=7} else {n=6}
  n.class.mlhc.arr = n.class.mlhc[2:n]
  names(n.class.mlhc.arr) = levels.class
  n.class.rf.arr = n.class.rf[2:n]
  names(n.class.rf.arr) = levels.class
  n.class.svm.arr = n.class.svm[2:n]
  names(n.class.svm.arr) = levels.class

  #---------------------------------------------------------------------------------------------
  #### Accuracy metrics (including the mask) # if a class is omitted olofsson will give an error!!!!
  ## See: https://urldefense.com/v3/__https://blogs.fu-berlin.de/reseda/accuracy-statistics-in-r/__;!!JYXjzlvb!hOLHnMhmrMM9FiGDneDIoKNWZSgAdlCDbAFaml3Ql0HdLPd4ooC8gHNGqWrubAdfPL4IR6bNULMH9K6525NzkgkZ8XCNXDohsk0$ 
  #---------------------------------------------------------------------------------------------
  #https://urldefense.com/v3/__http://www.cookbook-r.com/Data_input_and_output/Writing_text_and_output_from_analyses_to_a_file/__;!!JYXjzlvb!hOLHnMhmrMM9FiGDneDIoKNWZSgAdlCDbAFaml3Ql0HdLPd4ooC8gHNGqWrubAdfPL4IR6bNULMH9K6525NzkgkZ8XCNNX0EIjI$ 
  #out_info_file = paste(outdir, img_id,run_id,".txt",sep='')
  #sink(out_info_file)
  pred.class.mlhc  = raster::extract(SC.mlhc$map,evalData)
  pred.class.mlhc  = get_value((pred.class.mlhc),c(class_list))
  olof.table.mlhc  = mapaccuracy::olofsson(as.character(pred.class.mlhc),as.character(evalData$class),n.class.mlhc.arr) # to remove 6th class/here NA
  accuracy.mlhc    = accuracy(pred.class.mlhc,as.character(evalData$class))
  error.mtx.mlhc  = confusionMatrix((as.factor(pred.class.mlhc)),as.factor(evalData$class), mode = "everything", positive="1")

  pred.class.rf = raster::extract(SC.rf$map,evalData)
  pred.class.rf = get_value((pred.class.rf),c(class_list))
  olof.table.rf = olofsson(pred.class.rf,as.character(evalData$class),n.class.rf.arr)
  accuracy.rf   = accuracy(pred.class.rf,as.character(evalData$class))
  error.mtx.rf  = confusionMatrix((as.factor(pred.class.rf)),as.factor(evalData$class), mode = "everything", positive="1")

  pred.class.svm = raster::extract(SC.svm$map,evalData)
  pred.class.svm = get_value((pred.class.svm),c(class_list))
  olof.table.svm = olofsson(pred.class.svm,as.character(evalData$class),n.class.svm.arr)
  accuracy.svm   = accuracy(pred.class.svm,as.character(evalData$class))
  error.mtx.svm  = confusionMatrix((as.factor(pred.class.svm)),as.factor(evalData$class), mode = "everything", positive="1")

  ##############################################################################
  #- Plotting 
  ##############################################################################
  #----------------------------
  # (A). Variable importance plot 
  #----------------------------
  raster::importance(rf.classifier)
  if (plot) {
    png(filename = file.path(fig_dir,paste("varImpRF_",img_id,run_id,".png",sep="")),width = 480, height = 480, 
        units = "px", pointsize = 12,bg = "white")
        raster::varImpPlot(rf.classifier)
    dev.off()
  }

  #----------------------------
  # (B1). Class spectra 
  #----------------------------
  if (plot) {
    png(filename = file.path(fig_dir,paste("spectra_",img_id,run_id,".png",sep="")),width = 1380, height = 680, 
        units = "px", pointsize = 12,bg = "white")
    par(mfrow = c(1, 2),cex = 2,mar = c(3, 4, 0, 0), oma = c(0.2, 0.2, 0.5, 0.5),tcl = -0.25,mgp = c(2, 0.6, 0))
    ave_spec_arr = matrix(data=NA,nrow=n.class,ncol=n.spec)
    std_spec_arr = matrix(data=NA,nrow=n.class,ncol=n.spec)
    x.val=c(1:n.spec)
    col_arr = c("dark green","red","light green","darkgoldenrod1","black","dark gray")
    for (c in 1:n.class) {
      index_class = which(train.arr[,n.bands+1] == c)
      tmp_arr = train.arr[index_class,c(1:n.spec)]
      ave_spec_arr[c,] = colMeans(tmp_arr)
      std_spec_arr[c,] = apply(tmp_arr,2,sd)
      if (c == 1) {
        plot(ave_spec_arr[c,],type="l",col=col_arr[c],ylim=c(0,1*max(train.arr[,nir_band])),
             ylab="Spectral value",xlab="Band",lwd=c(6))
        arrows(x.val,(ave_spec_arr[c,]+(std_spec_arr[c,])),x.val,(ave_spec_arr[c,]-(std_spec_arr[c,])),
               col=col_arr[c],length=0,lwd=1.5)
      } else {
        lines(x.val+(randu[c+5,1]/10)-0.05,ave_spec_arr[c,],col=col_arr[c],lwd=c(6))
        arrows(x.val+(randu[c+5,1]/10)-0.05,(ave_spec_arr[c,]+(std_spec_arr[c,])),x.val+(randu[c+5,1]/10)-0.05,(ave_spec_arr[c,]-(std_spec_arr[c,])),
               col=col_arr[c],length=0,lwd=1.5)
      }
    }
    class_txt_arr = c("Healthy trees","Red trees","Herbaceous","Bare ground","Shadow","Gray trees")
    legend(1,max(train.arr[,nir_band]),class_txt_arr,col=col_arr,lty=c(1,1,1,1,1),lwd=c(5,5,5,5,5),cex=0.8)
    #----------------------------
    # (B2). Scatter plot 
    #----------------------------
    min.xy = c(min(train.arr[,band1]),min(train.arr[,band2]))
    max.xy = c(max(train.arr[,band1]),max(train.arr[,band2]))
    for (c in 1:n.class) {
      index_class = which(train.arr[,n.bands+1] == c)
      tmp_arr = train.arr[index_class,]
      if (c==1) {
        plot(tmp_arr[,band1],tmp_arr[,band2],col=col_arr[c],ylim=c(min.xy[2],max.xy[2]),
             xlim=c(min.xy[1],max.xy[1]),xlab=band_names_arr[band1],ylab=band_names_arr[band2],
             cex =1, lwd = 3)
      } else {
        points(tmp_arr[,band1],tmp_arr[,band2],col=col_arr[c],cex = 1, lwd = 3)
      }
    }
    dev.off()
  }

  #----------------------------
  # (C). Plotting MASKS images
  #----------------------------
  if (plot) {
    png(filename = file.path(fig_dir,paste("mask_",img_id,run_id,".png",sep="")),width = 1380, height = 680, 
        units = "px", pointsize = 12,bg = "white")
    par(mfrow = c(1, 3),cex = 1,mar = c(1, 2, 4, 4), oma = c(1.8, 0, 1.5, 1.2),tcl = -0.25,mgp = c(2, 0.6, 0))
    plot(index_mask1,main=c("Mask shadow/water"), zlim=c(-1,1))
    plot(index_mask2,main=c("Mask clouds"), zlim=c(-1,1))
    plot(mask_all,main=c("Mask all"), zlim=c(-1,1))
    #plot(mask_focal,main=c("Mask focal"))
    dev.off()
  }

  #----------------------------
  # (D). Plot MLHC classification 
  #----------------------------
  breakpoints <- c(-0.1,0.9,1.9,2.9,3.9,4.9,5.9,6.9)
  colors <- c("black","yellow","dark green","light green","red","dark gray")
  if (plot) {
    png(filename = file.path(fig_dir,paste("mlhc_",img_id,run_id,".png",sep="")),width = 1380, height = 680, 
        units = "px", pointsize = 12,bg = "white")
    par(mfrow = c(1, 2),cex = 1,mar = c(0, 2, 0, 0), oma = c(1.8, 2.5, 1.5, 1.2),tcl = -0.25,mgp = c(2, 0.6, 0))
    #plot RGB
    raster::plotRGB(img,r=3,g=2,b=1,scale=9000,stretch="hist",ext=ext_shp,xaxt = "n",yaxt = "n")
    plot(extent(tmp_img),add=T,lwd=2)
    plot(shp,add=T,col=c("blue"))
    mtext("RGB image",side=3,cex=1)
    # plot classification
    plot(SC.mlhc$map,breaks=breakpoints,col=colors,legend=F,xaxt = "n",bty="n",yaxt = "n",xaxt = "n",yaxt = "n",box=FALSE)
    plot(SC.mlhc$map,breaks=breakpoints,col=colors,legend.only=TRUE,legend.shrink=0.75,
         axis.args=list(at=seq(0.5, 5.5, 1),
                      labels=c("NoData","GT","RT","HB","BG","SH"),
                      cex.axis=1.6))
    mtext("Maximum likelihood classification",side=3,cex=1)
    plot(extent(tmp_img),add=T,lwd=2)
    dev.off()
  }

  #----------------------------
  # (E). Plot RF classification
  #----------------------------
  if (plot) {
    png(filename = file.path(fig_dir,paste("rf_",img_id,run_id,".png",sep="")),width = 1380, height = 680, 
      units = "px", pointsize = 12,bg = "white")
    par(mfrow = c(1, 2),cex = 1,mar = c(0, 2, 0, 0), oma = c(1.8, 2.5, 1.5, 1.2),tcl = -0.25,mgp = c(2, 0.6, 0))
    #plot RGB
    raster::plotRGB(img,r=3,g=2,b=1,scale=9000,stretch="hist",ext=ext_shp,xaxt = "n",yaxt = "n")
    plot(extent(tmp_img),add=T,lwd=2)
    plot(shp,add=T,col=c("blue"))
    mtext("RGB image",side=3,cex=1)
    # plot classification
    plot(SC.rf$map,breaks=breakpoints,col=colors,legend=F,xaxt = "n",bty="n",yaxt = "n",xaxt = "n",yaxt = "n",box=FALSE)
    plot(SC.rf$map,breaks=breakpoints,col=colors,legend.only=TRUE,legend.shrink=0.75,
         axis.args=list(at=seq(0.5, 5.5, 1),
                      labels=c("NoData","GT","RT","HB","BG","SH"),
                      cex.axis=1.6))
    mtext("RandomForest classification",side=3,cex=1)
    plot(raster::extent(tmp_img),add=T,lwd=2)
    dev.off()
  }

  #----------------------------
  # (F). Plot Vector Support Machine Classification
  #----------------------------
  if (plot) {
    png(filename = file.path(fig_dir,paste("svm_",img_id,run_id,".png",sep="")),width = 1380, height = 680, 
        units = "px", pointsize = 12,bg = "white")
    par(mfrow = c(1, 2),cex = 1,mar = c(0, 2, 0, 0), oma = c(1.8, 2.5, 1.5, 1.2),tcl = -0.25,mgp = c(2, 0.6, 0))
    #plot RGB
    raster::plotRGB(img,r=3,g=2,b=1,scale=9000,stretch="hist",ext=ext_shp,xaxt = "n",yaxt = "n")
    plot(extent(tmp_img),add=T,lwd=2)
    plot(shp,add=T,col=c("blue"))
    mtext("RGB image",side=3,cex=1)
    # plot classification
    plot(SC.svm$map,breaks=breakpoints,col=colors,legend=F,xaxt = "n",bty="n",yaxt = "n",xaxt = "n",yaxt = "n",box=FALSE)
    plot(SC.svm$map,breaks=breakpoints,col=colors,legend.only=TRUE,legend.shrink=0.75,
         axis.args=list(at=seq(0.5, 5.5, 1),
                      labels=c("NoData","GT","RT","HB","BG","SH"),
                      cex.axis=1.6))
    mtext("RandomForest classification",side=3,cex=1)
    plot(extent(tmp_img),add=T,lwd=2)
    dev.off()
  }
  message("---Done Plotting----")

  #------------------------------------------------------------------------
  # Write rasters after mask implementation
  #------------------------------------------------------------------------
  #- Saving training data
  save(eval.arr,train.arr,trainData,evalData,cor.train.mtx,fit_msel,index_list_arr, sort.vardat,include2,
       accuracy.mlhc,error.mtx.mlhc,accuracy.rf,error.mtx.rf,accuracy.svm,error.mtx.svm,
       file= file.path(outdir, paste(img_id,run_id,".Rdata",sep="")))
  # saving images
  raster::writeRaster(SC.mlhc$map,file.path(outdir, paste("class_mlhc_",img_id,run_id,".tif",sep='')),overwrite=TRUE)
  raster::writeRaster(SC.rf$map,  file.path(outdir, paste("class_rf_",img_id,run_id,".tif",sep='')),overwrite=TRUE)
  raster::writeRaster(SC.svm$map, file.path(outdir, paste("class_svm_",img_id,run_id,".tif",sep='')),overwrite=TRUE)

  ##################################################################################################
  # Write out text file with results
  ##################################################################################################
  ptm.end <- Sys.time()
  message(paste(Sys.time()))
  time.it.took = print(ptm.end-ptm)

  #https://urldefense.com/v3/__http://www.cookbook-r.com/Data_input_and_output/Writing_text_and_output_from_analyses_to_a_file/__;!!JYXjzlvb!hOLHnMhmrMM9FiGDneDIoKNWZSgAdlCDbAFaml3Ql0HdLPd4ooC8gHNGqWrubAdfPL4IR6bNULMH9K6525NzkgkZ8XCNNX0EIjI$ 
  out_info_file = file.path(outdir, paste(img_id,run_id,"_test.txt",sep=''))
  sink(out_info_file)
  cat(paste("##################################################################################################\n"))
  cat(paste("###################### Start of code ######################\n"))
  cat(paste("##################################################################################################\n"))
  cat(paste("-> Image ID:",img_id,"\n"))
  cat(paste("Start time:",ptm,"\n"))
  cat(paste("End time:  ",ptm.end,"\n"))
  cat(paste("Run time:  ",time.it.took[],"\n"))
  cat(paste("#--------------------------------------#"),"\n")
  cat(paste("-> Run ID:                   ",run_id),"\n")
  cat(paste("-> Shapefile ID (eval):      ",shp_id),"\n")
  cat(paste("-> Seed:                     ",seed),"\n")
  cat(paste("-> Fig dir:                  ",fig_dir),"\n")
  cat(paste("-> Out dir:                  ",outdir),"\n")
  cat(paste("##################################################################################################\n"))
  cat(paste("########## Parameters ##########"),"\n")
  cat(paste("-> Cloud theshold (data<Band2):         ",cld_trhd),"\n")
  cat(paste("-> Cloud shadow theshold (data>Band3):  ",cld_shadow_trhd),"\n")
  cat(paste("-> Number of indices:                   ",no.bands.class),"\n")
  cat(paste("########## Image info ##########"),"\n")
  cat(paste("Extent: ",ext_shp),"\n")
  print("All Image / NA / Mask / Data")
  print(p.area)
  print(paste("Total Area (Km2):",area.img))
  print(paste("Total NA   (Km2):",area.na))
  print(paste("Total Mask (Km2):",area.mask))
  print(paste("Total Data (Km2):",area.data))
  cat(paste("##################################################################################################\n"))
  cat(paste("########## (1) Initial RF results ##########"),"\n")
  cat(print(rf.classifier))
  #get_value(pred.rf.maj$layer,class_list)
  print(class_txt_arr)
  cat(print(accuracy(eval.arr[,n.bands+1],pred.arr)),"\n")
  print("---correlation matrix train arr---")
  print(cor.train.mtx)
  print("---Variable selection----")
  print(fit_msel)
  print(paste(index_list_arr))
  cat(paste("##################################################################################################\n"))
  print("########## (2) MLHC results ##########")
  print("class percentages")
  print(p.class.mlhc)
  print("Mask(0) # Bare Ground (1) # Green Trees (2) # Herbaceous (3) # Red trees (4) # Shadow (5) # Gray trees (6)  # EXtra class (7) # No data (8)")
  print("-------Confusion Matrices--------")
  print(accuracy.mlhc)
  print(error.mtx.mlhc)
  print("----Olof cfm----")
  print(olof.table.mlhc)
  cat(paste("##################################################################################################\n"))
  print("########## (3) RandomForest results ##########")
  print("class percentages")
  print(p.class.rf)
  print("Mask(0) # Bare Ground (1) # Green Trees (2) # Herbaceous (3) # Red trees (4) # Shadow (5) # Gray trees (6)  # EXtra class (7) # No data (8)")
  print("-------Confusion Matrices--------")
  print(accuracy.rf)
  print(error.mtx.rf)
  print("----Olof cfm----")
  print(olof.table.rf)
  cat(paste("##################################################################################################\n"))
  print("########## (4) SVM results ##########")
  print("class percentages")
  print(p.class.svm)
  print("Mask(0) # Bare Ground (1) # Green Trees (2) # Herbaceous (3) # Red trees (4) # Shadow (5) # Gray trees (6)  # EXtra class (7) # No data (8)")
  print("-------Confusion Matrices--------")
  print(accuracy.svm)
  print(error.mtx.svm)
  print("----Olof cfm----")
  print(olof.table.svm)
  cat(paste("##################################################################################################\n"))
  cat(paste("###################### End of code ######################"),"\n")
  cat(paste("##################################################################################################\n"))
  #sink()
  closeAllConnections()   # Close Sink file


  #--------------------------------------------------------------------------------------------
  # Write out Confusion matrices to CSV
  # -- Move to FUNCTION
  #--------------------------------------------------------------------------------------------
  # Error matrix -- MLHC (not maj)
  PA = error.mtx.mlhc$byClass[,3] # Pos Pred Value >> not sure if this is correct?
  UA = error.mtx.mlhc$byClass[,4] # Neg Pred Value >> not sure if this is correct?
  mat = error.mtx.mlhc$table
  order.txt = c("gt","rt","yt","hb","bg","sh")
  order = match(order.txt,colnames(mat))
  mat = mat[c(order),c(order)]
  mat = rbind(mat,colSums(mat),UA[order[1:6]])
  mat = cbind(mat,rowSums(mat),c(PA[order[1:6]],NA,NA))
  rownames(mat) = c(order.txt,"sum","UA")
  colnames(mat) = c(order.txt,"sum","PA")
  mat
  info.mat = paste("error.mtx.mlhc$overall"," Overall Acc:",round(as.numeric(error.mtx.mlhc$overall[1]*100),3),"%")
  write.table(info.mat, file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = FALSE, sep = ',', col.names = F, row.names = F)
  write.table(mat, file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = T, row.names = TRUE)
  write.table(c("-------"), file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)
  #- - - - - - - - - - - - - - - - - 
  # Ollofson Error matrix -- MLHC (not maj)
  mat = olof.table.mlhc$matrix
  order.txt = c("gt","rt","yt","hb","bg","sh","sum")
  order = match(order.txt,colnames(mat))
  mat = mat[c(order),c(order)]
  mat = rbind(mat,c(olof.table.mlhc$PA[order[1:6]],NA))
  mat = cbind(mat,c(olof.table.mlhc$UA[order[1:6]],NA,NA))
  rownames(mat) = c(order.txt,"UA")
  colnames(mat) = c(order.txt,"PA")
  mat
  info.mat = paste("olof.table.mlhc$matrix, Overall Acc:",round(olof.table.mlhc$OA*100,3),"%")
  write.table(info.mat, file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = F, row.names = F)
  write.table(mat, file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = c(order.txt,"PA"), row.names = TRUE)
  write.table(c("-------"), file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)
  #--------------------------------------------------------------------------------------------
  # Error matrix -- rf (not maj)
  PA = error.mtx.rf$byClass[,3] # Pos Pred Value >> not sure if this is correct?
  UA = error.mtx.rf$byClass[,4] # Neg Pred Value >> not sure if this is correct?
  mat = error.mtx.rf$table
  order.txt = c("gt","rt","yt","hb","bg","sh")
  order = match(order.txt,colnames(mat))
  mat = mat[c(order),c(order)]
  mat = rbind(mat,colSums(mat),UA[order[1:6]])
  mat = cbind(mat,rowSums(mat),c(PA[order[1:6]],NA,NA))
  rownames(mat) = c(order.txt,"sum","UA")
  colnames(mat) = c(order.txt,"sum","PA")
  mat
  info.mat = paste("error.mtx.rf$overall"," Overall Acc:",round(as.numeric(error.mtx.rf$overall[1]*100),3),"%")
  write.table(info.mat, file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = T, sep = ',', col.names = F, row.names = F)
  write.table(mat, file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = T, row.names = TRUE)
  write.table(c("-------"), file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)
  #- - - - - - - - - - - - - - - - - 
  # Ollofson Error matrix -- rf (not maj)
  mat = olof.table.rf$matrix
  order.txt = c("gt","rt","yt","hb","bg","sh","sum")
  order = match(order.txt,colnames(mat))
  mat = mat[c(order),c(order)]
  mat = rbind(mat,c(olof.table.rf$PA[order[1:6]],NA))
  mat = cbind(mat,c(olof.table.rf$UA[order[1:6]],NA,NA))
  rownames(mat) = c(order.txt,"UA")
  colnames(mat) = c(order.txt,"PA")
  mat
  info.mat = paste("olof.table.rf$matrix, Overall Acc:",round(olof.table.rf$OA*100,3),"%")
  write.table(info.mat, file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = F, row.names = F)
  write.table(mat, file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = c(order.txt,"PA"), row.names = TRUE)
  write.table(c("-------"), file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)
  #--------------------------------------------------------------------------------------------
  # Error matrix -- svm (not maj)
  PA = error.mtx.svm$byClass[,3] # Pos Pred Value >> not sure if this is correct? 
  UA = error.mtx.svm$byClass[,4] # Neg Pred Value >> not sure if this is correct?
  mat = error.mtx.svm$table
  order.txt = c("gt","rt","yt","hb","bg","sh")
  order = match(order.txt,colnames(mat))
  mat = mat[c(order),c(order)]
  mat = rbind(mat,colSums(mat),UA[order[1:6]])
  mat = cbind(mat,rowSums(mat),c(PA[order[1:6]],NA,NA))
  rownames(mat) = c(order.txt,"sum","UA")
  colnames(mat) = c(order.txt,"sum","PA")
  mat
  info.mat = paste("error.mtx.svm$overall"," Overall Acc:",round(as.numeric(error.mtx.svm$overall[1]*100),3),"%")
  write.table(info.mat, file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = T, sep = ',', col.names = F, row.names = F)
  write.table(mat, file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = T, row.names = TRUE)
  write.table(c("-------"), file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)
  #- - - - - - - - - - - - - - - - - 
  # Ollofson Error matrix -- svm (not maj)
  mat = olof.table.svm$matrix
  order.txt = c("gt","rt","yt","hb","bg","sh","sum")
  order = match(order.txt,colnames(mat))
  mat = mat[c(order),c(order)]
  mat = rbind(mat,c(olof.table.svm$PA[order[1:6]],NA))
  mat = cbind(mat,c(olof.table.svm$UA[order[1:6]],NA,NA))
  rownames(mat) = c(order.txt,"UA")
  colnames(mat) = c(order.txt,"PA")
  mat
  info.mat = paste("olof.table.svm$matrix, Overall Acc:",round(olof.table.svm$OA*100,3),"%")
  write.table(info.mat, file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = F, row.names = F)
  write.table(mat, file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = c(order.txt,"PA"), row.names = TRUE)
  write.table(c("-------"), file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)
  #--wRITE VARIABLE SELECTION----------------------
  write.table(paste("vARIABLE SELECTION (mir):",varsel.treshold,"Tresholds"), file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)
  write.table(print(sort(vardat,decreasing = T)), file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)
  write.table(c("-------"), file.path(outdir, paste(img_id,"_cfm",run_id,".csv",sep='')), append = TRUE, sep = ',', col.names = TRUE, row.names = TRUE)

  #--------------------------------------------------------------------------------------------
  # Write out metrics to CSV
  #--------------------------------------------------------------------------------------------
  header = c("img","Total_area(km2)","Total_NA","Total_mask","Total_data",
             "OA_mlhc","OA_mlhc_olof","OA_rf","OA_rf_olof","OA_svm","OA_svm_olof",
             names(n.class.rf.arr),names(n.class.rf.arr),names(n.class.rf.arr))
  OA_metrics = c(img_id,area.img,area.na,area.mask,area.data,
                    as.numeric(error.mtx.mlhc$overall[1]),olof.table.mlhc$OA,
                    as.numeric(error.mtx.rf$overall[1]),olof.table.rf$OA,
                    as.numeric(error.mtx.svm$overall[1]),olof.table.svm$OA,
                    as.numeric(c(n.class.mlhc.arr[],n.class.rf.arr[],n.class.svm.arr[])))
  write.table(rbind(header,OA_metrics), file.path(outdir, paste(img_id,"_metrics",run_id,".csv",sep='')), append = FALSE, sep = ',', col.names = F, row.names = F)

  message(paste("End running.......",img_id,"...Number:",imgno))
} # End for IMG LOOP


#------------------------------------------------------------------------
#- End 
#------------------------------------------------------------------------

