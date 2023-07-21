# Code to compute river planform change detection and extract relevant statistics in GEE
# Anya Leenman
# 2022
#
# Note: requires you to install rgee/GEE first. Follow instructions here: https://github.com/r-spatial/rgee
#
# ----
# housekeeping
#----
rm(list = ls())

library(reticulate) # to manage python envs 
use_condaenv('rgee', required = TRUE) # set correct python env

# load packages
library(rgee)
library(dplyr)
library(sf)
library(cptcity)
library(lubridate)

ee_Initialize() # initialize GEE

# edit as needed
gettinginfo = F # get statistics calculated over region? Note this slows things down; 
# set to false if you just want to make plots
make_plots = T # make pretty interactive maps of change within the time/space filters you
# set below.

window_length <- 25 # n days to search for imagery (before + after flood)
start_date <- '2020-08-20' # flood start date
end_date <- '2020-12-16' # flood end date
# set vertices of polygon around your Area of Interest (AOI) (or import a polygon)
AOI <- st_polygon(list(matrix(c(-73.93,7.31243, 
                                -73.93,7.50282,
                                -73.885,7.50078,
                                -73.885,7.31243,
                                -73.93,7.31243),
                              ncol=2, 
                              byrow=TRUE)))

channel_belt = T # do you want to compute channel maps for channel belt (including bars)?
area_keep_threshold = 0.2 # minimum proportion of polygon that must have cloud-free data
snowthresh <- 0.6 # what proportion of cloud-free AOI must be snow covered before AOI is rejected? must be high as water often misclassified as snow
MAX_CLOUD_PROBABILITY <- 10 # max probability of cloud retained when masking clouds
# lower values = less clouds but some misclassification (and therefore masking out) of water/sediment as clouds
# higher values = less mis-masking of water, but more mis-classification of mid-river clouds as land (and therefore spurious change detection)
nMths <- 1 # window length to make mosaics for "after flood" (to check for permanent vs transient)
months_to_check <- 1:24 # how many times after the flood should be checked for whether transitions are permanent/transient?
# e.g. for '1', the first month to be checked STARTS 1 month after the end of the flood (And ends 2 mths after).

# classification parameters (from Boothroyd et al)
ndsi_param <- 0.5  # for snow detection - threshold of 0.4 used by Hofmeister et al 2022 JoHX
mndwi_param <- -0.40;
ndvi_param <- 0.20;
cleaning_pixels <- 10; # n connected pixels removed
run_cleaning <- T # should neighborhood-based noise filter be run? 

# don't edit below this line
#------------------------------------------------------------------------------
# timestamp cleaning:
pre_end <- start_date # end of pre-flood window
pre_srt <- as.character(ymd(pre_end) - days(window_length)) # start of pre-flood window
post_srt <- end_date # start of post-flood window
post_end <- as.character(ymd(post_srt) + days(window_length)) # end of post-flood window

# Earth engine datasets:
GEE_dataset <- 'COPERNICUS/S2_HARMONIZED'
cloud_dataset <- 'COPERNICUS/S2_CLOUD_PROBABILITY'
GRWL <- "projects/sat-io/open-datasets/GRWL/water_mask_v01_01"
bands_inc_cloud <- c('B2', 'B3', 'B4', 'B8', 'B8A', 'B9', 'B11', 'B12') # bands to use (for cloud masking, and DeepwaterMap2 if using)

# -----
# Def functions:
#------
# function to mask clouds out of sentinel data:
maskClouds <- function(img){
  clouds <- ee$Image(img$get('cloud_mask'))$select('probability')
  isNotCloud <- clouds$lt(MAX_CLOUD_PROBABILITY)
  return(img$updateMask(isNotCloud)$divide(10000)) # need to add /10000 scaling parameter
}

# // The masks for the 10m bands sometimes do not exclude bad data at
# // scene edges, so we apply masks from the 20m and 60m bands as well.
maskEdges <- function(s2_img){
  return(s2_img$updateMask(
    s2_img$select('B8A')$mask()$updateMask(s2_img$select('B9')$mask())))
}

# function to map water from chosen input bands:
mapWater <- function(img){
  # band ratios to use as inputs:
  ndvi_vals <- img$normalizedDifference(c("B8","B4"))
  # lswi_vals <- img$normalizedDifference(c("B8", "B11"))
  mndwi_vals <- img$normalizedDifference(c("B3", "B11"))
  evi_vals <- img$expression('2.5 * (Nir - Red) / (1 + Nir + 6 * Red - 7.5 * Blue)', 
                             c('Nir' = img$select('B8'), 
                               'Red' = img$select('B4'), 
                               'Blue'= img$select('B2')))
  # water (raw)
  water <- (mndwi_vals$gt(ndvi_vals)$
              Or(mndwi_vals$gt(evi_vals)))$
    And(evi_vals$lt(0.1)) # default 0.1 from Boothroyd et al 2020 code
  # channel belt (counts bars as 'river')
  if(channel_belt == T){
    activebelt <- (mndwi_vals$gte(mndwi_param))$And(ndvi_vals$lte(ndvi_param))
    water <- water$Or(activebelt)
  }
  return(water)
}

# function to detect snow (only apply if (!country %in% snowfree_countries) so 
# it runs as default unless country explicitly in snowfree list)
# following methods here: https://www.sciencedirect.com/science/article/pii/S2589915522000050
mapSnow <- function(img){
  # band ratios to use as inputs:
  ndsi_vals <- img$normalizedDifference(c("B3", "B11"))
  # snow
  snow <- ndsi_vals$gt(ndsi_param)
  return(snow)
}

# function to create the image collection
my_imcol <- function(startdate, enddate, aoi){
  imcol <- ee$  # create collection
    ImageCollection(GEE_dataset)$ # query Earth Engine database
    select(bands_inc_cloud)$ # filter to bands of interest
    filterDate(startdate, enddate)$ # filter to time window of interest
    filterBounds(aoi)$ # filter to aoi for this gauge
    map(maskEdges) # apply func to mask out bad data at edges
}

# function to mask out clouds, reduce to min for each pixel:
rem_clouds <- function(imcol, startdate, enddate, aoi){
  # get image collection of s2cloudless data:
  s2Clouds <- ee$
    ImageCollection(cloud_dataset)$ # query GEE servers
    filterDate(startdate, enddate)$ # filter to time window of interest
    filterBounds(aoi) # filter to aoi for this gauge
  # Join S2 SR with cloud probability dataset to add cloud mask *for each image* :
  s2SrWithCloudMask <- ee$
    Join$saveFirst('cloud_mask')$
    apply('primary' = imcol, 
          'secondary' = s2Clouds,
          'condition' = ee$Filter$equals(
            'leftField'= 'system:index', 
            'rightField'= 'system:index'))
  # mask clouds and reduce to min
  imcol <- ee$ImageCollection(s2SrWithCloudMask)$map(maskClouds)$min()
  # note "min" is a reducer that turns ImageCollection to Image and also finds minimum 
  # was using median but it misclassified cloud as water too often
  
  # note also that this method >> filtering image collection by % cloud as that would 
  # eliminate more cloudy scenes, even if part of a scene is not cloudy and might 
  # be useful.
}

#---------
# for each combo of a flood event and an aoi around the gauge that measured it:
#----------
AOI <- AOI %>% # get aoi polygon
  sf_as_ee() # convert to EE format for query

# get image collection covering whole AOI to get n. pixels in AOI: 
# NOTE: if running this for multiple floods at one gauge/AOI, calculate area first and then apply 
# code below to all the floods, rather than calculating area again for each flood!
imcol_mean <- ee$ImageCollection('COPERNICUS/S2_HARMONIZED')$ # query Earth Engine database
  filterDate('2020-07-01', '2021-01-01')$ # filter to 6 mth to ensure complete coverage (overkill)
  filterBounds(AOI)$ # filter to aoi for this gauge
  select('B2')$ # only one band
  mean() # force imcol to one image
A_poly <- imcol_mean$reduceRegion( # calc n cells (incl cloud) in AOI
  reducer = ee$Reducer$count(), 
  geometry = AOI, 
  scale = 10,
  maxPixels = 1e30)
A_poly <- A_poly$get("B2")$getInfo() # extract count of cells within AOI
area_m2 <- A_poly * 10 * 10 # conv cells to m

# image collections for mapping
pre_flood2 <- my_imcol(pre_srt, pre_end, AOI)
post_flood2 <- my_imcol(post_srt, post_end, AOI)

# get number of images in ImageCollections:
(n_images <- pre_flood2$size()$getInfo())
(n_images_post <- post_flood2$size()$getInfo()) 

# if either image collection is empty:
if((n_images < 1) | (n_images_post < 1)){
  print(paste0('At i = ', i, ' no images in date range for at least one of pre or post'))
  next # move on to next flood event
}

# remove clouds and reduce to min for each pixel:
pre_flood2 <- rem_clouds(pre_flood2, pre_srt, pre_end, AOI)
post_flood2 <- rem_clouds(post_flood2, post_srt, post_end, AOI)

# copy cloud mask from one to the other so both have matching cloud mask:
pre_flood2 <- pre_flood2$updateMask(post_flood2$mask())
post_flood2 <- post_flood2$updateMask(pre_flood2$mask())

# get total polygon area in number of pixels, but ACCOUNT FOR CLOUDS, 
# i.e. need to extract total number of non-NULL cells within AOI:
A <- pre_flood2$reduceRegion(
  reducer = ee$Reducer$count(), # gets total num of non-NULL cells
  geometry = AOI, 
  scale = 10,
  maxPixels = 1e30
)

A <- A$get("B2")$getInfo() # extract count of non-null cells within AOI
dat <- data.frame(cloudfree_area_px = A) #  assign to output df
A_raw <- A * 10 * 10 # convert to area in m2 by multiplying by cell size
dat$cloudfree_area_m2 <- A_raw # assign total area to data frame
if(A_raw < (area_keep_threshold * area_m2)){
  print('Not enough cloudfree data')
  next
}# code below skipped if cloudfree area < some threshold of total AOI area

# test for snow cover:
snow_pre <- mapSnow(pre_flood2)
snow_pre_count <- snow_pre$reduceRegion(
  reducer = ee$Reducer$sum(), # gets total num of snow cells
  geometry = AOI, 
  scale = 10,
  maxPixels = 1e30
)$get('nd')$getInfo() # extract count of non-null cells within AOI
if(snow_pre_count / A > snowthresh ){
  print('SNOW ALERT at preflood')
  next
}
snow_post <- mapSnow(post_flood2)
snow_post_count <- snow_post$reduceRegion(
  reducer = ee$Reducer$sum(), # gets total num of non-NULL cells
  geometry = AOI, 
  scale = 10,
  maxPixels = 1e30
)$get('nd')$getInfo() # extract count of non-null cells within AOI
if(snow_post_count / A > snowthresh ){
  print('SNOW ALERT at postflood')
  next
}

dat$cloudfree_percentage <- (A_raw / area_m2) * 100 # save cloudfree %

# Water segmentation/extraction/mapping (from Zou et al 2018; 
# code adapted from Boothroyd et al 2020 WIRES water paper):
#--------
water_pre <- mapWater(pre_flood2) # apply function defined at beginning
water_post <- mapWater(post_flood2)

#------------------------------------------------------
# generate after-flood water mask using months specified (for filtering out transient change):
#------------------------------------------------------
water_after_sum <- ee$Image(0) # empty layers
cloudfree_count <- ee$Image(1) # start with empty layer of 1s because we divide is_wetting_permanent by the number of afterflood timestamps + 1 
# loop through after-flood months:
for(j in 1:length(months_to_check)){
  after_flood <- my_imcol(startdate = as.character(as.Date(post_end) %m+% months(j)), 
                          enddate = as.character(as.Date(post_end) %m+% months(j + nMths)),
                          aoi = AOI)
  after_flood <- rem_clouds(imcol = after_flood, startdate = as.character(as.Date(post_end) %m+% months(j)), 
                            enddate = as.character(as.Date(post_end) %m+% months(j + nMths)),
                            aoi = AOI)
  # error trap:
  if(after_flood$bandNames()$getInfo()[1] != 'B2'){
    print('No data in this after_flood imcol; breaking out of loop')
    break
  }
  water_after <- mapWater(after_flood) # map water
  # copy of the cloud mask
  cloudfree_after <- ee$Image(1)$ # blank image of 1s
    updateMask(water_after$mask())$ # give it the cloud mask of the water mask layer
    unmask() # zeros in cloud areas
  # touch up the cloud mask for the water map:
  water_after <- water_after$unmask()$ # apply zeros in cloud holes
    updateMask(water_post$mask()) # apply cloud mask from flood layers so all data use same mask
  water_after_sum <- water_after$add(water_after_sum) # add to total after-flood stack water-presence stack
  cloudfree_count <- cloudfree_after$add(cloudfree_count) # add to total after-flood cloud stack
}
if(j <12 ){
  print('not enough post-flood data')
  next
}
dat$n_afterflood_mths <- j

# -----------------------------------------------------------------------
# compute change detection for channel mask and water mask
# -----------------------------------------------------------------------

# compute difference:
D_water_map_raw <- water_post$subtract(water_pre)
D_water_map <- D_water_map_raw$
  abs() # take absolute vals
# so that 0 = no change (either land, or water); 1 = change (either to or from water)

# should cleaning-pixels be run using code above?
if(run_cleaning == F){
  water_pre_corrected <- water_pre
  water_post_corrected <- water_post
}
if(run_cleaning == T){
  #--------------------
  # noise filtering:
  #--------------------
  # because of cloud holes, running the 'cleaning pixels' GEE function removes legitimate
  # patches of water unless the pre/post-flood water masks are first amalgamated with 
  # the spatially continuous GRWL mask. This step does so, so that 'cleaning pixels' only 
  # removes patches of 'water' that are known to be outside the channel belt.
  
  # GRWL:
  GRWL_mask <- ee$ 
    ImageCollection(GRWL)$ # query Earth Engine database
    filterBounds(AOI)$ # filter to aoi for this gauge 
    mean()$ # imcol to image
    divide(255) # rescale to binary (But this converts data type to float...)
  GRWL_unmasked <- GRWL_mask$unmask() # add zeros in non-water areas
  
  # mask using a giant water map (only run cleaningpixels code on 1 layer, not 2)
  # same pixels cleaned from both layers 
  giantwatermask <- water_pre$add(water_post)$
    add(GRWL_unmasked)$
    unmask(GRWL_unmasked)
  giantwatermask <- giantwatermask$
    where(giantwatermask$gt(0),1)$ # reclassify all to 1
    updateMask(giantwatermask$gt(0))$ # select only water
    uint8() # convert back to integer
  giantwatermask <- giantwatermask$
    updateMask(giantwatermask$ 
                 connectedPixelCount(cleaning_pixels, F)$ # calculate size of clumps
                 gte(cleaning_pixels)) # remove clumps smaller than 'cleaning_pixels '
  
  
  D_water_map2 <- D_water_map$updateMask(giantwatermask$mask())$ # mask out lakes/noise/crap
    unmask()$ # re-add zeros in all holes created by above mask
    updateMask(water_pre$mask()) # add the cloud mask back in
  cleaned_cells <- D_water_map2$subtract(D_water_map)$
    add(1)
  D_water_map <- D_water_map2 # rename 
  # correct water masks so that cells removed from change map are marked 'dry' in both water masks:
  water_pre_corrected <- water_pre$multiply(cleaned_cells)
  water_post_corrected <- water_post$multiply(cleaned_cells)
}

# mask so that only cells showing change are left ( reduces risk of hitting maxPixels, 
# even though the count() reducer is slightly slower than sum()).
D_water_map_change <- D_water_map$updateMask(D_water_map$gt(0))

# Update mask for raw change map with only changed pixels:
D_water_map_raw <- D_water_map_raw$updateMask(D_water_map_change$mask())

# then, for cleaned raw change map:
# if pixel = 1 (i.e. became wet):
is_wetting_permanent <- D_water_map_raw$updateMask(D_water_map_raw$gt(0))$ # mask to cells that became wet
  #   is pixel wet in afterflood data? add to check: 
  add(water_after_sum)$
  divide(cloudfree_count)
#     if pixel remained wet, e.g. (1+1+1+1...1) / n = n/n so yes, 
#     if pixel became dry again e.g. (1+0+0+0...0) / n = 1/n so no
#     if pixel became dry again but THEN got wet again a lot later? (1+1+1+0+1...1) / n = (n-1)/n so not permanent
# mask to ONLY cells that became permanently wet
permanently_wetted <- is_wetting_permanent$updateMask(is_wetting_permanent$eq(1))
# make layer of cells that are transiently wetted so these can be marked as 'dry' in the corrected post-flood water map:
# convert to vector for plotting
permanently_wetted_int <- permanently_wetted$uint8()
permanently_wetted_vector <- permanently_wetted_int$reduceToVectors(
  scale= 10,
  geometry= AOI,
  geometryType= 'polygon',
  maxPixels= 1e8,
  eightConnected = T,
  labelProperty = 'zone')
transiently_wetted_step <- is_wetting_permanent$
  updateMask(is_wetting_permanent$lt(1))
transiently_wetted <- transiently_wetted_step$
  where(transiently_wetted_step$neq(0), 0)$ # reclassify so transiently wetted = 0 
  unmask(1) # and everything else = 1. Then multiply this by post-flood water mask to set transiently 'wet' cells to 'dry'
water_post_corrected <- water_post_corrected$multiply(transiently_wetted)

# if pixel = -1 (i.e. became dry due to channel abandonment or transient stage reduction):
is_drying_permanent <- D_water_map_raw$updateMask(D_water_map_raw$lt(0))$ # mask to cells that became dry
  add(water_after_sum)
#   is pixel dry in all afterflood data? multiply to check: if pixel remained dry it is permanently abandoned, so 
#     -1 + 0 + 0 + 0 + 0 = -1
#   is pixel became wet again at any point (e.g. re-inundated due to stage change) it is only transiently abandoned 
# e.g. due to stage reduction during flood, so
#     -1 + 0 + 0 + 0 + 1 = 0 (or higher)
permanently_dried <- is_drying_permanent$updateMask(is_drying_permanent$eq(-1))
# make layer of cells that are transiently dried so these can be marked as 'dry' in the corrected pre-flood water map:
transiently_dried_step <- is_drying_permanent$updateMask(is_drying_permanent$gt(-1))
transiently_dried <- transiently_dried_step$
  where(transiently_dried_step$neq(-1),0)$ # reclassify so transiently dried = 0
  unmask(1) # and everything else = 1. Then multiply this by pre-flood water mask to set transiently 'dried' cells to 'dry'
water_pre_corrected <- water_pre_corrected$multiply(transiently_dried)

if(gettinginfo == T){ # getting data for each flood:
  # water surface area before flood (as metric for 'river size'):
  water_area_pre <- water_pre_corrected$reduceRegion(
    reducer = ee$Reducer$sum(), # adds all ones, i.e. water cells ()
    geometry = AOI,
    scale = 10,
    maxPixels = 1e30
  )
  dat$water_area_preflood_in_px <- water_area_pre$get('nd')$getInfo()
  
  # water surface area after flood
  water_area_post <- water_post_corrected$reduceRegion(
    reducer = ee$Reducer$sum(), # adds all ones, i.e. water cells ()
    geometry = AOI,
    scale = 10,
    maxPixels = 1e30
  )
  dat$water_area_postflood_in_px <- water_area_post$get('nd')$getInfo()
  
  # average water surface area:
  dat$water_area_av_in_px <- mean(c(dat$water_area_preflood_in_px, 
                                    dat$water_area_av_in_px), 
                                  na.rm = T)
  
  # sum permanently dried/wetted cells + save to dataframe
  dried_count <- permanently_dried$reduceRegion(
    reducer = ee$Reducer$count(),
    geometry = AOI,
    scale = 10,
    maxPixels = 1e30
  )
  dat$perm_dried_px <- dried_count$get('nd')$getInfo()
  
  wetted_count <- permanently_wetted$reduceRegion(
    reducer = ee$Reducer$count(),
    geometry = AOI,
    scale = 10,
    maxPixels = 1e30
  )
  dat$perm_wet_px <- wetted_count$get('nd')$getInfo()
  # normalise by count of cloudfree pre-flood water surface area pixels for metric of change
  dat$perm_wetting_norm_by_WSA <- dat$perm_wet_px / dat$water_area_preflood_in_px 
}

print('Calculations done, mapping started')

#-------------
# Interactive map visualisations w/ leaflet:
#-------------
if(make_plots == T){
  # visualising change outlines:
  # Create an empty image into which to paint the features, cast to byte.
  empty <- ee$Image()$byte()
  
  # Paint all the polygon edges with the same number and 'width', display.
  outline_wetted <- empty$paint(
    featureCollection = permanently_wetted_vector,
    color = 1,
    width = 2
  )
  outline_AOI <- empty$paint(
    featureCollection = AOI,
    color = 1,
    width = 2
  )
  
  Map$centerObject(AOI)
  m <- Map$addLayer(
      eeObject = pre_flood2, # preflood mosaic
      visParams = list(
        bands = c("B4", "B3", "B2"),
        max = 0.2),
      shown = T,
      name = "rgb_pre"
    ) +
    Map$addLayer(
      eeObject = post_flood2, # postflood mosaic
      visParams = list(
        bands = c("B4", "B3", "B2"),
        max = 0.2),
      shown = F,
      name = "rgb_post"
    ) +
    Map$addLayer(
      eeObject = water_pre_corrected,
      visParams = list(min = 0, max = 1,
                       palette = c('000000', 'FFFF00')),
      shown = F,
      name = "water_mask_pre"
    ) +
    Map$addLayer(
      eeObject = water_post_corrected,
      visParams = list(min = 0, max = 1,
                       palette = c('000000', '0000FF')),
      shown = F,
      name = "water_mask_post"
    ) +
    Map$addLayer(permanently_wetted, visParams = list(palette = c('FFFF00')),
                 shown = F,
                 name = 'permanently eroded, raster') +
    Map$addLayer(transiently_wetted_step, visParams = list(palette = c('FFB6C1')),
                 shown = F,
                 name = 'transiently wetted (not eroded), raster') +
    Map$addLayer(outline_wetted, 
                 list(palette = "FFFF00"), 
                 "Permanently eroded, vector") +
    Map$addLayer(outline_AOI, 
                 list(palette = "FFFFFF"), 
                 "AOI")
  
  print(m)
}
