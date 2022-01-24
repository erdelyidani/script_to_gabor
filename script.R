all(sapply( c("parallel","doParallel","foreach","caret","pipeR","tidyr","glmnet","randomForestSRC","raster","rasterVis","rworldxtra","magrittr","sp","mapview", "elevatr", "sf", "readxl", "GSIF", "rnaturalearth","rnaturalearthhires","rgeos","rgdal","geosphere"),
            require, character.only = T))

path = c("C:/Users/Daniel/Google Drive/Erdélyi Dániel/Doktori/Alaa")
path2 = c("C:/Users/Daniel/Google Drive/Erdélyi Dániel/Doktori/Alaa/Map_KG-Global")

setwd(path)


#Create the Grid cell
##################################################
#ori <- SpatialPoints(cbind(-1000000, 3500000), proj4string =  CRS("+init=epsg:3857"))

#Small version:
ori <- SpatialPoints(cbind(-835000, 3710000), proj4string =  CRS("+init=epsg:3857"))

#version istvan:
#ori <- SpatialPoints(cbind(-960000, 3500000), proj4string =  CRS("+init=epsg:3857"))

##################################################
# Convert the projection of ori
# Use EPSG: 3857 (Spherical Mercator)
ori_t <- spTransform(ori, CRSobj = CRS("+init=epsg:3857"))
# The origin has been rounded to the nearest 100

x_ori <- round(coordinates(ori_t)[1, 1]/100) * 100
y_ori <- round(coordinates(ori_t)[1, 2]/100) * 100

# Define how many cells for x and y axis

x_cell <- 484
y_cell <- 218

##################################################

# Define the resolution to be 10000 meters
cell_size <- 10000
# Create the extent
ext <- extent(x_ori, x_ori + (x_cell * cell_size), y_ori, y_ori + (y_cell * cell_size))
# Initialize a raster layer
ras <- raster(ext)
# Set the resolution to be
res(ras) <- c(cell_size, cell_size)
# fill up de ras layer with a random value, without it, the code will not run.
# ras[] <- rnorm(length(ras))
ras[] <- 1
# Project the raster
projection(ras) <- CRS("+init=epsg:3857")
# Save the raster layer 
# Convert to spatial pixel
GRID <- rasterToPoints(ras, spatial = TRUE)
gridded(GRID) <- TRUE
GRID <- as(GRID, "SpatialPixelsDataFrame")

GRID_coordinates = data.frame(coordinates(GRID))

#some correction
GRID_coordinates$y <- GRID_coordinates$y - 5000
GRID_coordinates$d <- 1

GRID <- GRID_coordinates
coordinates(GRID) <- ~ x + y
proj4string(GRID) <- CRS("+init=epsg:3857")
gridded(GRID) <- TRUE
GRID <- as(GRID, "SpatialPixelsDataFrame")

GRID_coordinates = data.frame(coordinates(GRID))


#get elevation data to grid cells (This predictor wont be used, so you can skip this part if you want):
rasDEM_europe2 = get_elev_raster(GRID_coordinates, z=7, prj ="+init=epsg:3857")
rasDEM_europe2[rasDEM_europe2 < 0] <- 0
mapview(rasDEM_europe2)

GRID$elevation = extract(rasDEM_europe2, GRID)
mapview(GRID)


#Get KG predictor values
################################################################################################################################
#KG codes:
###########################################################################################
##
## R source code to read and visualize KÃ¶ppen-Geiger fields (Version of 27 December 2019)                                                                                    
##
## Climate classification after Kottek et al. (2006), downscaling after Rubel et al. (2017)
##
## Kottek, M., J. Grieser, C. Beck, B. Rudolf, and F. Rubel, 2006: World Map of the  
## KÃ¶ppen-Geiger climate classification updated. Meteorol. Z., 15, 259-263.
##
## Rubel, F., K. Brugger, K. Haslinger, and I. Auer, 2017: The climate of the 
## European Alps: Shift of very high resolution KÃ¶ppen-Geiger climate zones 1800-2100. 
## Meteorol. Z., DOI 10.1127/metz/2016/0816.
##
## (C) Climate Change & Infectious Diseases Group, Institute for Veterinary Public Health
##     Vetmeduni Vienna, Austria
##
###########################################################################################
 data(countriesHigh)

setwd(path2)

# Read raster files
period='1986-2010'
r <- raster(paste('KG_', period, '.grd', sep=''))

# Color palette for climate classification
climate.colors=c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC", "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64", "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00", "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF", "#64FFFF", "#F5FFFF")

# Legend must correspond to all climate classes, insert placeholders
r0 <- r[1:32]; r[1:32] <- seq(1,32,1)

# Converts raster field to categorical data
r <- ratify(r); rat <- levels(r)[[1]]

# Legend is always drawn in alphabetic order
rat$climate <- c('Af', 'Am', 'As', 'Aw', 'BSh', 'BSk', 'BWh', 'BWk', 'Cfa', 'Cfb','Cfc', 'Csa', 'Csb', 'Csc', 'Cwa','Cwb', 'Cwc', 'Dfa', 'Dfb', 'Dfc','Dfd', 'Dsa', 'Dsb', 'Dsc', 'Dsd','Dwa', 'Dwb', 'Dwc', 'Dwd', 'EF','ET', 'Ocean')

# Remove the placeholders
r[1:32] <- r0; levels(r) <- rat

# Select region (Australia)
# x1=80; x2=180; y1=-50; y2=20; xat=5; yat=5	
# Select region (Europe)
# x1=-20; x2=80; y1=30; y2=75; xat=5; yat=5		
# Select region (US)
# x1=-130; x2=-60; y1=20; y2=60; xat=5; yat=5
# Select region (Global)
x1=-10; x2=40; y1=28; y2=48; xat=5; yat=5

r <- crop(r, extent(x1, x2, y1, y2))

# Visualization		
if(.Platform$OS.type=="windows") {quartz<-function(...) windows(...)}
quartz(width=13, height=10, dpi=100)

print(levelplot(r, col.regions=climate.colors, xlab="", ylab="", 
                scales=list(x=list(limits=c(xmin(r), xmax(r)), at=seq(xmin(r), xmax(r), xat)), 
                            y=list(limits=c(ymin(r), ymax(r)), at=seq(ymin(r), ymax(r), yat))), colorkey=list(space="top", tck=0, maxpixels=ncell(r)))
      +layer(sp.polygons(countriesHigh, lwd=0.25)))

out=paste('KG_', period,'_5m.pdf', sep='')
dev.copy2pdf(file=out)

# Find the climate class for Vienna, Europe (or another location)
lon=16.375; lat=48.210
KG=r[cellFromXY(r, c(lon, lat))]
print(KG)

# Output of ASCCI-Data (numbers correspond to the climate classes of the legend)
r <- crop(r, extent(x1, x1+100, y1, y1+100)) # extent must be within the selected region
z <- rasterToPoints(r, spatial=T); z <- spTransform(z, CRS=projection(r))
z <- as.data.frame(z); print(length(t(z[1]))); z = subset(z, z[1]!=32); print(length(t(z[1])))
names(z)=c('KG', 'lon', 'lat')
pts <- data.frame(lat=format(z[2], digits=4), lon=format(z[3], digits=7), KG=format(z[1], digits=3))
write.csv(pts, file=paste('KG_', period,'_5m.csv', sep=''), row.names=F)


#Convert KG to XY then extract to grid points
setwd(path) #change it



z.grid <- rasterToPoints(r, spatial=T); z.grid <- spTransform(z.grid, CRS=projection(r))


#

GRID_latlon <- spTransform(GRID, CRSobj = CRS("+init=epsg:4326"))
GRID_latlon_coordinates = data.frame( coordinates(GRID_latlon))
colnames(GRID_latlon_coordinates)<-c("lon", "lat")
GRID_latlon_df<- cbind(GRID_latlon, GRID_latlon_coordinates)

GRID_KG = extract(r, GRID_latlon_df)
nrow(GRID_KG)
GRID$KG <- GRID_KG

GRID_df <- as.data.frame(GRID)

#change KG values to text:

GRID_df$KG[GRID_df$KG == 1] <-'Af'
GRID_df$KG[GRID_df$KG == 2] <-'Am' 
GRID_df$KG[GRID_df$KG == 3] <-'As' 
GRID_df$KG[GRID_df$KG == 4] <-'Aw' 
GRID_df$KG[GRID_df$KG == 5] <-'BSh' 
GRID_df$KG[GRID_df$KG == 6] <-'BSk' 
GRID_df$KG[GRID_df$KG == 7] <-'BWh'
GRID_df$KG[GRID_df$KG == 8] <-'BWk' 
GRID_df$KG[GRID_df$KG == 9] <-'Cfa'
GRID_df$KG[GRID_df$KG == 10] <-'Cfb'
GRID_df$KG[GRID_df$KG == 11] <-'Cfc' 
GRID_df$KG[GRID_df$KG == 12] <-'Csa'
GRID_df$KG[GRID_df$KG == 13] <-'Csb' 
GRID_df$KG[GRID_df$KG == 14] <-'Csc'
GRID_df$KG[GRID_df$KG == 15] <-'Cwa'
GRID_df$KG[GRID_df$KG == 16] <-'Cwb'
GRID_df$KG[GRID_df$KG == 17] <-'Cwc' 
GRID_df$KG[GRID_df$KG == 18] <-'Dfa' 
GRID_df$KG[GRID_df$KG == 19] <-'Dfb'
GRID_df$KG[GRID_df$KG == 20] <-'Dfc'
GRID_df$KG[GRID_df$KG == 21] <-'Dfd' 
GRID_df$KG[GRID_df$KG == 22] <-'Dsa' 
GRID_df$KG[GRID_df$KG == 23] <-'Dsb' 
GRID_df$KG[GRID_df$KG == 24] <-'Dsc'
GRID_df$KG[GRID_df$KG == 25] <-'Dsd'
GRID_df$KG[GRID_df$KG == 26] <-'Dwa' 
GRID_df$KG[GRID_df$KG == 27] <-'Dwb' 
GRID_df$KG[GRID_df$KG == 28] <-'Dwc' 
GRID_df$KG[GRID_df$KG == 29] <-'Dwd'
GRID_df$KG[GRID_df$KG == 30] <-'EF'
GRID_df$KG[GRID_df$KG == 31] <-'ET' 
GRID_df$KG[GRID_df$KG == 32] <-'Ocean'

#convert it to GRID
coordinates(GRID_df) <- ~ x + y
proj4string(GRID_df) <- CRS("+init=epsg:3857")
# coerce to SpatialPixelsDataFrame
gridded(GRID_df) <- TRUE
GRID <- as(GRID_df, "SpatialPixelsDataFrame")

GRID <- GRID_df



#convert the grid to raster to extract the KG values for the stations:

GRID_KG_df <- data.frame(GRID_coordinates, GRID_KG)

coordinates(GRID_KG_df) <- ~ x + y
proj4string(GRID_KG_df) <- CRS("+init=epsg:3857")
# coerce to SpatialPixelsDataFrame
gridded(GRID_KG_df) <- TRUE
# coerce to raster
GRID_raster <- raster(GRID_KG_df)


################################################################################################################################

#Ez csak a random forestnek kell
nodesizeparam = 6
mtryparam = 20


{
  #Import the station data:
  full_Data  <-readxl::read_excel("122058.xlsx",col_names=TRUE, na="NA", sheet="set1_v2")
  colnames(full_Data) <- c(  "Valid. Categ.","Country","Site","Latitude","Longitude","Y","X","Altitude",        
                             "RMA_slope","RMA_intercept","r2","p","Nr of pairs (O,H)","Data from how many years",
                             "QC"  )
  
  #Extract KG values: 
  full_Data_KG  <- full_Data
  coordinates(full_Data_KG ) <- ~ X + Y
  proj4string(full_Data_KG ) <- CRS("+init=epsg:3857")
  full_Data_KG_values <- extract (GRID_raster, full_Data_KG)
  full_Data$KG <- full_Data_KG_values
  
  
  full_Data$KG[full_Data$KG == 1] <-'Af'
  full_Data$KG[full_Data$KG == 2] <-'Am' 
  full_Data$KG[full_Data$KG == 3] <-'As' 
  full_Data$KG[full_Data$KG == 4] <-'Aw' 
  full_Data$KG[full_Data$KG == 5] <-'BSh' 
  full_Data$KG[full_Data$KG == 6] <-'BSk' 
  full_Data$KG[full_Data$KG == 7] <-'BWh'
  full_Data$KG[full_Data$KG == 8] <-'BWk' 
  full_Data$KG[full_Data$KG == 9] <-'Cfa'
  full_Data$KG[full_Data$KG == 10] <-'Cfb'
  full_Data$KG[full_Data$KG == 11] <-'Cfc' 
  full_Data$KG[full_Data$KG == 12] <-'Csa'
  full_Data$KG[full_Data$KG == 13] <-'Csb' 
  full_Data$KG[full_Data$KG == 14] <-'Csc'
  full_Data$KG[full_Data$KG == 15] <-'Cwa'
  full_Data$KG[full_Data$KG == 16] <-'Cwb'
  full_Data$KG[full_Data$KG == 17] <-'Cwc' 
  full_Data$KG[full_Data$KG == 18] <-'Dfa' 
  full_Data$KG[full_Data$KG == 19] <-'Dfb'
  full_Data$KG[full_Data$KG == 20] <-'Dfc'
  full_Data$KG[full_Data$KG == 21] <-'Dfd' 
  full_Data$KG[full_Data$KG == 22] <-'Dsa' 
  full_Data$KG[full_Data$KG == 23] <-'Dsb' 
  full_Data$KG[full_Data$KG == 24] <-'Dsc'
  full_Data$KG[full_Data$KG == 25] <-'Dsd'
  full_Data$KG[full_Data$KG == 26] <-'Dwa' 
  full_Data$KG[full_Data$KG == 27] <-'Dwb' 
  full_Data$KG[full_Data$KG == 28] <-'Dwc' 
  full_Data$KG[full_Data$KG == 29] <-'Dwd'
  full_Data$KG[full_Data$KG == 30] <-'EF'
  full_Data$KG[full_Data$KG == 31] <-'ET' 
  full_Data$KG[full_Data$KG == 32] <-'Ocean'
  
  
  #Focus on the variables, we need:
  Train_Data = full_Data %>% dplyr::select(RMA_slope,X,Y, KG)
  
  colnames(Train_Data) = c("RMA_slope","x","y", "KG")

  #Random forest
  X1 = data.frame( coordinates(GRID),
                   #elevation = GRID$elevation,
                   KG = GRID$KG)
  X1 = X1 %>% dplyr::select(x,y, KG) # column order should be the same with the data above execpt for the response variable
  Data0.1 = data.frame(RMA_slope = NA, X1);
  Data0.1$KG <- as.factor(Data0.1$KG)
  
  Train_Data$KG <- as.factor(Train_Data$KG)
  levels(Train_Data$KG)  <- levels(Data0.1$KG)
  
  TuneRF2 = tune( RMA_slope ~ ., data = as.data.frame(Train_Data), ntreeTry = 1000)
  RF2 = rfsrc( RMA_slope ~ ., data = as.data.frame(Train_Data),
               nodesize = nodesizeparam,
               mtry = mtryparam, importance = TRUE)
  print(TuneRF2$rf)
  plot(RF2)
  
  #partial dependance plot:
  plot.variable(RF2, target=KG,show.plots = TRUE,  oob = TRUE, partial=FALSE, plots.per.page = 4)
  
  Y0_hat_RF2 = predict(RF2, newdata = Data0.1, na.action = "na.impute")$predicted
  Y0_hat2 = X1 %>% cbind.data.frame( Pred_RF = Y0_hat_RF2)
  
  
  
  ###############################################################################################
  
  #Validate - imput validation stations:
  valid_data  <-readxl::read_excel("122058.xlsx",col_names=TRUE, na="NA", sheet="valid1v2")
  valid_data<- valid_data  %>% drop_na(Site)
  
  #convert to spatial object:
  valid_data2 <- valid_data
  
  colnames(valid_data2) <- c(   "Valid. Categ.",            "Country" ,                 "Site"  ,                   "Latitude"   ,             
                                "Longitude"   ,             "Y"  ,                      "X"     ,                   "Altitude"  ,              
                                "RMA_slope",                "RMA_intercept"     ,       "r2"       ,                "p"         ,              
                                "Nr of pairs"     ,         "Data from how many years" ,"QC"          ,             "intercept"  ,             
                                "...17"    ,                "...18"   ,                 "...19"       ,             "slope"      ,             
                                "...21"            ,        "...22"    ,               "...23"   )
  
  
  ##############################################################################################
  #valid_data2 <- valid_data2[-c(2), ]
  ##############################################################################################
  
  #Extract KG values: 
  valid_data2_KG  <- valid_data2
  coordinates(valid_data2_KG ) <- ~ X + Y
  proj4string(valid_data2_KG ) <- CRS("+init=epsg:3857")
  valid_data2_KG_values <- extract (GRID_raster, valid_data2_KG)
  valid_data2$KG <- valid_data2_KG_values
  
  
  valid_data2$KG[valid_data2$KG == 1] <-'Af'
  valid_data2$KG[valid_data2$KG == 2] <-'Am' 
  valid_data2$KG[valid_data2$KG == 3] <-'As' 
  valid_data2$KG[valid_data2$KG == 4] <-'Aw' 
  valid_data2$KG[valid_data2$KG == 5] <-'BSh' 
  valid_data2$KG[valid_data2$KG == 6] <-'BSk' 
  valid_data2$KG[valid_data2$KG == 7] <-'BWh'
  valid_data2$KG[valid_data2$KG == 8] <-'BWk' 
  valid_data2$KG[valid_data2$KG == 9] <-'Cfa'
  valid_data2$KG[valid_data2$KG == 10] <-'Cfb'
  valid_data2$KG[valid_data2$KG == 11] <-'Cfc' 
  valid_data2$KG[valid_data2$KG == 12] <-'Csa'
  valid_data2$KG[valid_data2$KG == 13] <-'Csb' 
  valid_data2$KG[valid_data2$KG == 14] <-'Csc'
  valid_data2$KG[valid_data2$KG == 15] <-'Cwa'
  valid_data2$KG[valid_data2$KG == 16] <-'Cwb'
  valid_data2$KG[valid_data2$KG == 17] <-'Cwc' 
  valid_data2$KG[valid_data2$KG == 18] <-'Dfa' 
  valid_data2$KG[valid_data2$KG == 19] <-'Dfb'
  valid_data2$KG[valid_data2$KG == 20] <-'Dfc'
  valid_data2$KG[valid_data2$KG == 21] <-'Dfd' 
  valid_data2$KG[valid_data2$KG == 22] <-'Dsa' 
  valid_data2$KG[valid_data2$KG == 23] <-'Dsb' 
  valid_data2$KG[valid_data2$KG == 24] <-'Dsc'
  valid_data2$KG[valid_data2$KG == 25] <-'Dsd'
  valid_data2$KG[valid_data2$KG == 26] <-'Dwa' 
  valid_data2$KG[valid_data2$KG == 27] <-'Dwb' 
  valid_data2$KG[valid_data2$KG == 28] <-'Dwc' 
  valid_data2$KG[valid_data2$KG == 29] <-'Dwd'
  valid_data2$KG[valid_data2$KG == 30] <-'EF'
  valid_data2$KG[valid_data2$KG == 31] <-'ET' 
  valid_data2$KG[valid_data2$KG == 32] <-'Ocean'
  
  
  
  valid_data3 = valid_data2 %>% dplyr::select(RMA_slope, X, Y, KG) # response variables (RMA slope / RMA intercept) and other predictors.
  colnames(valid_data3 ) = c("RMA_slope", "x", "y", "KG")
  
  valid_data_spatial <- valid_data3
  
  coordinates(valid_data_spatial) <- ~ x + y
  proj4string(valid_data_spatial) <- CRS("+init=epsg:3857")
  
  
  
  #Random forest prediction to the validation points:
  
  X1_v = data.frame( coordinates(valid_data_spatial),
                     #elevation = valid_data_spatial$elevation,
                     KG=valid_data_spatial$KG)
  X1_v = X1_v %>% dplyr::select(x,y, KG) # column order should be the same with the data above execpt for the response variable
  Data0.1_v = data.frame(RMA_slope = NA, X1_v);
  
  
  Data0.1_v$KG <- as.factor(Data0.1_v$KG)
  
  Y0_hat_RF2_v = predict(RF2, newdata = Data0.1_v, na.action = "na.impute")$predicted
  
  
  Y0_hat2_v = X1_v %>% cbind.data.frame( Pred_RF = Y0_hat_RF2_v)
  
  
  #calculate the MSE:
  
  se = ((valid_data3$RMA_slope-Y0_hat2_v$Pred_RF)^2)
  MSE = mean(se)
  
  
  valid_data3$Predicted_RMA_slope_X_Y_KG <- Y0_hat2_v$Pred_RF
  
  #Random forest with distances /Hengl Method

  mapstation <- Train_Data
  coordinates(mapstation) <- ~x + y
  proj4string(mapstation) <- CRS("+init=epsg:3857")
 
  # if you get error like the raster values has NA values, bla bla ; then delete those stations, that are not inside the GRID boundry ;)
  grid.dist0 <- GSIF::buffer.dist(mapstation["RMA_slope"], GRID[1],  as.factor(1:nrow(mapstation)))
  
  #grid.dist0$elevation = extract(rasDEM_europe2, grid.dist0)
  grid.dist0$KG = extract(GRID_raster, grid.dist0)
  
  #convert the KG values to factors:
  
  
  grid.dist0$KG[grid.dist0$KG == 1] <-'Af'
  grid.dist0$KG[grid.dist0$KG == 2] <-'Am' 
  grid.dist0$KG[grid.dist0$KG == 3] <-'As' 
  grid.dist0$KG[grid.dist0$KG == 4] <-'Aw' 
  grid.dist0$KG[grid.dist0$KG == 5] <-'BSh' 
  grid.dist0$KG[grid.dist0$KG == 6] <-'BSk' 
  grid.dist0$KG[grid.dist0$KG == 7] <-'BWh'
  grid.dist0$KG[grid.dist0$KG == 8] <-'BWk' 
  grid.dist0$KG[grid.dist0$KG == 9] <-'Cfa'
  grid.dist0$KG[grid.dist0$KG == 10] <-'Cfb'
  grid.dist0$KG[grid.dist0$KG == 11] <-'Cfc' 
  grid.dist0$KG[grid.dist0$KG == 12] <-'Csa'
  grid.dist0$KG[grid.dist0$KG == 13] <-'Csb' 
  grid.dist0$KG[grid.dist0$KG == 14] <-'Csc'
  grid.dist0$KG[grid.dist0$KG == 15] <-'Cwa'
  grid.dist0$KG[grid.dist0$KG == 16] <-'Cwb'
  grid.dist0$KG[grid.dist0$KG == 17] <-'Cwc' 
  grid.dist0$KG[grid.dist0$KG == 18] <-'Dfa' 
  grid.dist0$KG[grid.dist0$KG == 19] <-'Dfb'
  grid.dist0$KG[grid.dist0$KG == 20] <-'Dfc'
  grid.dist0$KG[grid.dist0$KG == 21] <-'Dfd' 
  grid.dist0$KG[grid.dist0$KG == 22] <-'Dsa' 
  grid.dist0$KG[grid.dist0$KG == 23] <-'Dsb' 
  grid.dist0$KG[grid.dist0$KG == 24] <-'Dsc'
  grid.dist0$KG[grid.dist0$KG == 25] <-'Dsd'
  grid.dist0$KG[grid.dist0$KG == 26] <-'Dwa' 
  grid.dist0$KG[grid.dist0$KG == 27] <-'Dwb' 
  grid.dist0$KG[grid.dist0$KG == 28] <-'Dwc' 
  grid.dist0$KG[grid.dist0$KG == 29] <-'Dwd'
  grid.dist0$KG[grid.dist0$KG == 30] <-'EF'
  grid.dist0$KG[grid.dist0$KG == 31] <-'ET' 
  grid.dist0$KG[grid.dist0$KG == 32] <-'Ocean'
  
  grid.dist0$KG <- as.factor(grid.dist0$KG)
  levels(mapstation@data[["KG"]])  <- levels(grid.dist0$KG)
  
  
  dn0 <- paste(names(grid.dist0), collapse="+")
  fm0 <- as.formula(paste("RMA_slope ~ ", dn0))
  ov.RMA_slope <- over(mapstation["RMA_slope"], grid.dist0)
  rm.RMA_slope <- cbind(mapstation@data["RMA_slope"], ov.RMA_slope)
  
  # Step 0-2. RFSRC
  rm.RMA_slope$KG <- as.factor(rm.RMA_slope$KG)
  levels(rm.RMA_slope$KG)  <- levels(grid.dist0$KG)
  
  #Since RFSRC gave back for Nodesize a 15, I used the Ranger package and manually test all Mtry from 1-53 and Node size from 1-15 to get the lowers OOB error.
  #The optimal parameters in the current dataset for Nodesize is 6 and Mtry is sqrt 54 = 7, and increasing it, there is a drop of OOB at Mtr 14 so I set it to 14
  #devtools::install_github("imbs-hl/ranger")
  
  
  library(ranger)
  RF_buffer_distance <- ranger(fm0, rm.RMA_slope, quantreg=TRUE, 
                               num.trees=1000, 
                               importance="impurity",
                               oob.error=TRUE,
                               mtry = mtryparam,
                               min.node.size = nodesizeparam,
                               write.forest = TRUE
  ) 
  
  library(xlsx)
  importance <- as.list(ranger::importance(RF_buffer_distance))
  print(t(data.frame(importance[order(unlist(importance), decreasing=TRUE)[1:52]]))) 
  RF_buffer_distance[["prediction.error"]]
  
  importance_df = data.frame((t(data.frame(importance[order(unlist(importance), decreasing=TRUE)[1:52]]))) )
  OOB_df = data.frame(RF_buffer_distance[["prediction.error"]])
  
  write.xlsx(importance_df,file="maps/final/v4_ranger/importance_df_slope_set1v2.xlsx")
  write.xlsx(OOB_df,file="maps/final/v4_ranger/OOB_df_slope_set1v2.xlsx")
  
  ####################try modify tune
  #print(TuneRF$rf)
  
  #TuneRF_results_df <- data.frame(TuneRF[["results"]])
  #plot(TuneRF_results_df$nodesize, TuneRF_results_df$err)
  
  #library(rgl)
  #plot3d( 
  #  x=TuneRF_results_df$nodesize, y= TuneRF_results_df$err, z= TuneRF_results_df$mtry,
  
  #  xlab="nodesize", ylab="OOB error", zlab="mtry")
  
  #htmlwidgets::saveWidget(rglwidget(width = 1000, height = 1000), 
  #                        file = "3dscatter.html",
  #                        libdir = "libs",
  #                        selfcontained = FALSE
  #)
  
  
  ##########################
  grid.dist0$KG <- as.factor(grid.dist0$KG)
  
  
  Y0_hat_RF_buffer_distance = predict(RF_buffer_distance, grid.dist0@data, na.action = "na.impute")
  
  
  #mapview(Y0_hat_buffer_distance, alpha=1) + mapview(mapstation, zcol = "RMA_slope" ,alpha = 1)
  
  
  #Validate the buffer distance method
  #We can not use the validation set's coordinates to validate since buffer distance require an equal grid array to calculate distances.
  #In this case we will search for the nearest predicted Gridpoint to the validation stations
  
  #first step to convert the grid to raster
  
  Y0_hat_buffer_distance_df <- as.data.frame(GRID_coordinates)
  Y0_hat_buffer_distance_df$RMA_slope <- Y0_hat_RF_buffer_distance[["predictions"]]
  coordinates(Y0_hat_buffer_distance_df) <- ~ x + y
  proj4string(Y0_hat_buffer_distance_df) <- CRS("+init=epsg:3857")
  gridded(Y0_hat_buffer_distance_df) <- TRUE
  
  Y0_hat_buffer_distance_df <- raster(Y0_hat_buffer_distance_df)
  Y0_hat_buffer_distance_validations <- raster::extract(Y0_hat_buffer_distance_df,             # raster layer
                                                        valid_data_spatial,   # SPDF
                                                        df=TRUE)         # return a dataframe? 
  
  #add results to the validation df
  valid_data3$Predicted_RMA_slope_X_Y_KG_Buffer_distance <- Y0_hat_buffer_distance_validations$RMA_slope
  
  Raster1 <- as.data.frame(coordinates(GRID))
  Raster1$RMA_Slope <- Y0_hat2$Pred_RF
  coordinates(Raster1) <- ~ x + y
  proj4string(Raster1) <- CRS("+init=epsg:3857")
  gridded(Raster1) <- TRUE
  Raster1 <- raster(Raster1)
  
  
  #save the results:
  write.xlsx(valid_data3,file="maps/final/v4_ranger/RMA_slope_x_y_KG_set1v2.xlsx")
  writeRaster(Raster1, filename="maps/final/v4_ranger/RMA_slope_x_y_KG_set1v2.tif", overwrite = TRUE)
  writeRaster(Y0_hat_buffer_distance_df, filename="maps/final/v4_ranger/RMA_slope_x_y_KG_Buffer_distance_set1v2.tif",overwrite = TRUE)
  
  #######################################################################################################################################
  #RMA_Intercept
  
  
  
  #Focus on the variables, we need:
  Train_Data = full_Data %>% dplyr::select(RMA_intercept,X,Y, KG) # response variables (RMA slope / RMA intercept) and other predictors.
  colnames(Train_Data) = c("RMA_intercept","x","y", "KG")
  
  #Random forest
  X1 = data.frame( coordinates(GRID),
                   #elevation = GRID$elevation,
                   KG = GRID$KG)
  X1 = X1 %>% dplyr::select(x,y, KG) # column order should be the same with the data above execpt for the response variable
  Data0.1 = data.frame(RMA_intercept = NA, X1);
  Data0.1$KG <- as.factor(Data0.1$KG)
  
  Train_Data$KG <- as.factor(Train_Data$KG)
  levels(Train_Data$KG)  <- levels(Data0.1$KG)
  
  
  # nodesize = TuneRF2$optimal['nodesize']
  #mtry = TuneRF2$optimal['mtry']
  TuneRF2 = tune( RMA_intercept ~ ., data = as.data.frame(Train_Data), ntreeTry = 1000)
  RF2 = rfsrc( RMA_intercept ~ ., data = as.data.frame(Train_Data),
               nodesize = nodesizeparam,
               mtry = mtryparam, importance = TRUE)
  plot(RF2)
  #plot(get.tree(RF2, 3))
  v.max <- max.subtree(RF2)
  print(round(v.max$order, 3))
  print(round(v.max$order[, 1], 3))
  print(v.max$threshold)
  print(vimp(RF2)$importance)
  plot.variable(RF2, target=KG,show.plots = TRUE,  oob = TRUE, partial=FALSE, plots.per.page = 4)
  #partial.rfsrc(RF2, oob = TRUE)
  
  Y0_hat_RF2 = predict(RF2, newdata = Data0.1, na.action = "na.impute")$predicted
  
  
  Y0_hat2 = X1 %>% cbind.data.frame( Pred_RF = Y0_hat_RF2)

  valid_data3 = valid_data2 %>% dplyr::select(RMA_intercept, X, Y, KG) # response variables (RMA slope / RMA intercept) and other predictors.
  colnames(valid_data3 ) = c("RMA_intercept", "x", "y","KG")
  
  valid_data_spatial <- valid_data3
  
  coordinates(valid_data_spatial) <- ~ x + y
  proj4string(valid_data_spatial) <- CRS("+init=epsg:3857")

  #Random forest prediction to the validation points:
  
  X1_v = data.frame( coordinates(valid_data_spatial),
                     #elevation = valid_data_spatial$elevation,
                     KG=valid_data_spatial$KG)
  X1_v = X1_v %>% dplyr::select(x,y, KG) # column order should be the same with the data above execpt for the response variable
  Data0.1_v = data.frame(RMA_intercept = NA, X1_v);
  
  
  Data0.1_v$KG <- as.factor(Data0.1_v$KG)
  
  Y0_hat_RF2_v = predict(RF2, newdata = Data0.1_v, na.action = "na.impute")$predicted
  
  
  Y0_hat2_v = X1_v %>% cbind.data.frame( Pred_RF = Y0_hat_RF2_v)
  
  
  #calculate the MSE:
  
  se = ((valid_data3$RMA_intercept-Y0_hat2_v$Pred_RF)^2)
  MSE = mean(se)
  
  
  valid_data3$Predicted_RMA_intercept_X_Y_KG <- Y0_hat2_v$Pred_RF

  mapstation <- Train_Data
  coordinates(mapstation) <- ~x + y
  proj4string(mapstation) <- CRS("+init=epsg:3857")

  grid.dist0 <- GSIF::buffer.dist(mapstation["RMA_intercept"], GRID[1],  as.factor(1:nrow(mapstation)))
  #grid.dist0$elevation = extract(rasDEM_europe2, grid.dist0)
  grid.dist0$KG = extract(GRID_raster, grid.dist0)
  
  #convert the KG values to factors:
 
  grid.dist0$KG[grid.dist0$KG == 1] <-'Af'
  grid.dist0$KG[grid.dist0$KG == 2] <-'Am' 
  grid.dist0$KG[grid.dist0$KG == 3] <-'As' 
  grid.dist0$KG[grid.dist0$KG == 4] <-'Aw' 
  grid.dist0$KG[grid.dist0$KG == 5] <-'BSh' 
  grid.dist0$KG[grid.dist0$KG == 6] <-'BSk' 
  grid.dist0$KG[grid.dist0$KG == 7] <-'BWh'
  grid.dist0$KG[grid.dist0$KG == 8] <-'BWk' 
  grid.dist0$KG[grid.dist0$KG == 9] <-'Cfa'
  grid.dist0$KG[grid.dist0$KG == 10] <-'Cfb'
  grid.dist0$KG[grid.dist0$KG == 11] <-'Cfc' 
  grid.dist0$KG[grid.dist0$KG == 12] <-'Csa'
  grid.dist0$KG[grid.dist0$KG == 13] <-'Csb' 
  grid.dist0$KG[grid.dist0$KG == 14] <-'Csc'
  grid.dist0$KG[grid.dist0$KG == 15] <-'Cwa'
  grid.dist0$KG[grid.dist0$KG == 16] <-'Cwb'
  grid.dist0$KG[grid.dist0$KG == 17] <-'Cwc' 
  grid.dist0$KG[grid.dist0$KG == 18] <-'Dfa' 
  grid.dist0$KG[grid.dist0$KG == 19] <-'Dfb'
  grid.dist0$KG[grid.dist0$KG == 20] <-'Dfc'
  grid.dist0$KG[grid.dist0$KG == 21] <-'Dfd' 
  grid.dist0$KG[grid.dist0$KG == 22] <-'Dsa' 
  grid.dist0$KG[grid.dist0$KG == 23] <-'Dsb' 
  grid.dist0$KG[grid.dist0$KG == 24] <-'Dsc'
  grid.dist0$KG[grid.dist0$KG == 25] <-'Dsd'
  grid.dist0$KG[grid.dist0$KG == 26] <-'Dwa' 
  grid.dist0$KG[grid.dist0$KG == 27] <-'Dwb' 
  grid.dist0$KG[grid.dist0$KG == 28] <-'Dwc' 
  grid.dist0$KG[grid.dist0$KG == 29] <-'Dwd'
  grid.dist0$KG[grid.dist0$KG == 30] <-'EF'
  grid.dist0$KG[grid.dist0$KG == 31] <-'ET' 
  grid.dist0$KG[grid.dist0$KG == 32] <-'Ocean'
  
  grid.dist0$KG <- as.factor(grid.dist0$KG)
  levels(mapstation@data[["KG"]])  <- levels(grid.dist0$KG)
  
  
  dn0 <- paste(names(grid.dist0), collapse="+")
  fm0 <- as.formula(paste("RMA_intercept ~ ", dn0))
  ov.RMA_intercept <- over(mapstation["RMA_intercept"], grid.dist0)
  rm.RMA_intercept <- cbind(mapstation@data["RMA_intercept"], ov.RMA_intercept)
  
  # Step 0-2. RFSRC
  rm.RMA_intercept$KG <- as.factor(rm.RMA_intercept$KG)
  levels(rm.RMA_intercept$KG)  <- levels(grid.dist0$KG)

  RF_buffer_distance <- ranger(fm0, rm.RMA_intercept, quantreg=TRUE, 
                               num.trees=1000, 
                               importance="impurity",
                               oob.error=TRUE,
                               mtry = mtryparam,
                               min.node.size = nodesizeparam,
                               write.forest = TRUE
  ) 
  
  
  importance <- as.list(ranger::importance(RF_buffer_distance))
  print(t(data.frame(importance[order(unlist(importance), decreasing=TRUE)[1:52]]))) 
  RF_buffer_distance[["prediction.error"]]
  
  importance_df = data.frame((t(data.frame(importance[order(unlist(importance), decreasing=TRUE)[1:52]]))) )
  OOB_df = data.frame(RF_buffer_distance[["prediction.error"]])
  
  write.xlsx(importance_df,file="maps/final/v4_ranger/importance_df_intercept_set1v2.xlsx")
  write.xlsx(OOB_df,file="maps/final/v4_ranger/OOB_df_intercept_set1v2.xlsx")
  
  
  
  grid.dist0$KG <- as.factor(grid.dist0$KG)
  
  
  Y0_hat_RF_buffer_distance = predict(RF_buffer_distance, grid.dist0@data, na.action = "na.impute")

  #Validate the buffer distance method
  #We can not use the validation set's coordinates to validate since buffer distance require an equal grid array to calculate distances.
  #In this case we will search for the nearest predicted Gridpoint to the validation stations
  
  #first step to convert the grid to raster
  
  Y0_hat_buffer_distance_df <- as.data.frame(GRID_coordinates)
  Y0_hat_buffer_distance_df$RMA_intercept <- Y0_hat_RF_buffer_distance[["predictions"]]
  coordinates(Y0_hat_buffer_distance_df) <- ~ x + y
  proj4string(Y0_hat_buffer_distance_df) <- CRS("+init=epsg:3857")
  gridded(Y0_hat_buffer_distance_df) <- TRUE
  
  Y0_hat_buffer_distance_df <- raster(Y0_hat_buffer_distance_df)
  Y0_hat_buffer_distance_validations <- raster::extract(Y0_hat_buffer_distance_df,             # raster layer
                                                        valid_data_spatial,   # SPDF
                                                        df=TRUE)         # return a dataframe? 
  
  #add results to the validation df
  valid_data3$Predicted_RMA_intercept_X_Y_KG_Buffer_distance <- Y0_hat_buffer_distance_validations$RMA_intercept
  
  Raster1 <- as.data.frame(coordinates(GRID))
  Raster1$RMA_intercept <- Y0_hat2$Pred_RF
  coordinates(Raster1) <- ~ x + y
  proj4string(Raster1) <- CRS("+init=epsg:3857")
  gridded(Raster1) <- TRUE
  Raster1 <- raster(Raster1)
  
  write.xlsx(valid_data3,file="maps/final/v4_ranger/RMA_intercept_results_x_y_KG_set1v2.xlsx")
  writeRaster(Raster1, filename="maps/final/v4_ranger/RMA_intercept_x_y_KG_set1v2.tif", overwrite = TRUE)
  writeRaster(Y0_hat_buffer_distance_df, filename="maps/final/v4_ranger/RMA_intercept_x_y_KG_Buffer_distance_set1v2.tif",overwrite = TRUE)
  
}