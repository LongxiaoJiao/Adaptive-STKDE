#DATA

library(spatstat)
library(ks)
library(parallel)
library(abind)
library(sparr)
library(rgl)
library(rglplus)
library(sf)
library(raster)
library(misc3d)

#============================== MAIN FUNCTIONS =================================


# Construct metric_list using two different resolutions, 
# recording all key parameters of the Spatiotemporal Cube.
DATA.generate.metric.list <- function(dataset, plan_flag=c("full","half")){
  
  boundary_x <- range(ZZBoundaryNEW2$X)
  boundary_y <- range(ZZBoundaryNEW2$Y)
  
  # obwindow option 1
  # obwindow <- owin(xrange=boundary_x,yrange = boundary_y)
  
  # obwindow option 2
  
  obwindow <- owin(poly=list(x=ZZBoundaryNEW2$X,y=ZZBoundaryNEW2$Y))
  
  # obwindow is either the study area (an irregular closed polygon) or a regular rectangle.
  
  time_range <- range(dataset$DATE) * para_c  # 10 years * scale
  
  STdataframe <- data.frame(
    
    X = dataset$X,
    Y = dataset$Y,
    Z = dataset$DATE * para_c # proper scale to T-axis
    
  )
  
  pointpattern <- ppp(x=STdataframe$X, y=STdataframe$Y,window = obwindow)
  
  W <- window(pointpattern)
  
  time_unit <- WM <- grz <- NULL
  
  if(plan_flag == "half"){
    
    time_unit <- 20 * para_c  # 20 days per voxel with eps = c(100,100)
    WM <- as.mask(W,eps=c(100,100))
    grz <- seq(time_range[1] + time_unit/2, time_range[2]+50, by = time_unit)
    
  }
  else if(plan_flag == "full"){
    
    time_unit <- 10 * para_c  # 20 days per voxel with eps = c(50,50)
    WM <- as.mask(W,eps=c(50,50))
    grz <- seq(time_range[1] + time_unit/2, time_range[2], by = time_unit)
    
  }
  else{
    
    stop("plan_flag should be 'half' or 'full'")
    
  }
  
  
  # "As stated in the paper, we have two differentiation schemes, with spatial resolutions of `c(50,50)` and `c(100,100)`, respectively."
  
  xres <- WM$dim[2]
  yres <- WM$dim[1]
  
  
  zres <- length(grz)
  
  
  numeric_01_matrix <- matrix(as.numeric(WM$m),nrow =yres)  # I_V  01 2D matrix 
  
  # The indicator function I_V determines whether the estimation location(x,y,t) falls within the domain V, returning 1 if true and 0 otherwise.
  
  metric_list <- list()
  
  metric_list$grx <- WM$xcol
  metric_list$gry <- WM$yrow
  metric_list$grz <- grz
  
  metric_list$xstep <- WM$xstep
  metric_list$ystep <- WM$ystep
  metric_list$time_unit <- time_unit
  
  metric_list$boundary_x <- boundary_x
  metric_list$boundary_y <- boundary_y
  
  # metric_list$numeric_01_matrix <- numeric_01_matrix
  # metric_list$numeric_01_cube <- numeric_01_cube
  
  # you will need I_V if your study area is irregular polygon
  
  metric_list$dimensionxyz <- c(xres,yres,zres) # 348,368,183
  metric_list$numeric_01_matrix <- numeric_01_matrix
  metric_list$obwindow <- obwindow
  
  metric_list$type = plan_flag
  
  return(metric_list)
  
}


# 2D surface integral.
DATA.integral2d <- function(f,dx,dy){
  
  return(sum(f,na.rm =TRUE)*dx*dy)
  
}


# 3D surface integral.
DATA.integral3d <- function(f,dx,dy,dz){
  
  return(sum(f,na.rm = TRUE)*dx*dy*dz)
  
}


#========================= Global Environment Constants ======================


# the root of the Square Root Law
order_2 <- -0.5  


# the proper scale at T-axis
para_c <- 5 


# dataframe of your real data, which contains X, Y, Z, at least.
realDataFrame <- data.frame(
  
  X = realDataset$X,
  Y = realDataset$Y,
  Z = realDataset$DATE * para_c
  
)


# For prediction only
# This dataframe stores all observation points, representing the complete set of both the training and test datasets.
ALLdataframe <- data.frame(
  
  X = ALLSTdataset$X,
  Y = ALLSTdataset$Y,
  DATE = ALLSTdataset$DATE * para_c
  
)



#=============================== Other Constants Need to Set  Before Start ===================================

# 1st step : "order_2", "para_c", "realDataFrame", "ALLdataframe" need to be loaded in the Global Environment.
#            The Excel file containing your study area needs to be loaded. It must include two columns: X and Y, representing a counterclockwise closed polygon.
#            you need use your own real trivariate data.

order_2
para_c
realDataFrame
ALLdataframe

# 2nd step : metric_list (full or half) must be generated in the Global Environment
# "As stated in the paper, we have two differentiation schemes, with spatial resolutions of ` full : c(50,50)` and `half : c(100,100)`, respectively."
metric_list_half <- DATA.generate.metric.list(realDataset, plan_flag = "half")
metric_list_full <- DATA.generate.metric.list(realDataset, plan_flag = "full")

