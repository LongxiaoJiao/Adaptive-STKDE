

# excel --> dataframe --> t_points_matrix --> TSTKDE & ASTKDE (joint model) --> TSTKDE & ASTKDE (spatial conditional model) --> hotspots & peaks & trajectories

#============================== MAIN FUNCTIONS =================================


# The user provides a t_points_matrix that contains the observation points’ X, Y, and T coordinates, 
# as well as the spatial and temporal bandwidths, which are used to compute STKDE.
STKDE.joint.constructor.with.bandwidths <- function(t_points_matrix,metric_list, parallel = FALSE){
  
  #INPUT t_points_matrix - training set points matrix
  
  #  1st column  2nd column 3rd column        4th column                5th column
  #       X          Y          Z         spatial bandwidths        temporal bandwidths
  
  #  OR break the isotropy in XY plane
  
  #  1st column  2nd column 3rd column      4th column      5th column          6th column
  #       X          Y          Z         X bandwidths       Y bandwidths     temporal bandwidth
  
  #      metric_list    - include all key parameters of the spatiotemporal cube
  #      parallel       - TRUE or FALSE
  
  
  #OUTPUT q            - edge correction factor
  #       result_STKDE - 3D array of ASTKDE OR TSTKDE
  
  
  n <- nrow(t_points_matrix)
  
  grx <- metric_list$grx
  gry <- metric_list$gry
  grz <- metric_list$grz
  xstep <- metric_list$xstep
  ystep <- metric_list$ystep
  time_unit <- metric_list$time_unit

  result_STKDE <- q <- NULL
  
  start_time <- Sys.time()
  
  if(!parallel) {
    
    # single thread
    
    numeric_01_cube <- outer(metric_list$numeric_01_matrix, rep(1,metric_list$dimensionxyz[3]), "*")
    # numeric_01_cube : The indicator function I_V determines whether the estimation location(x,y,t) falls within the domain V, returning 1 if true and 0 otherwise.
    
    STKDE <- array(0,metric_list$dimensionxyz[c(2,1,3)])
    
    q <- apply(t_points_matrix,1,function(row) {
      
      mx <- dnorm(grx,row[1],row[4])
      my <- dnorm(gry,row[2],row[4])
      mz <- dnorm(grz,row[3],row[5])
      
      my.mx <- outer(my,mx,"*")
      my.mx.mz <- outer(my.mx,mz,"*")
      
      my.mx.mz <- my.mx.mz * numeric_01_cube
      
      q <- DATA.integral3d(my.mx.mz,xstep,ystep,time_unit)
      
      STKDE <<- STKDE + my.mx.mz / q
      
      return(q)
      
    })
    
    result_STKDE <- STKDE/n
    
  }
  
  else{
    
    # double thread
    
    mid <- floor(n/2)
    
    t_points_matrix_A <- t_points_matrix[1:mid, , drop = FALSE]
    t_points_matrix_B <- t_points_matrix[(mid+1):n, , drop = FALSE]
    
    cl <- makeCluster(2)
    
    clusterExport(cl, varlist = c("DATA.integral3d","makeCluster","clusterExport","parLapply","stopCluster"))
    
    result_list <- parLapply(cl,list(t_points_matrix_A, t_points_matrix_B), 
                        
                        function(matrix){
                          
                          STKDE <- array(0,metric_list$dimensionxyz[c(2,1,3)])
                          
                          numeric_01_cube <- outer(metric_list$numeric_01_matrix, rep(1,metric_list$dimensionxyz[3]), "*")
                          
                          q <- apply(t_points_matrix,1,function(row) {
                            
                            mx <- dnorm(grx,row[1],row[4])
                            my <- dnorm(gry,row[2],row[4])
                            mz <- dnorm(grz,row[3],row[5])
                            
                            my.mx <- outer(my,mx,"*")
                            my.mx.mz <- outer(my.mx,mz,"*")
                            
                            my.mx.mz <- my.mx.mz * numeric_01_cube  # numeric_01_cube : The indicator function I_V determines whether the estimation location(x,y,t) falls within the domain V, returning 1 if true and 0 otherwise.
                            
                            q <- DATA.integral3d(my.mx.mz,xstep,ystep,time_unit)
                            
                            STKDE <<- STKDE + my.mx.mz / q
                            
                            return(q)
                            
                          })
                          
                          return(list(q=q, result_STKDE=STKDE))
                          
                        })
    stopCluster(cl)
    
    q <- c(result_list[[1]]$q, result_list[[2]]$q)
    
    result_STKDE <- (result_list[[1]]$result_STKDE + result_list[[2]]$result_STKDE)/n
    
  }
  
  result_STKDE[result_STKDE == 0] <- NA
  
  cat("Done!", difftime(Sys.time(),start_time,units="mins"),"\n")
  
  return(list(q=q,result_STKDE=result_STKDE))
  
}


# The user provides a t_points_matrix, which is then used with the plug-in bandwidth selector
# to compute the corresponding TSTKDE and ASTKDE.
STKDE.PI.constructor <- function(dataframe, metric_list, parallel = FALSE){
  
  # INPUT - dataframe includes X,Y,Z
  
  # OUTPUT - 3D array of TSTKDE and ASTKDE (PI)
  
  # You may replace the PI selector to other bandwidth selectors you prefer.
  
  n <- nrow(dataframe)
  
  TSTKDE <- ASTKDE <- array(0,metric_list$dimensionxyz[c(2,1,3)]) # Y, X, Z
  
  t_points_matrix <- cbind(dataframe$X, dataframe$Y, dataframe$Z)
  
  # use ks::Hpi.diag as the plug-in selector to get fixed bandwidths
  
  bw_PI <- sqrt(diag((Hpi.diag(t_points_matrix,nstage = 2,pilot="dscalar",verbose=FALSE,
                               optim.fun="optim",binned = FALSE))))
  hp_PI <- mean(bw_PI[c(1,2)])
  lambda_PI<- bw_PI[3]
  
  t_points_matrix <- cbind(t_points_matrix, rep(hp_PI,n), rep(lambda_PI,n))

  TSTKDE_result <- STKDE.joint.constructor.with.bandwidths(t_points_matrix, metric_list, parallel)
  
  q <- TSTKDE_result$q
  
  # for adaptive bandwidths
  
  fhat_vector <- STKDE.q.pointsmatrix.to.fhatvector(t_points_matrix,t_points_matrix,q)
  
  gamma <- exp(mean(log(fhat_vector^order_2)))
  h_adaptive <- hp_PI * (fhat_vector^order_2)/gamma
  lambda_adaptive <- lambda_PI * (fhat_vector^order_2)/gamma
  
  t_points_matrix <- cbind(t_points_matrix[,c(1,2,3)],h_adaptive,lambda_adaptive)
  
  ASTKDE_result <- STKDE.joint.constructor.with.bandwidths(t_points_matrix, metric_list, parallel)
  
  return(list(TSTKDE = TSTKDE_result$result_STKDE, ASTKDE = ASTKDE_result$result_STKDE))
  
}


#============================== SUB FUNCTIONS =================================


# Using q and t_points_matrix, compute the fhat value for every point in v_points_matrix.
STKDE.q.pointsmatrix.to.fhatvector <- function(t_points_matrix, v_points_matrix, q) {
  

  
  #INPUT  q - edge correction factor
  #       t_points_matrix - observation points matrix (training set points matrix)
  #       v_points_matrix - test set points matrix
  
  #OUTPUT fhat_vector - vector of fhat at each points in test set.
  
  # Sometimes, it’s quite useful to set v_points_matrix equal to t_points_matrix.
  
  tn <- nrow(t_points_matrix)  # number of training set points
  
  mx.my.mz <- matrix(NA, nrow = tn, ncol = tn) # tn * tn matrix
  matrix_x <- apply(t_points_matrix,1,function(row) dnorm(v_points_matrix[,1],row[1],row[4]))
  matrix_y <- apply(t_points_matrix,1,function(row) dnorm(v_points_matrix[,2],row[2],row[4]))
  matrix_z <- apply(t_points_matrix,1,function(row) dnorm(v_points_matrix[,3],row[3],row[5]))
  
  if(is.matrix(q)) {
    
    mx.my <- matrix_x * matrix_y
    mx.my <- t(apply(mx.my,1,function(row) row/q[,1]))
    
    mz <- t(apply(matrix_z, 1, function(row) row/q[,2]))
    
    mx.my.mz <- mx.my * mz
    
  }
  else if(is.vector(q)) {
    
    mx.my.mz <- matrix_x * matrix_y * matrix_z
    mx.my.mz <- t(apply(mx.my.mz,1,function(row) row/q))
    
  }
  else(stop("INPUT q in not a matrix nor a vector."))
  
  fhat_vector <- apply(mx.my.mz,1,function(row) sum(row)/tn)
  
  return(fhat_vector)
  
}


# Given a spatiotemporal joint model (represented by a 3D array in RAM), 
# compute the spatial conditional model by normalizing each layer individually.
STKDE.cond.constructor <- function(STKDE_joint, metric_list) {

  
  #INPUT :  STKDE_joint (joint model)
  
  #OUTPUT : STKDE_cond (spatial conditional model)
  
  xstep <- metric_list$xstep
  ystep <- metric_list$ystep
  
  STKDE_cond <- apply(STKDE_joint,3,function(x) {
    
    x/DATA.integral2d(x, xstep, ystep)
    
  })
  
  STKDE_cond <- array(STKDE_cond, dim(STKDE_joint))
  
  return(STKDE_cond)
  
}

# Using the provided STKDE_cond, compute the temporal marginal density.
STKDE.Temporal.Marginal.density <- function(STKDE_cond, metric_list, alpha_vector = c(0.95,0.99), plot_flag = TRUE) {
  
  
  
  # temporal marginal density of the spatial conditional model
  # plot_flag = TRUE, "Display the density using a red-to-white gradient."， FALSE, "no plot".
  # alpha_vector, vector of significant level, which will draw contourlines.
  
  time_range <- range(metric_list$grz)
  
  result_matrix <- matrix(0, nrow = metric_list$dimensionxyz[2], ncol = metric_list$dimensionxyz[1]) # Y x
  
  apply(STKDE_cond,3,function(z) {
    
    result_matrix <<- result_matrix + z
    
  })

  
  if(plot_flag){
    
    mycolfunc <- colorRampPalette(c("white","red"))
    
    cols <- mycolfunc(100)
    
    density_min <- min(result_matrix,na.rm = TRUE)
    density_max <- max(result_matrix,na.rm = TRUE)
    
    color_indices <- findInterval(result_matrix, seq(density_min, density_max, length.out=101))
    bottom_colors <- matrix(cols[color_indices],nrow = metric_list$dimensionxyz[2], ncol = metric_list$dimensionxyz[1])
    
    z_matrix <- matrix(time_range[1]+10, nrow = metric_list$dimensionxyz[2], ncol = metric_list$dimensionxyz[1]) # Y X
    surface3d(x=metric_list$grx, y=metric_list$gry, z=t(z_matrix), color = t(bottom_colors), lit = FALSE, alpha = 0.5)
    
  }
  
  threshold_value <- quantile(result_matrix, alpha_vector, na.rm = TRUE)

  
  clines <- contourLines(x=metric_list$grx, y=metric_list$gry, z = t(result_matrix), levels = threshold_value)
  
  for(line in clines){
    
    lines3d(line$x,line$y, z = rep(time_range[1]+15,length(line$x)),lwd = 3)
    
  }

}


# Given different user-specified temporal domains, calculate q, which is the edge correction factor.
STKDE.pointsmatrix.to.q <- function(grz_this, metric_list, t_points_matrix){
  
  # INPUT : grz_this, a vector of temporal domain.
  #         metric_list, generated in DATA.R
  #         t_points_matrix, matrix of training points, with bandwidths
  
  # OUTPUT : q, edge correction factor, a vector.
  
  grx <- metric_list$grx
  gry <- metric_list$gry
  
  xstep <- metric_list$xstep
  ystep <- metric_list$ystep
  time_unit <- metric_list$time_unit
  
  numeric_01_cube <- outer(metric_list$numeric_01_matrix,rep(1,length(grz_this)),"*")
  
  q <- apply(t_points_matrix,1,function(row) {
    
    mx <- dnorm(grx,row[1],row[4])
    
    my <- dnorm(gry,row[2],row[4])
    
    mz <- dnorm(grz_this,row[3],row[5])
    
    my.mx <- outer(my,mx,"*")
    
    my.mx.mz <- outer(my.mx,mz,"*")
    
    my.mx.mz <- my.mx.mz * numeric_01_cube
    
    q <- DATA.integral3d(my.mx.mz,xstep,ystep,time_unit)
    
    return(q)
    
  })
  
  return(q)
  
}


#=============================== RUN CODE HERE ===================================

STKDE_result <- STKDE.PI.constructor(realDataFrame,metric_list_half,parallel = TRUE)

TSTKDE_joint <- STKDE_result$TSTKDE
ASTKDE_joint <- STKDE_result$ASTKDE

ASTKDE_cond <- STKDE.cond.constructor(STKDE_result$ASTKDE, metric_list_half)
