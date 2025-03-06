
#============================== MAIN FUNCTIONS =================================


# Construct a list to store all results, during one synthetic pdf experiment.
Synthetic.result.list.constructor <- function(metric_list,n){
  
  # n is the number of the sample points in one stochastic process
  
  result_list <- list()
  result_list$iteration <- numeric(0)
  result_list$method <- "PI & LIK & SCV"
  result_list$n <- n
  
  result_list$bw_PI <- numeric(0)
  result_list$bw_LIK <- numeric(0)
  result_list$bw_SCV <- numeric(0)
  
  
  result_list$bw_adaptive_PI <- array(dim=c(n,2,0))
  result_list$bw_adaptive_LIK <- array(dim=c(n,2,0))
  result_list$bw_adaptive_SCV <- array(dim=c(n,2,0))
  
  
  result_list$FIX_PI_ISE_matrix <- 1:100
  result_list$FIX_LIK_ISE_matrix <- 1:100
  result_list$FIX_SCV_ISE_matrix <- 1:100
  
  
  result_list$ADA_PI_ISE_matrix <- 1:100
  result_list$ADA_LIK_ISE_matrix <- 1:100
  result_list$ADA_SCV_ISE_matrix <- 1:100
  
  
  result_list$points_matrix <- array(dim=c(n,3,0))
  
  result_list$density_list <- Synthetic.pdf(metric_list)
  
  return(result_list)
  
}


# Based on the generalized formula, construct a specific synthetic PDF (spatiotemporal pdf).
Synthetic.pdf <- function(metric_list){
  
  #  generate synthetic pdf from the generic synthetic pdf formula
  
  xy_density <- list() #spatial density of each time slice
  
  z_density <- z_density_gaussian_mixture  <- numeric(length(metric_list$grz))  # temporal density
  
  time_layer <- 1:length(metric_list$grz)
  
  peak1 <- c(0.25, 0.04, 1/3)
  peak2 <- c(0.50, 0.02, 1/3)
  peak3 <- c(0.75, 0.03, 1/3)
  para <- c(1.0989e-9, -0.0002188)
  
  { #temporal density first
    
    time_range <- range(metric_list$grz)
    
    peaks_matrix <- rbind(peak1, peak2, peak3)
    
    apply(peaks_matrix, 1, function(row) {
      
      z_pdf <- dnorm(metric_list$grz, row[1] * diff(time_range) + time_range[1], row[2]*diff(time_range))
      
      z_pdf <- z_pdf / sum(z_pdf * metric_list$time_unit)
      
      z_density_gaussian_mixture <<- z_density_gaussian_mixture + z_pdf * row[3]
      
    })
    
    #linearly increasing density
    
    z_density_linear <- para[1] * metric_list$grz + para[2]
    
    z_density_linear <- z_density_linear / sum(z_density_linear * metric_list$time_unit)
    
    z_density <- 0.5 * z_density_gaussian_mixture + 0.5* z_density_linear
    
  } #temporal density first
  
  {
    
    boundary_x <- metric_list$boundary_x
    boundary_y <- metric_list$boundary_y
    diff_x <- diff(metric_list$boundary_x)
    diff_y <- diff(metric_list$boundary_y)
    
    xy_background <- xy_noise <- NULL
    
    {#  simulates the macroscopic state of central density and peripheral sparsity
      
      background_peaks <- c(0, 0, 1/6, 1/6) 
      x_background <- dnorm(metric_list$grx, mean(boundary_x) + background_peaks[1]*diff_x, background_peaks[3]*diff_x/2)
      y_background <- dnorm(metric_list$gry, mean(boundary_y) + background_peaks[2]*diff_y, background_peaks[4]*diff_y/2)
      xy_background <- outer(y_background, x_background, "*")
      xy_background <- xy_background / DATA.integral2d(xy_background, metric_list$xstep, metric_list$ystep)
      # cat("xy_background的维数是", dim(xy_background),"\n")
      
    }#  simulates the macroscopic state of central density and peripheral sparsity
    
    {# uniform noise background
      
      x_noise <- dunif(metric_list$grx, boundary_x[1]+3000, boundary_x[2]-3000)
      y_noise <- dunif(metric_list$gry, boundary_y[1]+3000, boundary_y[2]-3000)
      
      xy_noise <- outer(y_noise, x_noise, "*")
      xy_noise <- xy_noise / DATA.integral2d(xy_noise, metric_list$xstep, metric_list$ystep)
      
      # cat("xy_noise的维数是", dim(xy_noise))
      
    }# uniform noise background
    
    xy_density <- lapply(time_layer, function(x) {
      
      theta <- numeric(3)
      r <-     numeric(3)
      sigma <- numeric(3)
      
      xy_density_this_layer <- matrix(0, nrow = length(metric_list$gry), ncol = length(metric_list$grx))
      
      if( x <= length(metric_list$grz) * (peak1[1] + 2 * peak1[2])) {
        
        theta <- c(60, 232, 310)
        r <-     c(0.55, 0.3, 0.4)
        sigma <- c(0.02, 0.03, 0.04)  #key parameters
        
      }
      else if( x > length(metric_list$grz) * (peak1[1] + 2 * peak1[2])  &&  x <= length(metric_list$grz) * (peak2[1] + 2 * peak2[2]) ){
        
        theta <- c(6, 185, 239)
        r <-     c(0.45, 0.2, 0.6)
        sigma <- c(0.038, 0.028, 0.018) #key parameters
        
      }
      else {
        
        theta <- c(113, 236, 331)
        r <-     c(0.52, 0.15, 0.35)
        sigma <- c(0.032,0.042,0.022) #key parameters
        
      }
      
      peaks_matrix <- rbind(c(r[1]*cos(theta[1]*pi/180), r[1]*sin(theta[1]*pi/180), sigma[1], sigma[1], 1/3),
                            c(r[2]*cos(theta[2]*pi/180), r[2]*sin(theta[2]*pi/180), sigma[2], sigma[2], 1/3),
                            c(r[3]*cos(theta[3]*pi/180), r[3]*sin(theta[3]*pi/180), sigma[3], sigma[3], 1/3))
      
      apply(peaks_matrix, 1, function(row) {
        
        x_pdf <- dnorm(metric_list$grx, mean(boundary_x) + row[1] * diff_x/2, row[3] * diff_x)
        y_pdf <- dnorm(metric_list$gry, mean(boundary_y) + row[2] * diff_y/2, row[4] * diff_y)
        
        xy_pdf <- outer(y_pdf, x_pdf, "*")
        
        xy_density_this_layer <<- xy_density_this_layer + (xy_pdf / DATA.integral2d(xy_pdf, metric_list$xstep, metric_list$ystep) * row[5])
        
      })
      
      xy_density_this_layer <- (0.5*xy_background + 0.5*xy_density_this_layer) *0.95 + xy_noise * 0.05
      
      return(xy_density_this_layer)
      
    })
    
  } #spatial density second
  
  return(list(xy_density = xy_density, z_density = z_density))
  
}


# Start this experiment using parallel computing.
Synthetic.parallel <- function(origin_list, metric_list, m) {
  
  cat("begin...... ",format(Sys.time(),"%H:%M:%S"),"\n")
  
  start_index <- ifelse(length(origin_list$iteration) < 1, 1, max(origin_list$iteration) + 1)
  
  cat("iterations from ",start_index,"to start...\n")
  
  time_array <- numeric(m)
  
  density_list <- origin_list$density_list
  
  if(is.null(density_list)) {
    
    stop("density_list is null!")
    
  }
  
  n <- origin_list$n
  
  raw_density <- array(0,dim = metric_list$dimensionxyz[c(2,1,3)]) # y, x, z
  
  if(is.list(density_list$xy_density)) {
    
    for(i in seq_along(density_list$xy_density)){
      
      raw_density[,,i] <- density_list$xy_density[[i]] * density_list$z_density[i]
      
    }
    
  }
  else {
    
    raw_density <- outer(density_list$xy_density, density_list$z_density, "*")
    
  }
  
  for(j in 1:m) {
    
    gc()
    
    last_start_time <- Sys.time()
    
    cl <- makeCluster(2)
    
    clusterExport(cl, varlist = c("metric_list","Synthetic.STKDE.PI.LIK.SCV","STKDE.q.pointsmatrix.to.fhatvector",
                                  "STKDE.constructor.with.bandwidths",
                                  "DATA.integral3d",
                                  "Hpi.diag","LIK.spattemp","OS.spattemp","Hscv.diag","ppp","owin","makeCluster","clusterExport","parLapply","stopCluster",
                                  "order_2"))
    
    list_parallel_A <- Synthetic.points(n, metric_list, density_list)
    list_parallel_B <- Synthetic.points(n, metric_list, density_list)
    
    synthetic_sub_list <- parLapply(cl,list(list_parallel_A, list_parallel_B), 
                                    
                                    function(matrix){
                                      
                                      Synthetic.STKDE.PI.LIK.SCV(matrix, metric_list)
                                      
                                    })
    stopCluster(cl)
    
    for(i in seq_along(synthetic_sub_list)){
      
      if(is.null(synthetic_sub_list[[i]])){
        
        end_time <- Sys.time()
        
        cat("The ",j,"-th in ",m," times ignored*. Finish time: ",format(end_time,"%H:%M:%S"),
            ". Time cost: ", difftime(end_time,last_start_time,units="mins"),"\n")
        
        next
        
      }
      
      variable_name <- paste0("list_parallel_", LETTERS[i])
      
      origin_list$points_matrix <- abind(origin_list$points_matrix, get(variable_name), along = 3)
      
      #========================# record the fixed and adaptive bandwidths ======================
      
      origin_list$bw_PI <- cbind(origin_list$bw_PI, synthetic_sub_list[[i]]$bw_PI)
      origin_list$bw_LIK <- cbind(origin_list$bw_LIK, synthetic_sub_list[[i]]$bw_LIK)
      origin_list$bw_SCV <- cbind(origin_list$bw_SCV, synthetic_sub_list[[i]]$bw_SCV)
      
      origin_list$bw_adaptive_PI <- abind(origin_list$bw_adaptive_PI, synthetic_sub_list[[i]]$bw_adaptive_PI, along = 3)
      origin_list$bw_adaptive_LIK <- abind(origin_list$bw_adaptive_LIK, synthetic_sub_list[[i]]$bw_adaptive_LIK, along = 3)
      origin_list$bw_adaptive_SCV <- abind(origin_list$bw_adaptive_SCV, synthetic_sub_list[[i]]$bw_adaptive_SCV, along = 3)
      
      #=========================#计算局部ISE=======================
      
      cl <- makeCluster(6)
      
      clusterExport(cl, varlist = c("Synthetic.local.ISE"))
      
      local_ISE_sub_list <- parLapply(cl,list(synthetic_sub_list[[i]]$FIX_PI_STKDE,  # 1 - FIX_PI
                                              synthetic_sub_list[[i]]$FIX_LIK_STKDE, # 2 - FIX_LIK
                                              synthetic_sub_list[[i]]$FIX_SCV_STKDE, # 3 - FIX_SCV
                                              
                                              synthetic_sub_list[[i]]$ADA_PI_STKDE,  # 4 - ADA-PI
                                              synthetic_sub_list[[i]]$ADA_LIK_STKDE, # 5 - ADA-LIK
                                              synthetic_sub_list[[i]]$ADA_SCV_STKDE),# 6 - ADA_SCV
                                      
                                      function(matrix){
                                        
                                        Synthetic.local.ISE(matrix, raw_density, metric_list$xstep, metric_list$ystep, metric_list$time_unit)
                                        
                                      })
      stopCluster(cl)
      
      origin_list$FIX_PI_ISE_matrix <- rbind(origin_list$FIX_PI_ISE_matrix, local_ISE_sub_list[[1]])
      origin_list$FIX_LIK_ISE_matrix <- rbind(origin_list$FIX_LIK_ISE_matrix, local_ISE_sub_list[[2]])
      origin_list$FIX_SCV_ISE_matrix <- rbind(origin_list$FIX_SCV_ISE_matrix, local_ISE_sub_list[[3]])
      
      origin_list$ADA_PI_ISE_matrix <- rbind(origin_list$ADA_PI_ISE_matrix, local_ISE_sub_list[[4]])
      origin_list$ADA_LIK_ISE_matrix <- rbind(origin_list$ADA_LIK_ISE_matrix, local_ISE_sub_list[[5]])
      origin_list$ADA_SCV_ISE_matrix <- rbind(origin_list$ADA_SCV_ISE_matrix, local_ISE_sub_list[[6]])
      
      
      #========================= iteration = iteration + 1 =========================
      
      origin_list$iteration <- c(origin_list$iteration,start_index + 2 * (j-1) + (i - 1))
      
    } 
    
    end_time <- Sys.time()
    
    time_array[j] <- difftime(end_time,last_start_time,units="mins")
    
    cat("The ",j,"-th in ",m," times complete. Finish time: ",format(end_time,"%H:%M:%S"),
        ". Time cost: ", time_array[j],"\n")
    
    
  }
  
  cat("finish...... ",format(Sys.time(),"%H:%M:%S"),".  ")
  cat("Total :",sum(time_array),". Avarage:",mean(time_array),"\n")
  return(origin_list)
  
}


#============================== SUB FUNCTIONS =================================


# Based on the input pdf, execute a stochastic process to generate n points (3D).
Synthetic.points <- function(n,metric_list,density_list){
  
  # stochastic process, generate n points from the given pdf
  
  xy_points <- z_points <- NULL
  
  if(is.matrix(density_list$xy_density)) {
    
    xy_plain_probability <- density_list$xy_density * xstep * ystep
    
    im_obj <- im(xy_plain_probability,xcol = metric_list$grx, yrow = metric_list$gry)
    
    xy_points <- rpoint(n,f = im_obj)
    
    z_probability <- density_list$z_density * metric_list$time_unit
    
    z_probability_matrix <- as.matrix(cbind(rep(0,length(metric_list$grz)), z_probability))
    
    im_obj <- im(z_probability_matrix, xcol = metric_list$grx[c(1,2)], yrow = metric_list$grz)
    
    z_points <- rpoint(n, f = im_obj)
    
    t_points_matrix <- cbind(xy_points$x,xy_points$y,z_points$y)
    
    return(t_points_matrix)
    
  }
  
  else if(is.list(density_list$xy_density)){
    
    z_probability <- density_list$z_density * metric_list$time_unit
    
    z_probability_matrix <- cbind(rep(0,length(z_probability)), z_probability)
    
    im_obj <- im(z_probability_matrix, xcol = metric_list$grx[c(1,2)], yrow = metric_list$grz)
    
    #1.temporal tochastic process
    
    z_points <- rpoint(n, f = im_obj)
    
    z_points_vector <- z_points$y
    
    z_index <- sapply(z_points_vector, function(x) {
      
      return(which.min(abs(x - metric_list$grz)))
      
    })
    
    z_points_table <- table(z_index)
    
    z_points_matrix <- cbind(as.integer(rownames(z_points_table)), z_points_table)
    
    #2.spatial tochastic process
    
    t_points_matrix <- do.call(rbind, apply(z_points_matrix, 1, function(row) {
      
      xy_image <- im(density_list$xy_density[[row[1]]], xcol = metric_list$grx, yrow = metric_list$gry)
      
      xy_points <- rpoint(row[2], f = xy_image)
      
      return(cbind(xy_points$x, xy_points$y, z_points_vector[which(z_index == row[1])]))
      
    }))
    
    return(t_points_matrix)
    
  }
  
}


# For the training set, compute TSTKDE and ASTKDE using the plug-in, LIK, and SCV methods separately.
Synthetic.STKDE.PI.LIK.SCV <- function(t_points_matrix,metric_list) {
  
  # INPUT observation points matrix  and metric_list
  
  # OUTPUT global and top n% voxels ISE

  n <- nrow(t_points_matrix)
  
  grx <- metric_list$grx
  gry <- metric_list$gry
  grz <- metric_list$grz
  
  xres <- length(grx)
  yres <- length(gry)
  zres <- length(grz)
  
  xstep <- metric_list$xstep
  ystep <- metric_list$ystep
  time_unit <- metric_list$time_unit
  
  t_points_matrix_A <- t_points_matrix_B <- t_points_matrix_C <- NULL
  
  bw_PI <-  bw_adaptive_PI <- NULL
  bw_LIK <- bw_adaptive_LIK <- NULL
  bw_SCV <- bw_adaptive_SCV <- NULL
  
  q_PI <- q_LIK <- q_LCV <- NULL
  
  FIX_PI_STKDE <- ADA_PI_STKDE <- NULL
  FIX_LIK_STKDE <- ADA_LIK_STKDE <- NULL
  FIX_SCV_STKDE <- ADA_SCV_STKDE <- NULL 
  
  hp_PI <- lambda_PI <- NULL
  hp_LIK <- lambda_LIK <- NULL
  hp_SCV <- lambda_SCV <- NULL
  
  # for fixed bandwidths
  
  {
    
    bw_PI <- sqrt(diag((Hpi.diag(t_points_matrix,nstage = 2,pilot="dscalar",verbose=FALSE,
                                 optim.fun="optim",binned = FALSE))))
    hp_PI <- mean(bw_PI[c(1,2)])
    lambda_PI<- bw_PI[3]
    
    point_pattern <- ppp(x=t_points_matrix[,1],y=t_points_matrix[,2],
                         window = owin(xrange = metric_list$boundary_x, yrange=metric_list$boundary_y))
    
    bw_LIK <- LIK.spattemp(point_pattern, tt = t_points_matrix[,3], sedge="uniform", verbose = FALSE)
    hp_LIK <- bw_LIK[1]
    lambda_LIK <- bw_LIK[2]
    
    bw_SCV <- sqrt(diag(Hscv.diag(t_points_matrix, nstage=2, pilot = "dscalar",verbose=FALSE, optim.fun="optim", binned = FALSE)))
    hp_SCV <- mean(bw_SCV[c(1,2)])
    lambda_SCV <- bw_SCV[3]
    
  }
  
  # leave-one cross validation is really unstable when inhomogeneity increases, need to defend here
  
  if(hp_SCV < 50 || lambda_SCV > 3000) {
    
    return(NULL)
    
  }
  
  if(lambda_LIK > 2500) {
    
    return(NULL)
    
  }
  
  t_points_matrix_A <- cbind(t_points_matrix,rep(hp_PI,n),rep(lambda_PI,n))
  t_points_matrix_B <- cbind(t_points_matrix,rep(hp_LIK,n),rep(lambda_LIK,n))
  t_points_matrix_C <- cbind(t_points_matrix,rep(hp_SCV,n),rep(lambda_SCV,n))
  
  # parallel FIX STKDE of PI & LIK $ SCV 
  
  {
    
    cl2 <- makeCluster(3)
    
    clusterExport(cl2, varlist = c("STKDE.constructor.with.bandwidths","DATA.integral3d","order_2","metric_list"))
    
    STKDE_list <- parLapply(cl2,list(t_points_matrix_A, t_points_matrix_B, t_points_matrix_C), 
                            
                            function(matrix){
                              
                              STKDE.constructor.with.bandwidths(matrix, metric_list)
                              
                            })
    stopCluster(cl2)
    
    FIX_PI_STKDE <- STKDE_list[[1]]$result_STKDE
    q_PI <- STKDE_list[[1]]$q
    
    FIX_LIK_STKDE <- STKDE_list[[2]]$result_STKDE
    q_LIK <- STKDE_list[[2]]$q
    
    FIX_SCV_STKDE <- STKDE_list[[3]]$result_STKDE
    q_SCV <- STKDE_list[[3]]$q
    
  }
  
  cat("F ")
  
  # adaptive bandwidths
  
  {
    # for PI adaptive
    
    fhat_vector <- STKDE.q.pointsmatrix.to.fhatvector(t_points_matrix_A,t_points_matrix_A,q_PI)
    
    gamma <- exp(mean(log(fhat_vector^order_2)))
    h_adaptive <- hp_PI * (fhat_vector^order_2)/gamma
    lambda_adaptive <- lambda_PI * (fhat_vector^order_2)/gamma
    
    bw_adaptive_PI <- cbind(h_adaptive, lambda_adaptive)
    
    t_points_matrix_A <- cbind(t_points_matrix,h_adaptive,lambda_adaptive)
    
    # for LIK adaptive
    
    fhat_vector <- STKDE.q.pointsmatrix.to.fhatvector(t_points_matrix_B,t_points_matrix_B,q_LIK)
    
    gamma <- exp(mean(log(fhat_vector^order_2)))
    h_adaptive <- hp_LIK* (fhat_vector^order_2)/gamma
    lambda_adaptive <- lambda_LIK * (fhat_vector^order_2)/gamma
    
    bw_adaptive_LIK <- cbind(h_adaptive, lambda_adaptive)
    
    t_points_matrix_B <- cbind(t_points_matrix,h_adaptive,lambda_adaptive)
    
    # for SCV ADAPTIVE
    
    fhat_vector <- STKDE.q.pointsmatrix.to.fhatvector(t_points_matrix_C,t_points_matrix_C,q_SCV)
    
    gamma <- exp(mean(log(fhat_vector^order_2)))
    h_adaptive <- hp_SCV* (fhat_vector^order_2)/gamma
    lambda_adaptive <- lambda_SCV * (fhat_vector^order_2)/gamma
    
    bw_adaptive_SCV <- cbind(h_adaptive, lambda_adaptive)
    
    t_points_matrix_C <- cbind(t_points_matrix,h_adaptive,lambda_adaptive)

    
  }
  
  # parallel ASTKDE with PI & LIK & SCV
  
  {
    
    cl2 <- makeCluster(3)
    
    clusterExport(cl2, varlist = c("STKDE.constructor.with.bandwidths","DATA.integral3d","order_2","metric_list"))
    
    STKDE_list <- parLapply(cl2,list(t_points_matrix_A, t_points_matrix_B, t_points_matrix_C), 
                            
                            function(matrix){
                              
                              STKDE.constructor.with.bandwidths(matrix,metric_list)
                              
                            })
    stopCluster(cl2)
    
    ADA_PI_STKDE <- STKDE_list[[1]]$result_STKDE
    ADA_LIK_STKDE <- STKDE_list[[2]]$result_STKDE
    ADA_SCV_STKDE <- STKDE_list[[3]]$result_STKDE
    
  }
  
  cat("A ")
  
  return(list(FIX_PI_STKDE = FIX_PI_STKDE,
              FIX_LIK_STKDE = FIX_LIK_STKDE,
              FIX_SCV_STKDE = FIX_SCV_STKDE,
              ADA_PI_STKDE = ADA_PI_STKDE,
              ADA_LIK_STKDE = ADA_LIK_STKDE,
              ADA_SCV_STKDE = ADA_SCV_STKDE,
              
              bw_PI = bw_PI, bw_LIK = bw_LIK,bw_SCV = bw_SCV,
              bw_adaptive_PI = bw_adaptive_PI, bw_adaptive_LIK = bw_adaptive_LIK,
              bw_adaptive_SCV = bw_adaptive_SCV))
  
}


# Compute the ISE values for voxels ranging from the top 1% to the top 100%.
Synthetic.local.ISE <- function(STKDE, raw_STKDE, xstep, ystep, time_unit) {
  
  q <- 1:100
  
  q <- quantile(STKDE, (100-q)/100)
  
  local_ISE_vector <- sapply(1:100, function(x) {
    
    A <- STKDE[STKDE >= q[x]]
    
    B <- raw_STKDE[STKDE >= q[x]]
    
    sum((A-B) ^ 2) * xstep * ystep * time_unit
    
  })
  
  return(local_ISE_vector)
  
}

#=============================== RUN CODE HERE ===================================

synthetic_result_list <- Synthetic.result.list.constructor(metric_list,500)

synthetic_result_list <- Synthetic.parallel(synthetic_result_list, metric_list, 5)
