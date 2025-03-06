

MC.parallel <- function(result_list,m, metric_list, TSTKDE_origin, ASTKDE_origin) {
  
  # "Conduct 999 Monte Carlo simulations using the parallel method and obtain the hotspot identification results 
  #  of TSTKDE and ASTKDE based on Equation (9)."
  
  # "Users can dynamically reduce or increase the number of parallel threads based on their PC's performance."
  
  cat("Begin...... ",format(Sys.time(),"%H:%M:%S"),"\n")
  
  start_index <- ifelse(length(result_list$iteration) < 1, 1, max(result_list$iteration) + 1)

  cat("iterations from ",start_index,"to start...\n")
  
  time_array <- numeric(m)
  
  all_voxels <- prod(metric_list$dimensionxyz)
  
  for(j in 1:m) {
    
    gc()
    
    last_start_time <- Sys.time()
    
    cl <- makeCluster(4)
    
    clusterExport(cl, varlist = c("metric_list","MC.compute","MC.random.points",
                                  "STKDE.joint.constructor.with.bandwidths","STKDE.q.pointsmatrix.to.fhatvector",
                                  "DATA.integral2d","DATA.integral3d",
                                  "Hpi.diag","makeCluster","clusterExport","parLapply","stopCluster",
                                  "TSTKDE_joint","ASTKDE_joint",
                                  "order_2"))
    
    list_parallel_A <- list(points_3d_matrix = MC.random.points(1257,metric_list), iteration = start_index + 1 * (j-1) + 0)
    list_parallel_B <- list(points_3d_matrix = MC.random.points(1257,metric_list), iteration = start_index + 2 * (j-1) + 1)
    list_parallel_C <- list(points_3d_matrix = MC.random.points(1257,metric_list), iteration = start_index + 3 * (j-1) + 2)
    list_parallel_D <- list(points_3d_matrix = MC.random.points(1257,metric_list), iteration = start_index + 4 * (j-1) + 3)
    
    MC_sub_list <- parLapply(cl,list(list_parallel_A, list_parallel_B, list_parallel_C, list_parallel_D), 
                             
                             function(list){
                               
                               MC.compute(metric_list, list, TSTKDE_origin, ASTKDE_origin)
                               
                             })
    stopCluster(cl)
    
    for(i in seq_along(MC_sub_list)){
      
      if(MC_sub_list[[i]]$iteration >= 1000) next  #defend code
      
      result_list$TSTKDE_result_matrix <- result_list$TSTKDE_result_matrix + MC_sub_list[[i]]$TSTKDE_indicator_array
      result_list$ASTKDE_result_matrix <- result_list$ASTKDE_result_matrix + MC_sub_list[[i]]$ASTKDE_indicator_array
    
      result_list$TSTKDE_SS_vector <- c(result_list$TSTKDE_SS_vector,sum(result_list$TSTKDE_result_matrix == 0,na.rm=TRUE)/all_voxels)
      result_list$ASTKDE_SS_vector <- c(result_list$ASTKDE_SS_vector,sum(result_list$ASTKDE_result_matrix == 0,na.rm=TRUE)/all_voxels)
      
      result_list$iteration <- c(result_list$iteration,MC_sub_list[[i]]$iteration)
      
    }
    
    end_time <- Sys.time()
    
    time_array[j] <- difftime(end_time,last_start_time,units="mins")
    
    cat("The ",j,"-th in ",m," times complete. Finish time: ",format(end_time,"%H:%M:%S"),
        ". Time cost: ", time_array[j],"\n")
    
    
  }
  
  cat("Finished...... ",format(Sys.time(),"%H:%M:%S"),".  ")
  cat("Total :",sum(time_array),". Avarage:",mean(time_array),"\n")
  return(result_list)
  
}

MC.random.points <- function(n, metric_list){
  
  # Uniformly random n points in valid spatiotemporal domain
  
  # INPUT  : n, the same points number in your real data.
  # OUTPUT : points_3d, matrix n * 3
  
  time_range <- range(metric_list$grz)
  
  points_2d <- runifpoint(n,win=metric_list$obwindow)
  #obwindow is generated in DATA.R
  #obwindow is either the study area (an irregular closed polygon) or a regular rectangle.
  
  time_1d <- sample(time_range[1]:time_range[2],n,replace=TRUE)
  
  points_3d <- matrix(c(points_2d$x,points_2d$y,time_1d),nrow=n,ncol=3)
  
  return(points_3d)
  
}

MC.compute <- function(metric_list, list, TSTKDE_origin, ASTKDE_origin ){ 
  
  #输入为蒙特卡洛随机点，输出为减去standard数组后的01数组
  
  cat("Begin......",format(Sys.time(),"%H:%M:%S"),"\n")
  
  points_3d_matrix <- list$points_3d_matrix
  
  n <- nrow(points_3d_matrix)
  
  grx <- metric_list$grx
  gry <- metric_list$gry
  grz <- metric_list$grz
  
  xstep <- metric_list$xstep
  ystep <- metric_list$ystep
  time_unit <- metric_list$time_unit
  
  xres <- metric_list$dimensionxyz[1]
  yres <- metric_list$dimensionxyz[2]
  zres <- metric_list$dimensionxyz[3]
  
  numeric_01_cube <- metric_list$numeric_01_cube
  
  # fixed bandwidth STKDE :
  
  pilot_bw_vector <- ceiling(sqrt(diag((Hpi.diag(points_3d_matrix,nstage = 2,pilot="dscalar",verbose=FALSE,
                                                 optim.fun="optim",binned = FALSE)))))
  
  hp <- mean(pilot_bw_vector[c(1,2)])
  lambda <- pilot_bw_vector[3]
  
  t_points_matrix_A <- cbind(points_3d_matrix,rep(hp,n),rep(lambda,n))
  
  TSTKDE_result <- STKDE.joint.constructor.with.bandwidths(t_points_matrix_A, metric_list)
  
  TSTKDE_star <- TSTKDE_result$result_STKDE
  
  q <- TSTKDE_result$q
  
  # adaptive bandwidth STKDE :
  
  fhat_vector <- STKDE.q.pointsmatrix.to.fhatvector(t_points_matrix_A, t_points_matrix_A, q)
  
  gamma <- exp(mean(log(fhat_vector^order_2)))
  h_adaptive <- hp * (fhat_vector^order_2)/gamma
  lambda_adaptive <- lambda * (fhat_vector^order_2)/gamma
  
  t_points_matrix_B <- cbind(points_3d_matrix,h_adaptive,lambda_adaptive)
  
  ASTKDE_result <- STKDE.joint.constructor.with.bandwidths(t_points_matrix_B, metric_list)
  
  ASTKDE_star <- ASTKDE_result$result_STKDE
  
  # "The secondary evaluation of TSTKDE and ASTKDE is based on Equation (9)."
  
  # TSTKDE_star[TSTKDE_star == 0] <- NA
  # ASTKDE_star[ASTKDE_star == 0] <- NA
  
  cat("TSTKDE_star & ASTKDE_star are done!",format(Sys.time(),"%H:%M:$S"),"\n")
  
  TSTKDE_indicator_array <- ASTKDE_indicator_array <- array(0,c(yres,xres,zres))
  
  {
    #compared with TSTKDE, OUTPUT : TSTKDE_indicator_array, see Equation (9)
    
    TSTKDE_indicator_array <- TSTKDE_star - TSTKDE_origin
    
    TSTKDE_indicator_array[ TSTKDE_indicator_array > 0 ] <- 1
    TSTKDE_indicator_array[ TSTKDE_indicator_array < 0 ] <- 0

    
  } #compared with TSTKDE, OUTPUT : TSTKDE_indicator_array,  see Equation (9)
  
  { #compared with ASTKDE, OUTPUT : ASTKDE_indicator_array, see Equation (9)
    
    ASTKDE_indicator_array <- ASTKDE_star - ASTKDE_origin
    
    ASTKDE_indicator_array[ ASTKDE_indicator_array > 0 ] <- 1
    ASTKDE_indicator_array[ ASTKDE_indicator_array < 0 ] <- 0
    
    
  } #compared with ASTKDE, OUTPUT : ASTKDE_indicator_array, see Equation (9)
  
  return(list(TSTKDE_indicator_array = TSTKDE_indicator_array,
              ASTKDE_indicator_array = ASTKDE_indicator_array,
              iteration = list$iteration))
  
} 

MC.result.constructor <- function(metric_list){
  
  init_matrix <- outer(metric_list$numeric_01_matrix, rep(1,metric_list$dimensionxyz[3]), "*")
  
  init_matrix[init_matrix==0] <- NA
  
  init_matrix[init_matrix==1] <- 0
  
  TSTKDE_result_matrix <- init_matrix
  ASTKDE_result_matrix <- init_matrix
  
  return(list(TSTKDE_result_matrix = TSTKDE_result_matrix,
              ASTKDE_result_matrix = ASTKDE_result_matrix,
              TSTKDE_SS_vector = NULL,
              ASTKDE_SS_vector = NULL,
              iteration = NULL))
  
  
}



#=======================RUN===========================

# 1st step : MC_result_list need to be generated
MC_result_list <- MC.result.constructor(metric_list)

# 2nd step : run MC.parallel to begin the MC simulation, make sure TSTKDE_joint & ASTKDE_joint are in the Global Environment
MC_result_list <- MC.parallel(MC_result_list, 5, metric_list, TSTKDE_joint, ASTKDE_joint)

