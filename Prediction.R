
#============================== MAIN FUNCTIONS =================================

# Compute the log-likelihood for each iteration (i.e., each prediction).
Prediction.LIK.iteration <- function(list,ALLdataframe,metric_list, parallel = FALSE){ #list包含 TTR  VTR  strategy
  
  # INPUT  : list, include TTR, VTR & name
  #          ALLdataframe, include all training points and test points
  #          metric_list_full is needed in this kind of compute
  
  # OUTPUT : prediction likelihood of TSTKDE & ASTKDE in the specific iteration

  start_time <- Sys.time()
  
  time_unit <- metric_list$time_unit

  time_range <- range(ALLdataframe$DATE)
  
  grz <- NULL
  
  if(metric_list$type == "full"){
    
    grz <- seq(time_range[1]+time_unit/2,227175,by=time_unit)  # if dim = 695,735,365,full
    
  }
  else{
    
    stop("metric_list$type needs to be 'full'")
    
  }

  TTR <- list$TTR
  VTR <- list$VTR
  
  Tdataframe <- subset(ALLdataframe,DATE >= TTR[1] & DATE < TTR[2],select=c(X,Y,DATE))
  Vdataframe <- subset(ALLdataframe,DATE >= VTR[1] & DATE < VTR[2],select=c(X,Y,DATE))
  
  tn <- nrow(Tdataframe)  # number of points in training set
  
  vn <- nrow(Vdataframe)  # number of points in validation set (test set)
  
  #PI bandwidth
    
  matrix <- matrix(c(Tdataframe$X,Tdataframe$Y,Tdataframe$DATE),nrow=nrow(Tdataframe),ncol=3) #XYT矩阵
  pilot_bw_vector <- ceiling(sqrt(diag((Hpi.diag(matrix,nstage = 2,pilot="dscalar",verbose=FALSE,
                                                 optim.fun="optim",binned = FALSE)))))
  
  h <- mean(pilot_bw_vector[c(1,2)])
  lambda <- pilot_bw_vector[3]
  
  cat("h,lambda :", h ," ---- ", lambda,"\n")
  
  reference_grz <- abs(grz - TTR[1])
  T_index_start <- which.min(reference_grz)
  
  reference_grz <- abs(grz - TTR[2])
  T_index_end <- which.min(reference_grz)
  
  reference_grz <- abs(grz - VTR[2])
  V_index_end <- which.min(reference_grz)
  
  grz_this <- grz[T_index_start:T_index_end]
  grz_extend_this <- grz[T_index_start:V_index_end]
  
  cat("grz_this        :",grz_this[1],"-----",grz_this[length(grz_this)],",",length(grz_this),"\n")
  cat("grz_extend_this :",grz_extend_this[1],"-----",grz_extend_this[length(grz_extend_this)],",",length(grz_extend_this),"\n")
  
  # for TSTKDE :
  
  t_points_matrix <- matrix(c(Tdataframe$X,Tdataframe$Y,Tdataframe$DATE,rep(h,tn),rep(lambda,tn)),nrow=tn,ncol=5)  #1257 * 5 matrix
  v_points_matrix <- matrix(c(Vdataframe$X,Vdataframe$Y,Vdataframe$DATE), nrow=vn,ncol=3)
  
  TSTKDE_result <- Prediction.on.testset.fhatvector(grz_extend_this,t_points_matrix,v_points_matrix, metric_list, parallel)
  TSTKDE_LIK <- sum(log(TSTKDE_result$fhat_vector))
  
  cat("TSTKDE_LIK 计算完毕!", difftime(Sys.time(),start_time,units="mins"),"\n")
  cat("TSTKDE_LIK = ",TSTKDE_LIK,"\n")
  start_time <- Sys.time()
  
  # for ASTKDE :
  
  TSTKDE_result <- Prediction.on.testset.fhatvector(grz_this, t_points_matrix, t_points_matrix, metric_list, parallel)
  
  fhat_vector <- TSTKDE_result$fhat_vector
  
  cat("fhat_vector = ", mean(fhat_vector),"\n")
  
  gamma <- exp(mean(log(fhat_vector^order_2)))
  
  h_adaptive <- h * (fhat_vector^order_2)/gamma
  
  lambda_adaptive <- lambda * (fhat_vector^order_2)/gamma
  
  cat("h_adaptive:", mean(h_adaptive))
  cat("lambda_adaptive:", mean(lambda_adaptive))

  t_points_matrix <- matrix(c(Tdataframe$X,Tdataframe$Y,Tdataframe$DATE,h_adaptive,lambda_adaptive),nrow=tn,ncol=5)
  v_points_matrix <- matrix(c(Vdataframe$X,Vdataframe$Y,Vdataframe$DATE), nrow=vn,ncol=3)
  
  ASTKDE_result <- Prediction.on.testset.fhatvector(grz_extend_this,t_points_matrix,v_points_matrix, metric_list, parallel)
  ASTKDE_LIK <- sum(log(ASTKDE_result$fhat_vector))
  
  cat("ASTKDE_LIK 计算完毕!", difftime(Sys.time(),start_time,units="mins"),"\n")
  cat("ASTKDE_LIK = ",ASTKDE_LIK,"\n")

  return(list(TSTKDE_LIK=TSTKDE_LIK, ASTKDE_LIK=ASTKDE_LIK,
              name = list$name))
}

#============================== SUB FUNCTIONS =================================


# Provide the fhat values derived from the training set for each point in the test set.
Prediction.on.testset.fhatvector <- function(grz_this,t_points_matrix,v_points_matrix, metric_list, parallel){
  
  #  INPUT : grz_this, a vector of temporal domain.
  #          t_points_matrix, matrix of training set with bandwidths.
  #          v_points_matrix  XYZ matrix of test set.
  #          "'t_points_matrix' can be identical to 'v_points_matrix', which is quite ingenious".
  #          parallel, TRUE ----> "Two-thread parallel computing".
  
  
  # OUTPUT : q, edge correction factor, a vector length = nrow(t_points_matrix).
  #          fhat_vector,  a vector, which length = nrow(v_points_matrix), Stored the fhat of each point in test set.
  
  tn <- nrow(t_points_matrix)
  vn <- nrow(v_points_matrix)

  
  # 1. need to get q vector first
  
  q <- NULL
  
  if(!parallel) {
    
    q <- STKDE.pointsmatrix.to.q(grz_this, metric_list, t_points_matrix)
    
  }
  
  else{
    
    mid <- floor(tn/2)
    t_points_matrix_A <- t_points_matrix[1:mid, , drop = FALSE]
    t_points_matrix_B <- t_points_matrix[(mid+1):tn, , drop = FALSE]
    
    cl <- makeCluster(2)
    
    clusterExport(cl, varlist = c("STKDE.pointsmatrix.to.q",
                                  "DATA.integral3d",
                                  "makeCluster","clusterExport","parLapply","stopCluster"))
    
    q_list <- parLapply(cl,list(t_points_matrix_A, t_points_matrix_B), 
                             
                             function(matrix){
                               
                               STKDE.pointsmatrix.to.q(grz_this, metric_list, matrix)
                               
                             })
    stopCluster(cl)
    
    q <- unlist(q_list)
    
  }
  
  # 2. generate matrix vn * tn , for all fhat at points in the validation set (test set)
  
  matrix_x <- apply(t_points_matrix,1,function(row) dnorm(v_points_matrix[,1],row[1],row[4]))
  matrix_y <- apply(t_points_matrix,1,function(row) dnorm(v_points_matrix[,2],row[2],row[4]))
  matrix_z <- apply(t_points_matrix,1,function(row) dnorm(v_points_matrix[,3],row[3],row[5]))
  
  mx.my.mz <- matrix_x * matrix_y * matrix_z
  
  mx.my.mz <- t(apply(mx.my.mz,1,function(row) row/q))
  
  fhat_vector <- apply(mx.my.mz,1,function(row) {
    
    sum(row)/tn
    
  })
  
  return(list(q=q,fhat_vector= fhat_vector))
  
}


#=============================== RUN CODE HERE ===================================
Prediction_11 <- Prediction.LIK.iteration(list(name = "Prediction 11 : 2014/11/1 - 2024/3/1",
                                                   TTR = c(41944 * 5, 45231 * 5),
                                                   VTR = c(45231 * 5, 45352 * 5)),ALLdataframe, metric_list_full,
                                               parallel = TRUE)


#     NAME                  TTR                     VTR


# Prediction 1 : c(41640 * 5, 44927 * 5)    c(44927 * 5, 45047 * 5)
# Prediction 2 : c(41671 * 5, 44958 * 5)    c(44958 * 5, 45078 * 5)
# Prediction 3 : c(41699 * 5, 44986 * 5)    c(44986 * 5, 45108 * 5)
# Prediction 4 : c(41730 * 5,	45017 * 5)    c(45017 * 5, 45139 * 5)
# Prediction 5 : c(41760 * 5,	45047 * 5)    c(45047 * 5, 45170 * 5)
# Prediction 6 : c(41791 * 5,	45078 * 5)    c(45078 * 5, 45200 * 5)
# Prediction 7 : c(41821 * 5,	45108 * 5)    c(45108 * 5, 45231 * 5)
# Prediction 8 : c(41852 * 5,	45139 * 5)    c(45139 * 5, 45261 * 5)
# Prediction 9 : c(41883 * 5,	45170 * 5)	  c(45170 * 5, 45292 * 5)
# Prediction 10: c(41913 * 5,	45200 * 5)    c(45200 * 5, 45323 * 5)
# Prediction 11: c(41944 * 5,	45231 * 5)	  c(45231 * 5, 45352 * 5) 
# Prediction 12: c(41974 * 5, 45261 * 5)	  c(45261 * 5, 45383 * 5) 
# Prediction 13: c(42005 * 5,	45292 * 5)  	c(45292 * 5, 45413 * 5)



