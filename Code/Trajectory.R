

#============================== MAIN FUNCTIONS =================================


# Identify all density peaks from the given STKDE, applicable to both the joint model and the conditional model.
# generate the trajectory_list.
Trajectory.Init <- function(STKDE, metric_list) {
  
  # INPUT  : STKDE, either joint or cond.
  # OUTPUT : trajectory_list, which includes all the trajectories.
  
  trajectory_list <- list()
  layer_list <- list()
  
  n <- dim(STKDE)[3]
  
  distance_threshold <- 200   # "Peaks smaller than this threshold will be combined into a trajectory."
  
  if(n == 183) distance_threshold <- 500 # "Peaks smaller than this threshold will be combined into a trajectory."
  
  kernel <- matrix(1, nrow=5, ncol=5)
  
  # 1. compute the peaks layer by layer
  
  sapply(1:n , function(z) {
    
    slice_raster <- raster(STKDE[,,z])
    
    local_max_raster <- focal(slice_raster, kernel, function(x) {
      
      center_value <- x[ceiling(length(x)/2)]
      
      if(is.na(center_value)){
        
        return(NA)
        
      }
      else if(center_value == max(x, na.rm = TRUE)){
        
        return(center_value)
        
      }
      else {
        
        return(NA)
        
      }
      
    }, pad = TRUE, padValue = NA)
    
    local_max_vector <- getValues(local_max_raster)
    local_max_index <- which(!is.na(local_max_vector))
    local_max_index_matrix <- rowColFromCell(slice_raster, local_max_index)
    local_max_index_matrix <- cbind(local_max_index_matrix,z) # with z coord.
    
    centroids_coords_matrix <- cbind(metric_list$grx[local_max_index_matrix[,2]],
                                     metric_list$gry[local_max_index_matrix[,1]],
                                     metric_list$grz[local_max_index_matrix[,3]])
    
    # 2.  link the peaks into "trajectory_list"
    
    apply(centroids_coords_matrix,1 , function(row) {
      
      centroid_3d <- st_point(c(row[1],row[2],row[3]))
      
      i <- z
      
      assigned.flag <- FALSE
      
      if(length(trajectory_list) == 0 ){
        
        trajectory_new <- list(centroid_3d)
        trajectory_list <<- append(trajectory_list,list(trajectory_new))
        
      }
      
      if(length(layer_list)==0){
        
        layers_new <- list(i)
        layer_list <<- append(layer_list,layers_new)
        
      }
      
      for( j in seq_along(trajectory_list)){
        
        tail_point <- tail(trajectory_list[[j]],1)[[1]]
        
        distance <- sqrt((tail_point[1]-centroid_3d[1])^2 + 
                           (tail_point[2]-centroid_3d[2])^2 + 
                           (tail_point[3]-centroid_3d[3])^2)
        
        if(distance < distance_threshold) {

          trajectory_list[[j]] <<- append(trajectory_list[[j]],list(centroid_3d))
          
          layer_list[[j]] <<- append(layer_list[[j]],i)
          
          assigned.flag <- TRUE
          
        }
      }
      
      if(!assigned.flag){
        
        trajectory_new <- list(centroid_3d)
        
        trajectory_list <<- append(trajectory_list,list(trajectory_new))
        
        layer_list <<- append(layer_list,list(i))
        
      }
      
    })
    
    cat("layer", z, " is done!\n")
    
    
  })
  
  return(trajectory_list)
  
}


# Iterate through trajectory_list to render trajectories.
Trajectory.plot <- function(trajectory_list, invisible_vector= NULL, type = "loess", text_flag = TRUE) {
  
  # INPUT : trajectory_list: initiated by "Trajectory.Init()“
  #         invisible_vector : trajectories in invisible_vector will not be rendered.
  #         type : "origin" or "loess", "origin" means link all the peaks directly. "loess" smooth each trajectory using LOESS.
  #         text_flag: plot the name of each trajectory. TRUE or FALSE.
  
  # OUTPUT: NULL
  
  threshold <- 1  # "The minimum number of peaks required in each trajectory; 
                  #trajectories with fewer peaks than this value will not be rendered."
  
  draw_matrix <- NULL
  
  if(type == "loess"){
    
    threshold <- 7 # "The minimum number of peaks required in each trajectory; 
                    #trajectories with fewer peaks than this value will not be rendered."
    
  }
  
  smooth <- 0.3 # the smooth of loess
  
  for(i in seq_along(trajectory_list)) {
    
    if( i %in% invisible_vector) next
    
    if(length(trajectory_list[[i]]) >= threshold){
      
      draw_matrix <- do.call(rbind,lapply(trajectory_list[[i]],st_coordinates))
      
      if(type == "origin") {
        
        lines3d(x=draw_matrix)
        
      }
      else if(type == "loess"){
        
        points_df <- as.data.frame(draw_matrix)
        colnames(points_df) <- c("X", "Y", "Z")
        loess_x <- loess(X ~ Z, data = points_df, span = smooth)
        loess_y <- loess(Y ~ Z, data = points_df, span = smooth)
        smooth_x <- predict(loess_x, newdata = data.frame(Z = points_df$Z))
        smooth_y <- predict(loess_y, newdata = data.frame(Z = points_df$Z))
        draw_matrix <- cbind(smooth_x,smooth_y,points_df$Z)
        lines3d(draw_matrix, col = "blue", lwd = 2.5)
        
        # 3D arrows
        n <- nrow(draw_matrix)
        
        point1 <- draw_matrix[n-1,]
        point2 <- draw_matrix[n,]
        height <- 100
        direction <- point2 - point1
        norm_direction <- sqrt(sum(direction^2))
        unit_direction <- direction / norm_direction
        cone_base <- point2 - height * unit_direction
        point2 <- point2 + 100 * unit_direction
        
        # arrow <- cylinder3d(center = rbind(cone_base,point2), radius=c(radius,0), sides=12, clsed = TRUE, color ="green")
        arrow3d(p0 = cone_base, p1 = point2, type = "rotation", s =1, theta = pi/10, n = 4, lit = FALSE, color = "blue")

      }
      else{
        
        stop("type is wrong!")
        
      }
      
      if(text_flag) {
        
        texts3d(x=draw_matrix[1,1],y=draw_matrix[1,2],z=draw_matrix[1,3],
                texts = as.character(i),pos=3,cex=1,col = "green")
      }
      
    }
    
  }
  
}


#=============================== RUN CODE HERE ===================================

trajectory_list <- Trajectory.Init(STKDE_result$ASTKDE, metric_list)

{
  
  open3d()
  Visualization.Cube(metric_list)
  Trajectory.plot(trajectory_list, invisible_vector = NULL, type = "loess", text_flag=TRUE)
  
}

{
  
  open3d()
  Visualization.Cube(metric_list)
  STKDE.Temporal.Marginal.density(ASTKDE_cond, metric_list, alpha_vector = c(0.95,0.99), plot_flag=TRUE)
  
}
