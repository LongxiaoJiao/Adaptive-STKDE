

#============================== MAIN FUNCTIONS =================================


# Plot the STKDE joint model at the specified statistical significance level, using contour lines.
Visualization.ContourLinesJoint <- function(STKDE_joint, alpha_vector, metric_list){
  
  dimensionxyz <- metric_list$dimensionxyz
  grx <- metric_list$grx
  gry <- metric_list$gry
  grz <- metric_list$grz
  
  value_vector <- quantile(STKDE_joint, alpha_vector, na.rm=TRUE)
  
  clines_red <- clines_yellow <- clines_blue <- list()
    
  for(i in seq(1,dimensionxyz[3],by=3)) {
    
    v_tranposed <- t(STKDE_joint[,,i])
    
    clines_red <- contourLines(x = grx, y = gry,z = v_tranposed, levels = value_vector[3])
    
    clines_yellow <- contourLines(x = grx, y = gry,z = v_tranposed, levels = value_vector[2])
    
    clines_blue <- contourLines(x = grx, y = gry,z = v_tranposed, levels = value_vector[1])
    
    for(line in clines_red){
      
      lines3d(x = line$x, y = line$y, z = grz[i], col = "red")

    }
    
    for(line in clines_yellow){
      
      lines3d(x = line$x, y = line$y, z = grz[i], col = "orange")
      
    }
    
    for(line in clines_blue){
      
      lines3d(x = line$x, y = line$y, z = grz[i], col = "blue")
      
    }
    
  }
  
}


# Plot the layer-by-layer STKDE spatial conditional model at the given statistical significance level, using contour lines.
Visualization.ContourLinesCond <- function(STKDE_cond, alpha_vector, metric_list){
  
  dimensionxyz <- metric_list$dimensionxyz
  grx <- metric_list$grx
  gry <- metric_list$gry
  grz <- metric_list$grz
  
  clines_outer <- clines_mid <- clines_inner <- NULL
  color_outer <- color_mid <- color_inner <- NULL
  
  for( i in seq(1, dimensionxyz[3], by=3)){
    
    matrix_this_layer <- t(STKDE_cond[,,i])
    
    threshold_this_layer <- quantile(matrix_this_layer, alpha_vector, na.rm=TRUE)
    
    color_vector <- c("blue","orange","red")
    line_width <- 0.5
    
    if(length(alpha_vector) == 3) {
      
      clines_outer <- contourLines(x=grx, y=gry, z=matrix_this_layer,levels=threshold_this_layer[1])
      clines_mid <- contourLines(x=grx, y=gry,z=matrix_this_layer,levels=threshold_this_layer[2])
      clines_inner <- contourLines(x=grx, y=gry, z=matrix_this_layer,levels=threshold_this_layer[3])
      
      color_outer <- color_vector[1]
      color_mid <- color_vector[2]
      color_inner <- color_vector[3]
      
    }
    else if(length(alpha_vector) == 2) {
      
      clines_outer <- contourLines(x=grx, y=gry, z=matrix_this_layer,levels=threshold_this_layer[1])
      clines_inner <- contourLines(x=grx, y=gry, z=matrix_this_layer,levels=threshold_this_layer[2])
      
      color_outer <- color_vector[2]
      color_inner <- color_vector[3]
    }
    else{
      
      clines_inner <- contourLines(x=grx, y=gry, z=matrix_this_layer,levels=threshold_this_layer[1])
      
      color_inner <- color_vector[3]
      
    }
    
    for (line in clines_outer) {
      
      lines3d(x=line$x,y=line$y,z=grz[i], col=color_outer,lwd=line_width)
      
    }
    for (line in clines_mid) {
      
      lines3d(x=line$x,y=line$y,z=grz[i], col=color_mid,lwd=line_width)
      
    }
    for (line in clines_inner) {
      
      lines3d(x=line$x,y=line$y,z=grz[i], col=color_inner,lwd=line_width)
        
    }

      
  }
  
}


# Plot the largest enclosing bounding rectangle of the spatiotemporal cube.
Visualization.Cube <- function(metric_list){
  
  boundary_x <- metric_list$boundary_x
  boundary_y <- metric_list$boundary_y
  time_range <- range(metric_list$grz)
  
  A <- c(boundary_x[1],boundary_y[1],time_range[1])
  B <- c(boundary_x[2],boundary_y[1],time_range[1])
  C <- c(boundary_x[2],boundary_y[2],time_range[1])
  D <- c(boundary_x[1],boundary_y[2],time_range[1])
  
  E <- c(boundary_x[1],boundary_y[1],time_range[2])
  F <- c(boundary_x[2],boundary_y[1],time_range[2])
  G <- c(boundary_x[2],boundary_y[2],time_range[2])
  H <- c(boundary_x[1],boundary_y[2],time_range[2])
  
  T <- c(boundary_x[1],boundary_y[1],time_range[2]+diff(time_range)/2)
  X <- c(boundary_x[2]+diff(boundary_x)/4,boundary_y[1],time_range[1])
  Y <- c(boundary_x[1],boundary_y[2]+diff(boundary_y)/4,time_range[1])
  
  segments3d(x=c(A[1],E[1],B[1],F[1],C[1],G[1],D[1],H[1]),y=c(A[2],E[2],B[2],F[2],C[2],G[2],D[2],H[2]),
             z=c(A[3],E[3],B[3],F[3],C[3],G[3],D[3],H[3]),col="black",lwd=0.1) #AE BF CG DH
  
  segments3d(x=c(A[1],B[1],D[1],C[1],H[1],G[1],E[1],F[1]),y=c(A[2],B[2],D[2],C[2],H[2],G[2],E[2],F[2]),
             z=c(A[3],B[3],D[3],C[3],H[3],G[3],E[3],F[3]),col="black", lwd=0.1) #AB DC HG EF
  
  segments3d(x=c(A[1],D[1],B[1],C[1],F[1],G[1],E[1],H[1]),y=c(A[2],D[2],B[2],C[2],F[2],G[2],E[2],H[2]),
             z=c(A[3],D[3],B[3],C[3],F[3],G[3],E[3],H[3]),col="black",lwd=0.1) #AD BC FG EH
  
  
}


# Plot 3D isosurfaces of the joint model based on the specified statistical significance level.
Visualization.Isosurface <- function(STKDE_joint, visible_range, alpha_vector, metric_list ){
  
  # "Render the 3D isosurface of the joint model at the specified significance levels in alpha_vector 
  #  within the given time period visible_range."
  
  grx <- metric_list$grx
  gry <- metric_list$gry
  grz <- metric_list$grz
  
  levels <- quantile(STKDE_joint, alpha_vector, na.rm = TRUE)
  
  if(is.null(visible_range)) visible_range <- 1:dim(STKDE_joint)[3]
  
  bottom_layer <- visible_range[1]
  bottom_matrix <- t(STKDE_joint[,,bottom_layer])
  bottom_grz <- grz[bottom_layer]
  
  STKDE_joint <- STKDE_joint[,,visible_range]
  grz <- grz[visible_range]
  
  STKDE_joint <- aperm(STKDE_joint,c(2,1,3))
  
  STKDE_joint[is.na(STKDE_joint)] <- 0
  
  if(length(alpha_vector) == 1){
    
    contour3d(f=STKDE_joint,level=levels,x=grx,y=gry,z=grz,color = "blue", alpha =0.2,
              add=TRUE,engine="rgl",smooth = 0 , lit = FALSE)
    
  }
  else if(length(alpha_vector) == 2){
    
    contour3d(f=STKDE_joint,level=levels,x=grx,y=gry,z=grz,color = c("orange","red"), alpha =c(0.35,0.5),
              add=TRUE,engine="rgl",smooth = 0 ,material = "dull", lit = FALSE)
    
  }
  else if(length(alpha_vector) == 3){
    
    contour3d(f=STKDE_joint,level=levels,x=grx,y=gry,z=grz,color = c("blue","orange","red"), alpha =c(0.2,0.3,0.4),
              add=TRUE,engine="rgl",smooth = 0, lit=FALSE) #material = "dull",alpha =c(0.2,0.35,0.5), 或 alpha =c(0.1,0.2,0.3)
    
  }
  else{
    
    stop(" levels should be a vector less than length 3")
    
  }
  
  bottom_lines <- contourLines(x = grx, y = gry, z = bottom_matrix, levels=levels[1])
  
  for(line in bottom_lines){
    # 
    # print(line$level)
    # 
    # lines3d(x=line$x,y=line$y,z=rep(bottom_grz,length(line$x)), col = "black")
    
    polygon3d(x=line$x,y=line$y,z=rep(bottom_grz,length(line$x)),
              fill =TRUE,col="blue",alpha = 0.1,lit=FALSE)
    
  }
  
}


# Plot the maximum projected envelope line for the hotspot region within the specified time range.
Visualization.ContourLines.Union <- function(STKDE_joint, layers, alpha, color, metric_list){
  
  grx <- metric_list$grx
  gry <- metric_list$gry
  grz <- metric_list$grz
  
  polygon_list <- new_polygon_list <- list()
  
  area <- 0
  
  if(length(layers) == 1) {
    
    return
    
  }
  
  else{
    
    level <- quantile(STKDE_joint, alpha, na.rm=TRUE)
    
    i <- 1
    
    layer_matrix <- t(STKDE_joint[,,layers[i]])
    
    contour_lines <- contourLines(x=grx,y=gry,z=layer_matrix,levels = level)
    
    polygon_list <- contour_lines
    
    for (i in seq_along(layers)[-1]){
      
      layer_matrix <- t(STKDE_joint[,,layers[i]])
      
      contour_lines <- contourLines(x=grx,y=gry,z=layer_matrix,levels = level)
      
      for(j in seq_along(contour_lines)){
        
        j_th_polygon <- st_polygon(x=list(cbind(contour_lines[[j]]$x,contour_lines[[j]]$y)))
        
        j_th_flag <- FALSE
        
        for(k in seq_along(polygon_list)) {
          
          k_th_polygon <- st_polygon(x=list(cbind(polygon_list[[k]]$x,polygon_list[[k]]$y)))
          
          if(as.logical(st_intersects(j_th_polygon,k_th_polygon,sparse = FALSE))){
            
            this_polygon <- st_union(j_th_polygon,k_th_polygon)
            
            x <- this_polygon[[1]][,1]
            
            y <- this_polygon[[1]][,2]
            
            polygon_list[[k]] <- list(x=x,y=y)
            
            j_th_flag <- TRUE
            
            break
            
          }
          
        }
        
        if(!j_th_flag) {

          polygon_list <- append(polygon_list,list(contour_lines[[j]]))
          
        }
        
      }
      
    }
  
  }
  
  #=============="Perform intersection checks within polygon_list."=============
  
  {
    
    polygons <- lapply(polygon_list, function(list) {
      
      st_polygon(list(cbind(list$x,list$y)))
      
    })
    
    sf_polygons <- st_sfc(polygons)
    
    intersects <- st_intersects(sf_polygons)
    
    # groups
    
    groups <- list()
    visited <- rep(FALSE, length(sf_polygons))
    
    for (i in seq_along(sf_polygons)) {
      if (!visited[i]) {
        group <- unlist(intersects[i], use.names = FALSE)
        groups <- append(groups, list(group))
        visited[group] <- TRUE
      }
    }
    
    merged_polygons <- lapply(groups, function(group) {
      st_union(sf_polygons[group])
    })
    
    
    new_polygon_list <- lapply(merged_polygons, function(polygon){
      
      coords <- st_coordinates(polygon)
      
      list(
        
        x=coords[,1],
        y=coords[,2]
      )
      
    })
    
  }
  
  #==============================="draw"===============================
  
  n <- length(layers)
  
  n <- n %/% 2
  
  for(lines in new_polygon_list){

    lines3d(x=lines$x,y=lines$y,z=ifelse(dim(STKDE_joint)[3] == 365,364 * 50 + 41640 *5 +25, 182 * 100 + 41640 * 5 + 50),
            col=color, lwd = 3)
    
    area<- area + 0.5 * abs(sum(lines$x[-1] * lines$y[-length(lines$y)] - lines$y[-1] * lines$x[-length(lines$x)]))
    
    #draw at the top
    
  }
  
  print(area)
  
  
}
  



#=============================== RUN CODE HERE ===================================

{
  open3d()
  Visualization.ContourLinesJoint(ASTKDE_joint, c(0.95,0.99,0.999), metric_list_half)
  Visualization.Cube(metric_list_half)
}

{
open3d()
Visualization.Cube(metric_list_half)
Visualization.Isosurface(ASTKDE_joint, visible_range = c(1:18), alpha_vector = c(0.95,0.99,0.999), metric_list_half)
Visualization.ContourLines.Union(ASTKDE_joint,
                                 layers = c(1:18), 
                                 alpha = 0.95, 
                                 color="blue", metric_list_half)
rgl.orthoview("xy",fov=0)
observer3d(0,-0,49000)
}
