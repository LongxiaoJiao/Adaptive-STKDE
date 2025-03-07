



#============================== MAIN FUNCTIONS =================================

# Export the density cube to an external file for use by other software.
Export.file <- function(matrix, metric_list, filename) {
  
  grx <- metric_list$grx
  gry <- metric_list$gry
  grz <- metric_list$grz
  
  matrix[matrix==0] <- NA
  
  con <- file(filename, "w")
  
  cat("X,Y,Z,C\n",file = con)
  
  dims <- dim(matrix)
  
  for (x in 1:dims[2]) {
    
    start_time <- Sys.time()
    
    for (y in 1:dims[1]) {
      
      for (z in 1:dims[3]) {
        
        cat(grx[x],",",gry[y],",",grz[z],",", matrix[y, x, z], "\n", file=con)
        
      }
    }
    
    end_time <- Sys.time()
    
    cat("x=",x,"is done, time cost: ",difftime(end_time,start_time,units="secs"),"\n")
    
  }
  
  close(con)
  
}








#=============================== RUN CODE HERE ===================================

Export.file(STKDE_joint, metric_list_half, "STKDE_joint.csv")