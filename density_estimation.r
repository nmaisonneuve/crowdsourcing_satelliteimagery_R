#############################  Global DENSITY  ESTIMATION ####################
# for each worker inside the dataframe we compute the density of his contribution 
# using the mean of the k nearest-point for each worker 
# inputs:
# data : the spatialdataframe
# 
democratic_global_density=function(data){   
  workers=unique(data$workerID)
  density=c()
 for (worker in workers){
   worker_points=data[data$workerID==worker,]
   density=append(density,global_density(worker_points))
  }
  return (mean(density))
}


global_density=function(ref, tolerance=1){
   dist=seq(0.001,0.02,0.0005)
   dist_matrix= spDists(ref, ref, longlat = TRUE)   
   dist_matrix[which(dist_matrix<0.00001)]<-NA
   nndist=apply(dist_matrix, 2, function(x) {min(x, na.rm=TRUE)})     
   freq=apply(as.matrix(dist),1,FUN=function(x){ return(length(which(nndist>x)))} )
   freq=freq/nrow(ref)
   min_idx=which(freq<tolerance)[1]   
   #result=cbind(dist, freq)
  return (nndist[min_idx])
}

#############################################################################


#############################  LOCAL DENSITY  ESTIMATION ####################
# for a given point in a space we try to estimate the local density 
# using the mean of the k nearest-points for each worker 
# 
# inputs:
# data : the spatialdataframe
# point_idx = idx of the referent point 
# k , the number of neightbor points to take into account 
democratic_local_density=function(data,point_idx,k){     
   workers=unique(data$workerID)
   local_density=c()
   for (worker in workers){
     worker_points=data[data$workerID==worker,]
     ref_point=data[point_idx,]
     local_density=append(local_density,(local_density(worker_points,ref_point,k)))
  }
  return (mean(local_density))
}


local_density=function(data,point,k){
  if (nrow(data)>0){
    dist = spDistsN1(data, point, longlat = TRUE)
    dist=dist[dist>0.000000001] # we remove the point_idx
    return (sort(dist)[1:k])
  }else
    return(NA)
}
#############################################################################