source ('density_estimation.r')

democratic_clustering=function (data,max_dist=0.007,min_volunteers=2){
  clustered_data=democratic_clustering_internal(data,max_dist,min_volunteers)  
  clustered_data=clustered_data[clustered_data$cluster!=0,]
  centroids=compute_centroid(clustered_data)  
  return (centroids)
}

democratic_clustering_internal=function (data,max_dist=0.007,min_volunteers=2){  
  data$cluster=0 # by default all are noises  
  cluster_idx=1  
  workers=unique(data$workerID)
  size=nrow(data)  
  cv=integer(size)    
  unclass=(1:size)    
  
  workers_idx=lapply(as.list(workers),FUN=function(x){ which(data$workerID==x)})  
  for (i in 1:size) {
      if (cv[i]==0){
      unclass=(1:size)[cv<1]  
      others=setdiff(workers,data$workerID[i])
      points_idx=c(i)            
      for (other in others){  
        idx=match(other,workers)
        selected=intersect(workers_idx[[idx]],unclass)
        dist = spDistsN1(data[selected,], data[i,], longlat = TRUE)
        idx=which.min(dist)
        
        if ((length(dist)>0) && (dist[idx]<max_dist)){
          points_idx=append(points_idx,selected[idx])
        }       
      }
      if (length(points_idx)>=min_volunteers){
        data$cluster[points_idx]=cluster_idx
        cluster_idx=cluster_idx+1
        cv[points_idx]=1
      }else{        
        cv[i]=1
      }             
    }
  }    
  return (data)
}

# a bit faster 
democratic_clustering2=function (data,min_dist=0.07,min_volunteers=2){  
  #data=as.matrix(as.data.frame(data))
  data$free=TRUE
  data$cluster=0 # by default all are noises  
  workers_idx=unique(data$workerID)
  
  size=length(workers_idx)-min_volunteers+1    
  clusters=list()
  
  for (v in 1:size) {

    my_idx=which(data$workerID==workers_idx[v] & data$free==TRUE)          
    my_markers=data[my_idx,]          
  
    matches=c()    
    clusters_v=c()
     if (nrow(my_markers)>0){
        others_idx=workers_idx[-c(1:v)]     #only the next      
        
      for (other in others_idx){ #for each other volunteers                        
            other_idx=which(data$workerID==other & data$free==TRUE)
            others_selected=data[other_idx,]                        
            match=zerodist2(others_selected, my_markers,min_dist/112)
            match2=subset(match, !duplicated(match[,2])) 
            if (nrow(match)!=nrow(match2)){
            }
            match=match2            
            if (nrow(match)>0){
              matches=rbind(cbind(other_idx[match[,1]], match[,2], other),matches)
            }
     }
    if (length(matches)>0 && nrow(matches)>0){
    matches[,2]=my_idx[matches[,2]]
    matches=rbind(matches,cbind(my_idx,my_idx,workers_idx[v]))   # we add the volunteer's points
    }else{
      matches=cbind(my_idx,my_idx,workers_idx[v])
    }
    clusters_v=aggregate(matches[,1], list(matches[,2]),function(x){
    if (length(x)>=(min_volunteers)){ 
        data$free[x]<<-FALSE
      return (c(mean(data@coords[x,1]),mean(data@coords[x,2]), length(x)))  
    }      
    else {
      return(NA)
    }}, simplify=FALSE)$x        

    if (length(clusters_v)>0){
       clusters_v=clusters_v[!is.na(clusters_v)]
    }
    if (length(clusters_v)>0){        
        clusters=append(clusters,clusters_v)       
    }
   }
  }

  if (length(clusters)>0){
    clusters=matrix(unlist(clusters), ncol=3, byrow=TRUE)        
    if (nrow(clusters)==1){
      output=SpatialPointsDataFrame(coords=t(clusters[,1:2]), data=as.data.frame(t(clusters[,3])), proj4string=CRS(ps))
    }else{
      output=SpatialPointsDataFrame(coords=clusters[,1:2], data=as.data.frame(clusters[,3]), proj4string=CRS(ps))  
    }
  }else{
    output=SpatialPointsDataFrame(coords=cbind(0,0), data=data.frame(0), proj4string=CRS(ps))  
  }  
  names(output)=c('support')  
  return (output)
}
