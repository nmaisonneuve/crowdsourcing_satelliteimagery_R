# EM Algorithm inspired from the article:
# Dawid and Skene "Maximum likelihood estimation of 
#  observer error-rates using the EM algorithm" (1979)
#  http://www.jstor.org/stable/10.2307/2346806

# input data= data frame with 
# - column: worker, with the index of the worker. 
# - column: item, with the index of the item. note  index>=1
# - column: label, with the index of the label. note  index>=1
#note all the index should start >=1

# example of input:
#  Worker,  item , label
#  1 1 1 
#  2 3 2
#  1 4 1
#  1 5 3

#  result=EM_classification(data,nb_iterations=10)
#

EM_classification=function(data,nb_iterations=10){

  nb_labels=length(unique(data$label))
  
  # in case the first category start at 0
  data$label[data$label==0]=nb_labels 
  
  categories=1:nb_labels

  # initialize the object probability matrix [nb_objects x nb_categories]
  # row idx = object idx, column idx= probability that the object belongs to the category idx
  items_pb  =init_items(data)

  # initialize the prior probability for each category 
  prior=init_prior(data)

  # initialize the error rates matrix 
  # row template:  workerID, category_from, category_to,  pb_error
  # by default  if from==to 0.9  else 0.1/(|Categories|-1)
  error_rates=init_error_rates(data)

  working=init_working(data)
  
  # for each iterations   
  for(i in 1:nb_iterations){
    
  working=fill_working(working,error_rates,categories)
  
  #estimate class proba for each item
  items_pb=estimate_item2(working,prior)
  
  #estimate prior
  prior=estimate_prior(items_pb)
  
  #estimate error rates for each worker
    for (worker in unique(error_rates$worker)){
      error_rates=estimate_individual_error_rates(worker,categories,working,prior,error_rates,data)
    }
  }
  
  #to remove if the category> binary class
  result=ifelse(items_pb[,1]>0.5,TRUE,FALSE)
  
  out=list(out = result, items_pb=items_pb, prior = prior, error_rates = error_rates, nb_iterations=nb_iterations)
  
return (out)
}


# initialize the object probability matrix [nb_objects x nb_categories]
# row idx = object idx, column idx= probability that the object belongs to the category idx
init_items=function(input){
  
  nb_items=length(unique(input$item))
  nb_labels=length(unique(input$label))  
 
  results= matrix(1/nb_labels,nrow=nb_items,ncol=nb_labels)
  return (results)
}

# initialize the prior probability for each category 
init_prior=function(input){
  size=length(unique(input$label))  
  return (rep(1/size,size))
}

# initialize the error rates matrix 
# row template:  workerID, category_from, category_to,  pb_error
# by default  if from==to 0.9  else 0.1/(|Categories|-1)
init_error_rates=function(input){
  
  nb_workers=length(unique(input$worker))
  nb_labels=length(unique(input$label))  
  
  #create it
  worker_id=rep(unique(input$worker), each = nb_labels*nb_labels)
  to=rep(1:nb_labels, times= nb_workers*nb_labels)
  from=rep(1:nb_labels, each=nb_labels, times= nb_workers)
  error_rates=data.frame(worker=worker_id, from=from, to=to, error=0)
  
  # fill it   
  error_rates$error[error_rates$to==error_rates$from]=0.9
  error_rates$error[error_rates$to!=error_rates$from]=0.1/(nb_labels-1)
  
  return(error_rates)
}

#init a temporary matrix from the data
#each row: itemIDX, workerIDX, labelIDX, true_labelIDX, error_rate
init_working=function(data){  
  nb_labels=length(unique(data$label))
  data2=data[rep(1:nrow(data),each=nb_labels),]
  data2$true_label=rep(1:nb_labels, times=nrow(data2)/nb_labels)
  data2$evidence=1
  return (data2)
}


#fill the temporary matrix with the new error_rates 
fill_working=function(working,error_rates,categories){
  for (worker in unique(working$worker)){
      for (label in categories){
        for (true_label in categories){
          error=error_rates$error[error_rates$worker==worker & error_rates$from==true_label & error_rates$to==label]          
          working$evidence[working$worker==worker & working$label==label & working$true_label==true_label]=error     
        }
      }
  }    
return (working)
}


# estimate the prior probability
estimate_prior=function(items_pb){
  return (colSums(items_pb)/nrow(items_pb))
}

# equation Dawid and Skene 2.5  (vectorized)
estimate_items=function(working,prior){
    
  # compute the product of the error rates of all the volunteers for all the items and the categorises
  detailed_items=aggregate(evidence~item+true_label,working, prod)  
  items=cast(detailed_items,item~true_label)  
  items$item=NULL #remove   
  
  # (we add the prior)
  for (i in 1:length(prior){
    items[,i]=items[,i]*prior[i]
  }  
  # we normalize
  items=items/rowSums(items)
    
  return(items)
}



# estimate the error rates for a given worker (worker)
estimate_individual_error_rates=function(worker,categories,working,prior,error_rates,data){
  
  # for the items labeled by the worker,we recompute locally 
  # the items pb matrix without the influence of the worker 
  # and that we use as gold standard
  items=data$item[data$worker==worker]
  working=working[working$item %in% items & working$worker!=worker,]
  items_worker=estimate_items(working,prior)
  
  # + the actual label given by the worker
  items_worker$to=data$label[data$worker==worker]
  
  for (from in categories){
    for (to in categories){
      # we update the error rate 
      error_rates$error[error_rates$worker==worker & error_rates$from==from & error_rates$to==to]=sum(items_worker[items_worker$to==to,from])
    }
    # we then normalize 
    idx=error_rates$worker==worker & error_rates$from==from
    error_rates$error[idx]=error_rates$error[idx]/sum(error_rates$error[idx])
  }  
  return (error_rates)
}





