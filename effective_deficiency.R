# Deficiency of networks with and without balanced complexes
setwd('~')
setwd('GitHub/balanced_complexes/effective_deficiency')

# load packages
library(igraph)
library(R.matlab)
library(Matrix)

source('deficiency')

files_irrev = Sys.glob('Results_irreversibility_considered/*any_kinetic.mat')
files_obj = Sys.glob('Results_objective/*any_kinetic.mat')

d_original = rep(-1,length(files_irrev))
d_eff_irrev = rep(-1,length(files_irrev))
d_eff_obj = rep(-1,length(files_irrev))

for(i in 1:length(files_irrev)){
  print(i)
  # load results  
  Results_irrev = readMat(files_irrev[i])
  Results_obj = readMat(files_obj[i])
  model_original = Results_irrev$MODEL.r[[1]]
  
  ###############################################################
  # calculate deficiency
  ###############################################################
  rev = model_original[[1]][,,1]$lb<0
  S = model_original[[1]][,,1]$S
  D_original = deficiency(S,rev)
  
  A=D_original$A
  Y=D_original$Y
  
  # number of complexes
  complexes=dim(A)[1]
  
  linkage_classes=count_components(D_original$g,mode="weak")
  s=rankMatrix(S)
  
  d_original[i] = complexes-linkage_classes-s[1]
  
  ##########################################################
  # calculate effective deficiency - irrev
  ##########################################################
  B = c(unlist(Results_irrev$B.out[[1]]), unlist(Results_irrev$B.r[[1]]))
  EB = matrix(0,length(B),dim(A)[1]) 
  EB[,B] = diag(1,length(B))
  
  EA = EB%*%A
  
  d_eff_irrev[i] = d_original[i] + s[1] - rankMatrix(rbind(S,EA)) 
  
  ##########################################################
  # calculate effective deficiency - obj
  ##########################################################
  B = c(unlist(Results_obj$B.out[[1]]), unlist(Results_obj$B.o[[1]]))
  EB = matrix(0,length(B),dim(A)[1]) 
  EB[,B] = diag(1,length(B))
  
  EA = EB%*%A
  
  d_eff_obj[i] = d_original[i] + s[1] - rankMatrix(rbind(S,EA)) 
  
}
