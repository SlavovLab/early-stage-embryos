####-------------------- SOURCE FILE FOR EMBRYO PROJECT ----------------------------------------#####
 
###---------------------------------------------------------------------------------------------------------.
###--------------------------------------- LIBRARIES -------------------------------------------------------
###---------------------------------------------------------------------------------------------------------.
library(beeswarm)
library(circlize)
library(cluster)
library(ComplexHeatmap)
library(corrplot)
library(data.table)
library(dplyr)
library(factoextra)
library(fgsea)
library(ggbeeswarm)
library(ggdendro)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(impute)
library(lsa)
library(matrixStats)
library(pals)
library(plyr)
library(psych)
library(RColorBrewer)
library(readxl)
library(reshape2)
library(robustbase)
library(stringi)
library(stringr)
library(wesanderson)
library(writexl)


my_col2<-c("blue",rgb(0,0,1,0.5),"white",rgb(1,0,0,0.5),"red")
my_col3 = c("purple", "white", "yellow")


###---------------------------------------------------------------------------------------------------------.
###--------------------------------------- FUNCTIONS -------------------------------------------------------
###---------------------------------------------------------------------------------------------------------.

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

remove.duplicates<-function(data,Cols){
  
  return(data[!duplicated(data[,Cols]),])
  
}

cv<-function(x){
  
  sd(x, na.rm=T) / mean(x, na.rm=T)
  
}

cvna<-function(x){
  
  sum(!is.na(x))
  
}


# K - nearest neighbors imputation from Harrison Specht, SCoPE2
hknn<-function(dat, k){
  
  # Create a copy of the data, NA values to be filled in later
  dat.imp<-dat
  
  # Calculate similarity metrics for all column pairs (default is Euclidean distance)
  dist.mat<-as.matrix( dist(t(dat)) )
  #dist.mat<- 1-as.matrix(cor((dat), use="pairwise.complete.obs"))
  
  #dist.mat<-as.matrix(as.dist( dist.cosine(t(dat)) ))
  
  # Column names of the similarity matrix, same as data matrix
  cnames<-colnames(dist.mat)
  
  # For each column in the data... 
  for(X in cnames){
    
    # Find the distances of all other columns to that column 
    distances<-dist.mat[, X]
    
    # Reorder the distances, smallest to largest (this will reorder the column names as well)
    distances.ordered<-distances[order(distances, decreasing = F)]
    
    # Reorder the data matrix columns, smallest distance to largest from the column of interest
    # Obviously, first column will be the column of interest, column X
    dat.reordered<-dat[ , names(distances.ordered ) ]
    
    # Take the values in the column of interest
    vec<-dat[, X]
    
    # Which entries are missing and need to be imputed...
    na.index<-which( is.na(vec) )
    
    # For each of the missing entries (rows) in column X...
    for(i in na.index){
      
      # Find the most similar columns that have a non-NA value in this row
      closest.columns<-names( which( !is.na(dat.reordered[i, ])  ) )
      
      #print(length(closest.columns))
      
      # If there are more than k such columns, take the first k most similar
      if( length(closest.columns)>k ){
        
        # Replace NA in column X with the mean the same row in k of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns[1:k] ] )
        
      }
      
      
      # If there are less that or equal to k columns, take all the columns
      if( length(closest.columns)<=k ){
        
        # Replace NA in column X with the mean the same row in all of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns ] )
        
      }
      
      
    }
    
    # Populate a the matrix with the new, imputed values
    dat.imp[,X]<-vec
    
  }
  
  return(dat.imp)
  
}






rowed <- function(x, start.num, end.num, static.num) {
  for ( i in start.num:end.num) { 
    x[,i] <- x[,i] / (x[,static.num])
  }
  return(x)
}

rowed_log2 <- function(x, start.num, end.num, static.num) {
  for ( i in start.num:end.num) { 
    x[,i] <- x[,i] - (x[,static.num])
  }
  return(x)
}



columned <- function(x, start.num, end.num) {
  for ( i in start.num:end.num) { 
    x[,i] <- x[,i] / median(as.matrix(x[,i]), na.rm=TRUE)
  }
  return(x)
}

columned_log2 <- function(x, start.num, end.num) {
  for ( i in start.num:end.num) { 
    x[,i] <- x[,i] - median(as.matrix(x[,i]), na.rm=TRUE)
  }
  return(x)
}


# sums up values, ignoring NA's 
sum.na <- function(x) {
  y <- sum(x, na.rm=T)
}

var.na <- function(x) {
  y <- var(x, na.rm=T)
}

mean.na <- function(x) {
  y <- mean(x, na.rm=T)
}


CorrFil <- function(prepp, obsThreshold){
  
  # count the number of shared observations between protein (i.e how many single cells have both proteins quantified)
  obser <- pairwiseCount(t(prepp))
  # obtain the upper triangle of shared observation matrix 
  obser[lower.tri(obser, diag = T)] <- 188 
  obser.m <- reshape2::melt(obser)
  obser.m$Var1Var2 <- paste(obser.m$Var1, " ", obser.m$Var2)
  
  # obtain the correlation matrix 
  ppCor <- cor(t(prepp), method = "spearman", use = "pairwise.complete.obs")
  # obtain the upper trianlge of the correlation matrix
  ppCor1 <- ppCor
  ppCor1[lower.tri(ppCor1, diag = T)] <- 188
  ppCor1.m <- reshape2::melt(ppCor1)
  ppCor1.m$Var1Var2 <- paste(ppCor1.m$Var1, " ",  ppCor1.m$Var2)
  
  # merge the two matrices together 
  ppCor1.m$obser <- obser.m$value[match(ppCor1.m$Var1Var2, obser.m$Var1Var2)]
  # filter for the number of minimum shared observations, and for the upper triangle 
  ppCor1.m.fil <- ppCor1.m %>% dplyr::filter(obser > obsThreshold & value != 188)
  
  ggplot(data = ppCor1.m.fil, aes(x = value)) + geom_histogram(color = "black", fill = "grey") +
    geom_vline(xintercept = 0) + 
    theme_classic() + 
    labs(title = "Correlations with at least x observations", 
         subtitle = "filtered for upper triangle")
  
  return(ppCor1.m.fil)
}

