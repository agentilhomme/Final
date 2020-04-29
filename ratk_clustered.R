# use Ratkowksy 1983 model to predict optimal growth temperatures
library(minpack.lm)
library(ape)
library(phytools)
library(DataCombine)
library(phylogram)
library(dendextend)
library(evobiR)
library(ggtree)
library(dplyr)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(adephylo)
library(phylobase)
library(phylosignal)
library(treeio)
library(DECIPHER)
setwd("C:/Users/Anais/OneDrive/Documents/UAF/Research3/ag_thesis-master/")

# read in data
grd = read.table("growthrates.tsv",sep="\t")
ClusterID = read.table("Cleaned Clusters.csv",sep="\t")


biospec = read.csv("pone.0153343.s004.csv")
biospec$rate.per.hour = 60*biospec$rate.per.minute

genome.list = as.character(unique(biospec$binomial.name[grepl("Colwellia",biospec$binomial.name)]))

biogrd = biospec[biospec$binomial.name %in% genome.list,]
biogrd = data.frame("strain"=biogrd$binomial.name,
                    "replicate"="OD1",
                    "temp"=biogrd$T.C,
                    "mumax"=biogrd$rate.per.hour,
                    "r2"=0.99)

grd = rbind(grd,biogrd)

colnames(grd) = c("strain","rep","temp","rate","r2")
grd$temp = as.numeric(as.character(grd$temp)) + 273 #calculate Kelvin temperature
grd$sqrate = sqrt(grd$rate) #Ratkowsky uses square root of rate
grd = na.omit(grd)
grd$sqperday = sqrt(grd$rate*24) #convert to 1/day from 1/hour
grd$weights = -1*log10(1-grd$r2) #calculate weights based on R^2 values


# rT = (cc * (T - T1)(1 - exp( k * (T - T2))))^2
# Ratkowsky et al. 1983 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC217594
# rT is the growth rate
# T the temperature in Kelvin
# T1 and T2 the minimum and maximum temperatures at which rate of growth is zero
# sqrtcc) * k1 is the slope of the regression 
# k is a constant

rat83 = function(x, cc, T1, k, T2) {
  rT = cc * (x - T1) * (1 - exp(k * (x - T2)))
  return(rT)
}

rat82a = function(x, b, T0) {
  r = b*(x - T0)
}

#Clean up data frame for matching names 
grd$strain = gsub("[.]","-",grd$strain)
for(i in 1:nrow(grd)){
  strain = as.character(grd[i,1])
  if(grepl("Colwellia",strain,fixed = T) ==F){
    grd[i,1] = paste0("Colwellia sp. ",strain)
  }
}

# use an example strain for each main type of Colwellia to plot
grd$strain = gsub("Mb3u","MB3u",grd$strain)
grd$strain = gsub("Colwellia psychrerythraea","Colwellia psychrerythraea 34H",grd$strain)
grd$strain = gsub("Colwellia hornerae","Colwellia hornerae strain ACAM 607",grd$strain)
grd$strain = gsub("Colwellia piezophila","Colwellia piezophila ATCC BAA-637",grd$strain)
grd$strain = gsub("Colwellia demingiae","Colwellia demingiae strain ACAM 459",grd$strain)

#remove strains with extremely large confidence intervals 
grd = grd[!grepl("MB02u-9",grd$strain),]
grd = grd[!grepl("MB3u-43",grd$strain),]
grd = grd[!grepl("BRX10-3",grd$strain),]
grd = grd[!grepl("MB02u-11",grd$strain),]

#Merge data frames into one
GR_ClusterID = data.frame(strain = row.names(ClusterID),ClusterID); rownames(GR_ClusterID) = NULL
Growth_Data_Clustered = merge.data.frame(GR_ClusterID,grd,by = "strain",all = T)
colnames(Growth_Data_Clustered)[2] = "ClusterID"
GR_Clusters = split.data.frame(Growth_Data_Clustered,Growth_Data_Clustered$ClusterID)
Cluster_names = c("Colwellia hornerae Cluster", "Colwellia chuckiensis Cluster", "Colwellia psychrerythreae Cluster")

#set up color ramp 
temp_values = round(seq(-5,25,.1),1)
brPal = colorRampPalette(c('blue','white','red'))(length(temp_values))



#pdf(file="ratkowsky83_fits_test.pdf", width=5,height=5)

# average values from a previous run to use as initial values when needed
init.avg = list("b.est"=0.01111670, "Tmin.est"=182.30633469, "c.est"=0.05744382, "Tmax.est"=323.47254692)

for (k in 1:length(GR_Clusters)) {
  grd = GR_Clusters[[k]]
  
  # new fits will go here
  jmax = 1
  #init.save = matrix(NA,nrow=jmax*length(unique(grd$strain)), ncol=7)
  #strain_rep_ls = list()
  plot(1,type = "n", xlim=c(-10,30),ylim = c(0,4),xlab = "Temperature (°C)", ylab = "Growth Rate (1/day)")
  for(i in 1:length(unique(grd$strain))) {
    for(j in 1:jmax) {
      
      #initialize
      rt = NULL
      topt = NA
      
      
      strn = sort(unique(grd$strain))[i]
      print(as.character(strn))
      
      #subset data
      sgrd = grd[grd$strain == strn,c("temp","sqperday","weights")]
      sgrd = sgrd[sample(1:nrow(sgrd),nrow(sgrd),replace=TRUE),]
      sgrd = rbind(sgrd,aggregate(sgrd,by = list(sgrd$temp), FUN=mean)[,c("temp","sqperday","weights")])
      mintemp = sort(unique(sgrd$temp))[1]
      midtemp = sort(unique(sgrd$temp))[3]
      maxtemp = max(grd$temp)
      
      #estimate parameters before fitting, using directions in Ratkowsky 1983
      Tmin.est = (sgrd$sqperday[sgrd$temp == midtemp] * sgrd$temp[sgrd$temp == mintemp] - 
                    sgrd$sqperday[sgrd$temp == mintemp] * sgrd$temp[sgrd$temp == midtemp])/
        (sgrd$sqperday[sgrd$temp == midtemp] - sgrd$sqperday[sgrd$temp == mintemp])
      
      b.est = (sgrd$sqperday[sgrd$temp == mintemp] - sgrd$sqperday[sgrd$temp == midtemp])/
        (mintemp - midtemp)
      
      lhs1 = log(1 - sgrd$sqperday[sgrd$temp == mintemp]/(b.est * (maxtemp - Tmin.est)))
      lhs2 = log(1 - sgrd$sqperday[sgrd$temp == mintemp]/(b.est * (midtemp - Tmin.est)))
      
      Tmax.est = (lhs1*midtemp - lhs2*maxtemp)/(lhs1-lhs2)
      
      c.est = (lhs1 - lhs2)/(maxtemp - midtemp)
      
      Tmin.est = median(Tmin.est[Tmin.est > 0 & Tmin.est < mintemp & !is.nan(Tmin.est)])
      if(is.na(Tmin.est)) Tmin.est = 250
      b.est = median(b.est[b.est > 0 & !is.nan(b.est)])
      Tmax.est = median(Tmax.est[Tmax.est>0 & !is.nan(Tmax.est)])
      c.est = median(c.est[c.est>0 & !is.nan(c.est)])
      
      # plot rates
      #points(sgrd$temp - 273,sgrd$sqperday^2, pch=20, col=colors[i], bg=1,cex = sgrd$weights)
      
      # do fitting using nonlinear weighted least squares
      # first try using calculated starting values
      try({
        rt = nlsLM(sqperday ~ b * (temp - Tmin) * (1 - exp(cc * (temp - Tmax))),
                   data = sgrd,
                   weights = sgrd$weights,
                   start = list(b = b.est, Tmin = Tmin.est, cc = c.est, Tmax = Tmax.est),
                   algorithm = "port", trace = FALSE, control = list(maxiter=1000, tol=1e-6, minFactor=1e-5, printEval = TRUE, warnOnly = TRUE))
      }, silent=TRUE,)
      
      try({
        # if that fails, try using previous good estimates as starting values  
        if(is.null(rt)) {
          rt = nlsLM(sqperday ~ b * (temp - Tmin) * (1 - exp(cc * (temp - Tmax))),
                     data = sgrd,
                     weights = sgrd$weights,
                     start = list(b = init.avg$b.est, Tmin = init.avg$Tmin.est, cc = init.avg$c.est, Tmax = init.avg$Tmax.est),
                     algorithm = "port", trace = FALSE, control = list(maxiter=1000, tol=1e-6, minFactor=1e-5, printEval = TRUE, warnOnly = TRUE))
          
        }
      }, silent=TRUE)
      
      if(!is.null(rt)) {
        # finish plotting
        rtco = as.data.frame(t(coef(rt)))
        
        xs = seq(0,1000,0.1)
        rt.fit = rat83(xs, rtco$b, rtco$Tmin, rtco$cc, rtco$Tmax)
        rt.fit[rt.fit < 0] = NA
        
        xs = xs - 273
        #calculating uncertainties of opt to propagate error
        rt.fit.data <- as.data.frame(rt.fit)
        rt.fit.data <- na.omit(rt.fit.data)
        temp.fit.data <- xs[as.numeric(row.names(rt.fit.data))]
        
        topt = xs[which.max(rt.fit)]
        topt = round(topt,1)
        topt.index = which(temp_values == topt)
      
        #points(xs, rt.fit^2, type="l", col=brPal[topt.index],lwd=3)
        points(xs,rt.fit^2,type = "l",col = brPal[topt.index],lwd=3)
        
         
        # save fitted parameters
        #init.save[(i-1)*jmax+j,] = c(as.character(strn), j, b.est, Tmin.est, c.est, Tmax.est, topt)
        #strain_rep_matrix[i,] = c(as.character(strn), b.est, Tmin.est, c.est, Tmax.est, topt)
        
      } else {
        topt = sgrd$temp[which.max(sgrd$sqperday)]-273
        #init.save[(i-1)*jmax+j,] = c(as.character(strn), j, NA, NA, NA, NA, topt)
        #strain_rep_matrix[i,] = c(as.character(strn), NA, NA, NA, NA, topt)
      }
      #print(init.save[(i-1)*jmax+j,])
      # add title to plot (keep outside of try{})
      
    }
    #strain_rep_ls[j] = data.frame(strain_rep_matrix)
    
  }
  title(paste0(Cluster_names[k]))

}
