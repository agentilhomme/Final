# ------------------Prepping packages and directory --------------------------------------------------------
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

# ------------------------------------- Loading files and functions------------------------------------------
treefile="colwellia80+_100_treeWithGenomeIds.nwk"
genomes = read.table(paste0("colwellia_patric_data.tsv"),fill=T,head=T,sep="\t",colClasses="character")
init.save = read.table("rat83_params_test.tsv",fill=T,head=T,sep="\t",colClasses="character")
colnames(init.save) = c("strn","j","b", "Tmin", "c", "Tmax", "topt") 
# read in data
grd = read.table("growthrates.tsv",sep="\t")

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
#-------------------------Calculate the CI and SE error -------------------------------------------------------------
init.save = as.data.frame(init.save)
strain_rep_ls = split(init.save,as.factor(init.save$strn))
strain_rep_bounds = list()

for (i in 1:length(strain_rep_ls)){
  data = strain_rep_ls[[i]]
  str = paste0(data[1,1])
  UB = as.numeric(as.character(data[1:50,"topt"]))+1.96*sd(as.numeric(as.character(data[1:50,"topt"])))
  LB = as.numeric(as.character(data[1:50,"topt"]))-1.96*sd(as.numeric(as.character(data[1:50,"topt"])))
  strain_rep_bounds_matrix = cbind.data.frame(str,UB,LB)
  colnames(strain_rep_bounds_matrix) = c("strain","UB","LB")
  strain_rep_bounds_matrix$CI = strain_rep_bounds_matrix$UB - strain_rep_bounds_matrix$LB
  strain_rep_bounds_matrix$avg = mean(as.numeric(as.character(data[1:50,"topt"])))
  strain_rep_bounds[[i]] = strain_rep_bounds_matrix
}
strain_rep_bounds = do.call(rbind,strain_rep_bounds)
strain_CI_data = strain_rep_bounds[,c("strain","CI","avg")]
strain_CI_data = unique.data.frame(strain_CI_data)
strain_CI_data$se = strain_CI_data$CI/2
strain_CI_data$strain = as.character(strain_CI_data$strain)
row.names(strain_CI_data) = NULL

#-------------------------Clean up strain names so that the structures all match and clean up dataframe  -------------------------------------------------------------
strain_CI_data$strain = gsub("[.]","-",strain_CI_data$strain)
for(i in 1:nrow(strain_CI_data)){
  strain = as.character(strain_CI_data[i,1])
  if(grepl("Colwellia",strain,fixed = T) ==F){
    strain_CI_data[i,1] = paste0("Colwellia sp. ",strain)
  }
}
# use an example strain for each main type of Colwellia to plot 
strain_CI_data$strain = gsub("Mb3u","MB3u",strain_CI_data$strain)
strain_CI_data$strain = gsub("Colwellia psychrerythraea","Colwellia psychrerythraea 34H",strain_CI_data$strain)
strain_CI_data$strain = gsub("Colwellia hornerae","Colwellia hornerae strain ACAM 607",strain_CI_data$strain)
strain_CI_data$strain = gsub("Colwellia piezophila","Colwellia piezophila ATCC BAA-637",strain_CI_data$strain)
strain_CI_data$strain = gsub("Colwellia demingiae","Colwellia demingiae strain ACAM 459",strain_CI_data$strain)
uni.CI.data = unique(strain_CI_data)
uni.CI.data$CI = as.numeric(uni.CI.data$CI)
uni.CI.data$avg = as.numeric(uni.CI.data$avg)
uni.CI.data$se = as.numeric(uni.CI.data$se)
uni.CI.data = uni.CI.data %>% mutate_if(is.numeric,~round(.,2))
uni.CI.data = unique(uni.CI.data)

# remove MB02u-9,MB3u-43,BRX10-3,MB2u-11 as the ratkowski models give greater than 50 degree error bar 
uni.CI.data = uni.CI.data[!grepl("MB02u-9",uni.CI.data$strain),]
uni.CI.data = uni.CI.data[!grepl("MB3u-43",uni.CI.data$strain),]
uni.CI.data = uni.CI.data[!grepl("BRX10-3",uni.CI.data$strain),]
uni.CI.data = uni.CI.data[!grepl("MB02u-11",uni.CI.data$strain),]

genomes[,3:10] = NULL
colnames(genomes) = c("PatricID","strain")
All.Data = merge.data.frame(genomes,uni.CI.data, by = "strain",all = T)
All.Data = All.Data[!is.na(All.Data$PatricID),]
row.names(All.Data) =All.Data$PatricID
write.csv(All.Data,"Final.OGT.Data.csv")

tre.phylo = read.tree(treefile)
dend = midpoint.root(tre.phylo)
dend = ladderize(dend)
dend = drop.tip(dend,"1380381.3") #drop long tip (poor genome assembly)

common_names = intersect(dend$tip.label,rownames(All.Data))

dend = keep.tip(dend, common_names) # keeping the tips that are in the data.frame 
All.Data = All.Data[common_names,]

dend = rename_taxa(dend,All.Data,PatricID,strain) # rename the ends. 
dend$node.label = NULL
dend2 <- ReadDendrogram(textConnection(write.tree(dend)))

All.Data = data.frame(id = dend$tip.label,avg = All.Data$avg,se = All.Data$se)

# fill in duplicate Bg.28 data
Bg.28.Data = c(12.71, 3.82)
All.Data[3,2:3] = Bg.28.Data
All.Data[4,2:3] = Bg.28.Data

# add max and min bounds
All.Data$semax = All.Data$avg + All.Data$se
All.Data$semin = All.Data$avg - All.Data$se

#Create color ramp
temp_values = round(seq(-5,25,.1),1)
brPal = colorRampPalette(c('blue','white','red'))
All.Data$color = brPal(81)[as.numeric(cut(All.Data$avg,breaks = 81))]

#------------------------------------Plotting-----------------------------------

p = ggtree(dend) + theme(legend.position = "none") #+ geom_cladelabel(node = 133,label = "Colwellia psychrerythrea Cluster",barsize = 2,color = "red",size = 2) 
  #geom_cladelabel(node = 157,label = "Colwellia Other",barsize = 2,color = "blue") +
  #geom_cladelabel(node = 83,label = "Colwelli hornerae Cluster",barsize = 2,color = "green") +
  #geom_cladelabel(node = 109,label = "Colwelli chuckchiensis Cluster",barsize = 2,color = "orange")
  

p + geom_facet(panel = "Average Temp",data = All.Data,geom = geom_pointrange,mapping = aes(x = avg, xmax = semax,xmin = semin, color = avg)) + 
  theme_tree2() + 
  theme(legend.position = "none") +
  xlim_expand(c(0,5),"Average Temp")


