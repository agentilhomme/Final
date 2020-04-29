setwd("C:/Users/Anais/OneDrive/Documents/UAF/Research3/Data/")

library(protr)
library(readxl)
library(seqinr)
library(stringr)

#gives table of the amino acid info text file 
aminos = read.table("amino-acid-info.txt", sep = "\t", header = T, stringsAsFactors = F) 
rownames(aminos) = aminos$abb1
aminos =  aminos[order(aminos$amino.acid),] # set the amino acid alphabetically

Patric_ID = read.table(paste0("All_Colwellia_patric_data.tsv"),fill = T,head = T,sep = "\t",colClasses = "character")
rownames(Patric_ID) = Patric_ID$All.Colwellia...all.my.strains...public.genome_id; Patric_ID$All.Colwellia...all.my.strains...public.genome_id = NULL; Patric_ID$genome.genome_id = NULL
Patric_ID[,2:9] = NULL
colnames(Patric_ID) = "strain"
Patric_ID$P.ID.New = rownames(Patric_ID)
rownames(Patric_ID) = Patric_ID$strain

OGT.Data = read.csv("Final.OGT.Data1.csv",header = T, fill = T)
OGT.Data = merge.data.frame(OGT.Data,Patric_ID,by = "strain", all = F)

#(0 = false, 1 = true,mw = molecular weigh
#hydropathicity = used to measure hydrophobic/hyrophilic

#Indices
arg_lys_ratio = function(x) {
  x["R",]/ (x["R",] + x["K",] )
}
aliphatic_index = function(x) {
  #normalize the protein to 1
  mole_fraction = sweep(x, MARGIN = 1, STATS = aminos$mw, FUN = `*`)
  total_mole_fraction = colSums(mole_fraction)
  mole_percent = (mole_fraction/total_mole_fraction)*100	
  mole_percent["A",] + 2.9*mole_percent["V",] + 3.9*(mole_percent["I",] + mole_percent["L",])
}
aromaticity = function(x) {
  colSums(sweep(x, MARGIN = 1, STATS = aminos$aromatic, FUN = `*`))	
}
acidic_residue = function(x) {
  x["D",] + x["E",]
}
gravy = function(x) {
  colSums(sweep(x, MARGIN = 1, STATS = aminos$hydropathicity, FUN = `*`))	
}
proline_residue = function(x) {
  x["P",]
}

#Create Histograms of Indices
orgfiles = list.files(".", pattern = "*.faa") 
output = data.frame()


for( file in orgfiles) {
  fasta_in = readFASTA(file)
  file_name = paste(file)
  file_name = str_remove(file_name,".faa")
  fasta_good = fasta_in[(sapply(fasta_in, protcheck))] # check if all the protein sequences are valid
  fasta_good = c2s(paste(fasta_good,sep = "",collapse = ""))
  aac = data.frame(extractAAC(fasta_good)) # gives amino acid composition in each amino acid sequence in genome, gives amino acid frequencies
  output[file_name,1] = arg_lys_ratio(aac)
  output[file_name,2] = aliphatic_index(aac)
  output[file_name,3] = aromaticity(aac)
  output[file_name,4] = acidic_residue(aac)
  output[file_name,5] = gravy(aac)
  output[file_name,6] = proline_residue(aac)
}
colnames(output) <- c("Arg-Lys.Ratio","Aliphatic.Index","Aromaticity","Acidic.Residue","GRAVY","Proline.Residue")
output$P.ID.New = rownames(output)
All.Data = merge(OGT.Data,output, by = "P.ID.New", all = F )
All.Data$PatricID = as.character(All.Data$PatricID)

par(mfrow = c(2,3))
for(i in 7:12) {
  plot.title = paste0(colnames(All.Data[i]))
  plot(All.Data$avg,All.Data[,i], main = plot.title,xlab = expression(paste("Temperature [",degree,"C]")),ylab = "Amino Acid Ratio")
}


