# Library Imports
library(data.table)
library(readxl)
library(rstatix)
library(jmv)

# Read voxel volumes w/ treatments excel file; remove the first 3 columns (index, filename, ID) and leave treatment and data columns; remove NaNs
data=read_excel('/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Statistical Analysis/voxelVolumes_treatments.xlsx')
data=data[  , -c(1,2,3)] 
dim(data)
data=na.omit(data) 

# Convert voxel counts to proportions of total brain volume
data[,2:dim(data)[2]]=100*data[,2:dim(data)[2]]/rowSums(data[,2:dim(data)[2]]) 

# Create blank matrix to represent p values (column) for each brain regions (rows)
pvalsresults=matrix(NA,(dim(data)[2]-1)  , 3 )
rownames(pvalsresults)=names(data)[2:dim(data)[2]]
colnames(pvalsresults)= c("treatment", "homogen", "Sha-Wilk norm")

# Populate cretaed matrix with ANOVA-generated p-values
len = dim(data)[2]
for (i in 1:(len-1))  {
   
  tempname=rownames(pvalsresults)[i]
  res.aov <- anova_test(get(tempname) ~ as.factor(Treatment), data = data)
  a = get_anova_table(res.aov)
  p = a$p
  pvalsresults[i,1] <- p
}

# Create new 'pvalsresultadjusted' variable to eventually populated with FDR-corrected values
pvalsresultsadjusted=pvalsresults

###adjust pvalues Benjamini & Hochberg
# for (j in 1:dim(pvalsresultsadjusted)[2]) {
pvalsresultsadjusted[,1] = p.adjust(pvalsresultsadjusted[,1], "fdr") #Benjamini & Hochberg
# }

# Filter 'pvalsresultesadjusted' table to display brain regions that have signficant p values (p<0.05)
sig = pvalsresultsadjusted[pvalsresultsadjusted[,1]<=0.05,1] 
posthoc=matrix(NA,length(sig),4)
posthoc[,1]=sig

# Loop through each brain region and conduct a Tukey test and report p-values for each comparison group
for (i in 1:length(sig)) {
  tempname=names(sig)[i]
  res.aov <- aov(get(tempname) ~ Treatment, data = data)
  tuk=tukey_hsd(res.aov)
  posthoc[i,2:4]=tuk$p.adj
}

# Construct output CSV and save
colnames(posthoc)=c("FDR","ST","SW","TW")
rownames(posthoc)=names(sig)
write.csv(posthoc, 'posthoc_standard.csv')
