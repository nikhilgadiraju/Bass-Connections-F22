# Library Imports
library(data.table)
library(readxl)
library(rstatix)
library(jmv)
library(ggplot2)
library(patchwork)
library(ggeasy)
library(hash)

# Read voxel volumes w/ treatments excel file; remove the first 3 columns (index, filename, ID) and leave treatment and data columns; remove NaNs
data=read.csv('/Users/nikhilgadiraju/Box Sync/Home Folder nvg6/Sharing/Bass Connections/Processed Data/Mean Intensity & Voxel Volumes/voxelvolumes.csv')
data=data[  , -c(1,2,3)] 
dim(data)
data=na.omit(data) 

# Convert voxel counts to proportions of total brain volume
data[,2:dim(data)[2]]=100*data[,2:dim(data)[2]] 

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

# Update Treatment names
gsub("wheel_only", "Voluntary", data[,"Treatment"])
gsub("treadmill", "Voluntary + Enforced", data[,"Treatment"])
gsub("sedentary", "Sedentary", data[,"Treatment"])

# Define figure title hash/dict
dict <- hash()
dict[["st"]] = "Sendentary vs. Voluntary + Forced Exercise"
dict[["sw"]] = "Sendentary vs. Voluntary Exercise"
dict[["tw"]] = "Voluntary vs. Voluntary + Forced Exercise"

# Read top regions CSVs
comparison = 'sw' # 'st', 'sw', or 'tw'
print(comparison)
top_comp=read.csv(paste('/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/User-generated Files/top_regions_',comparison,'.csv',sep=""))

# Basic Violin Plots
sig_reg = top_comp$Abbreviation
reg_struc = top_comp$Structure
pvals_regs = formatC(top_comp$p.values, format = "e", digits = 2)

p1 <- ggplot(data, aes_string(x="Treatment", y=sig_reg[1])) + 
  geom_violin() + geom_boxplot(width=0.1) +
  labs(title=reg_struc[1], subtitle=paste("p-value of ",toString(pvals_regs[1])), y="Normalized Regional Proportion (%)", x="") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.title.y = element_text(margin = margin(r = 10)))

p2 <- ggplot(data, aes_string(x="Treatment", y=sig_reg[2])) + 
  geom_violin() + geom_boxplot(width=0.1) + 
  labs(title=reg_struc[2], subtitle=paste("p-value of ",toString(pvals_regs[2])), y="", x="Treatment Conditions") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.title.x = element_text(margin = margin(t = 10)))

p3 <- ggplot(data, aes_string(x="Treatment", y=sig_reg[3])) + 
  geom_violin() + geom_boxplot(width=0.1) + 
  labs(title=reg_struc[3], subtitle=paste("p-value of ",toString(pvals_regs[3])), y="", x="") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

# Add total figure titles
patchwork <- p1 + p2 + p3
full_plot <- patchwork + plot_annotation(
  title = 'Regions of Significance following Post-Hoc Analysis',
  subtitle = dict[[comparison]],
  caption = 'DRAFT',
  theme = theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
)

# Saving Plots
File <- paste("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Statistical Analysis/Output Figures/",comparison,'.png',sep="")
ggsave(File, plot = full_plot, width=1213, height=514, dpi = 150, units='px', scale=2)
