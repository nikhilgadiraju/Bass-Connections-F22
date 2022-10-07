# Library Imports
library(data.table)
library(readxl)
library(rstatix)
library(jmv)
library(ggplot2)
library(patchwork)
library(ggeasy)
library(hash)
library(emmeans)
library(effectsize)
library(reticulate)
library(car)
library(stringr)

# Read voxel volumes w/ treatments excel file; remove the first 3 columns (index, filename, ID) and leave treatment and data columns; remove NaNs
data=read.csv('/Users/nikhilgadiraju/Box Sync/Home Folder nvg6/Sharing/Bass Connections/Processed Data/Mean Intensity & Voxel Volumes/voxelvolumes.csv')

  # Replace appended "X" to region names due to read.csv wrapper
old_colnames = colnames(data)[substr(colnames(data),1,1)=="X"][-1]
new_colnames = sub('.','',old_colnames)
colnames(data)[colnames(data) %in% old_colnames] <- new_colnames

data=data[  , -c(1,2,3)] 
dim(data)
data=na.omit(data) 

# Convert voxel counts to proportions of total brain volume
data[,2:dim(data)[2]]=100*data[,2:dim(data)[2]] 

# Create blank matrix to represent p values (column) for each brain regions (rows)
colnames_vec = c("FDR corrected Pvalue", "Effect Size Eta^2", 
                 "CI lower bound", "CI upper bound", 
                 "Mean sedentary group", "Mean voluntary group", "Mean voluntary + enforced group",
                 "SD sedentary group", "SD voluntary group", "SD voluntary + enforced group", 
                 "F-value", "Shapiro-Wilk Pvalue (norm)", "Levene Test Pvalue (homog)")
pvalsresults=matrix(NA,(dim(data)[2]-1), length(colnames_vec))
rownames(pvalsresults)=names(data)[2:dim(data)[2]]
colnames(pvalsresults)=colnames_vec

# Populate created matrix with ANOVA-generated p-values
# Returning Partial Eta Squared (PES) effect size since we're using ANOVAs for each given brain region. Note
# that Cohen's d can be used however those are typically utilized when doing a t-test and evaluating a difference 
# between two groups
# Reference: https://www.youtube.com/watch?v=e5od9tH_QUo

# Comparison between ANOVA and linear model: 
# https://stats.libretexts.org/Bookshelves/Applied_Statistics/Book%3A_Learning_Statistics_with_R_-_A_tutorial_for_Psychology_Students_and_other_Beginners_(Navarro)/16%3A_Factorial_ANOVA/16.06%3A_ANOVA_As_a_Linear_Model
len = dim(data)[2]
for (i in 1:(len-1))  {
  tempname=rownames(pvalsresults)[i]
  
  res.aov <- anova_test(get(tempname) ~ as.factor(Treatment), data = data)
  aov_table = get_anova_table(res.aov)
  
  lm <- lm(get(tempname) ~ as.factor(Treatment), data=data) 
  eff=eta_squared(lm, partial = FALSE)
  
  # Ho: data come from a normal distribution, H1: data do not come from a normal distribution
  # If p > 0.05, do NOT reject null, and thus data is normal
  normality = shapiro.test(lm$residuals)
  
  # Ho: variances are equal, H1: at least one variance is different
  # If p > 0.05, do NOT reject null, and thus data has equal variances (homogeneity == equality of variances)
  homogeneity = leveneTest(get(tempname) ~ as.factor(Treatment), data=data)

  means=by(data[,i+1],as.factor(data$Treatment), mean)
  sds=by(data[,i+1],as.factor(data$Treatment), sd)
  
  val_list = c(aov_table$p, eff$Eta2, eff$CI_low, eff$CI_high, means[1], means[2], means[3], sds[1], sds[2], sds[3], aov_table$F, normality$p.value, homogeneity$`Pr(>F)`[1]) #normality$p.value>0.05, homogeneity$`Pr(>F)`[1]>0.05 
  for (j in seq_along(val_list)){
    pvalsresults[i,j] <- val_list[j]
  }
}

# Create new 'pvalsresultadjusted' variable to eventually populated with FDR-corrected values
pvalsresultsadjusted=pvalsresults

###adjust pvalues Benjamini & Hochberg
# To understand what we're doing here, we are comparing all the p-values determined from the ANOVAs conducted on all
# 332 brain regions. Because we can treat each ANOVA as an individual hypothesis test, we need to account for the 
# multiple comparisons effect in the type 1 error rate. This can be done using various P-value (or significance) 
# correction methods, but we will be using the Benjamini & Hochberg Correction method (specified by either 'BH or 'fdr')
# References: https://www.youtube.com/watch?v=rZKa4tW2NKs&t=483s
#             https://www.youtube.com/watch?v=K8LQSvtjcEo&t=43s
pvalsresultsadjusted[,1] = p.adjust(pvalsresultsadjusted[,1], "fdr") #Benjamini & Hochberg


# Filter 'pvalsresultesadjusted' table to display brain regions that have significant p values (p<0.05)
#sig = pvalsresultsadjusted[pvalsresultsadjusted[,1]<=0.05,] 
sig = pvalsresultsadjusted[pvalsresultsadjusted[,1]<=0.05 & 
                           pvalsresultsadjusted[,ncol(pvalsresultsadjusted)]>0.05 &
                           pvalsresultsadjusted[,ncol(pvalsresultsadjusted)-1]>0.05,]

posthoc = matrix(NA,dim(sig)[1],7)
posthoc[,1]=sig[,1]

# Loop through each brain region and conduct a Tukey test and report p-values for each comparison group
for (i in 1:dim(sig)[1]) {
  tempname=rownames(sig)[i]
  res.aov <- aov(get(tempname) ~ Treatment, data = data)
  tuk=tukey_hsd(res.aov)
  treatment = c('sedentary', 'sedentary', 'wheel_only')
  control = c('wheel_only', 'treadmill', 'treadmill')
  for (j in 1:length(control)){
    hedges_out = hedges_g(get(tempname) ~ factor(Treatment, levels=c(control[j], treatment[j])), data=data)
    posthoc[i,j+4]=hedges_out$Hedges_g #Go through columns 5, 6, 7
  }
  posthoc[i,2:4]=tuk$p.adj
}

# Construct output CSV and save
colnames(posthoc)=c(colnames(sig)[1],"ST Comparison Group Pvalue","SW Comparison Group Pvalue","TW Comparison Group Pvalue",
                    "ST Comparison Group Effect Size","SW Comparison Group Effect Size","TW Comparison Group Effect Size")
rownames(posthoc)=rownames(sig)

# Write output
write.csv(sig, '/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/User-generated Files/posthoc_group_stats.csv')
write.csv(posthoc, '/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/User-generated Files/posthoc_comparison_stats.csv')

# Run read_posthoc.py
setwd("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22")
system("source BassVenv/bin/activate && python 'Statistical Analysis'/read_posthoc.py", wait=TRUE)

# Update Treatment names
data[,"Treatment"] = data[,"Treatment"] %>% str_replace_all(c("wheel_only" = "Voluntary", "treadmill" = "Voluntary + Enforced", "sedentary" = "Sedentary"))

# Define figure title hash/dict
dict <- hash()
dict[["st"]] = c("Sedentary", "Voluntary + Enforced")
dict[["sw"]] = c("Sedentary", "Voluntary")
dict[["tw"]] = c("Voluntary", "Voluntary + Enforced")

comp_groups = c('st','sw','tw')
plot_list = list()

# Read top regions CSVs
for (j in c('positive', 'negative')) {
  for (i in 1:length(comp_groups)) {
    comparison = comp_groups[i] # 'st', 'sw', or 'tw'
    top_comp=read.csv(paste('/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/User-generated Files/',comparison,'_regs/top_',substr(j,1,3),'_regions_',comparison,'.csv',sep=""))
    # data_comp = data[data$Treatment %in% dict[[i]],]
    
    # Basic Violin Plots
    sig_reg = top_comp$Abbreviation
    reg_struc = top_comp$Structure
    pvals_regs = formatC(top_comp$P.value, format = "e", digits = 2)
    eff_sizes = formatC(top_comp$Effect.Size, format = "e", digits = 2)
    
    p1 <- ggplot(data, aes_string(x="Treatment", y=sig_reg[1])) + 
      geom_violin() + geom_boxplot(width=0.1) +
      labs(title=reg_struc[1], subtitle=paste("P-value of ",toString(pvals_regs[1])," | Effect size of ",toString(eff_sizes[1])), y="", x="") + theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.title.y = element_text(margin = margin(r = 10)))
    
    p2 <- ggplot(data, aes_string(x="Treatment", y=sig_reg[2])) + 
      geom_violin() + geom_boxplot(width=0.1) + 
      labs(title=reg_struc[2], subtitle=paste("P-value of ",toString(pvals_regs[2])," | Effect size of ",toString(eff_sizes[2])), x="", y="") + theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.title.x = element_text(margin = margin(t = 10)))
    
    p3 <- ggplot(data, aes_string(x="Treatment", y=sig_reg[3])) + 
      geom_violin() + geom_boxplot(width=0.1) + 
      labs(title=reg_struc[3], subtitle=paste("P-value of ",toString(pvals_regs[2])," | Effect size of ",toString(eff_sizes[3])), y="", x="") + theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
    
    # Add total figure titles
    patchwork <- p1 + p2 + p3 + plot_annotation(
      subtitle = paste(dict[[comparison]][1],'vs.',dict[[comparison]][2],'Exercise'),
      theme = theme(plot.subtitle = element_text(hjust = 0.5))
    )
    plot_list[[i]] = patchwork
    
    # Saving Plots
    File <- paste("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Statistical Analysis/Output Figures/",j,"_eff/",comparison,'_',substr(j,1,3),'.png',sep="")
    # ggsave(File, plot = full_plot, width=1213, height=514, dpi = 150, units='px', scale=2)
  }
}

comp_plot = plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], nrow=3, ncol=1)
full_plot <- comp_plot + plot_annotation(
  title = paste('Regions of Significance following Post-Hoc Analysis (',str_to_title(j),' Effect Size)',sep=""),
  caption = 'DRAFT',
  theme = theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
)