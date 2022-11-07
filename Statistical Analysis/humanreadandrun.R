### LIBRARY IMPORTS
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
library(cowplot)
library(ggpubr)

### WHOLE BRAIN ANALYSIS
# Specify and read folder with all of Nariman's data (individual_label_statistics/)
path_vol="/Users/nikhilgadiraju/Box Sync/Home Folder nvg6/Sharing/Bass Connections/Data/individual_label_statistics/"
file_list=list.files(path_vol)

# Specify and read ID_Treatment.CSV - this CSV associates contains metadata for all the mice brain
# data files: Original ID, Modified ID, N-number, Treatment
path_metadata = "/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/User-generated Files/ID_Treatment.csv"
metadata = read.csv(path_metadata)

# *******
# Remove 'A19120503_T1' ID dataset since we don't have that file
metadata = metadata[!(metadata$Modified.ID=='A19120503_T1'),]
# *******

# Create empty matrix to hold total brain volume for each brain dataset
temp=read.delim( paste0(path_vol,file_list[1]) )
len=length(temp$volume_mm3)
vol_tab = matrix(NA,length(file_list),3)

# Populate matrix with brain dataset filename and total volume (mm3)
for (i in 1:length(file_list)) {
  temp=read.delim(paste0(path_vol,file_list[i]))
  vol_tab[i,3]=sum(temp$volume_mm3[2:len])
  vol_tab[i,2]=substr(file_list[i], 1, 6)
}

# Set column names and filter vol_tab to only display CVN mice we are analyzing (by N-numbers)
colnames(vol_tab) <- c('ID', 'N-number', 'Volume')
vol_tab <- vol_tab[which(metadata$N.number %in% vol_tab[,'N-number']),]

# Replace N-numbers with each dataset's treatment condition for downstream ANOVA
for (g in 1:length(vol_tab[,'N-number'])) {
  treat_row = which(metadata$N.number == metadata$N.number[g])
  vol_tab[,'N-number'][g] = metadata$Treatment[treat_row]
  vol_tab[,'ID'][g] = metadata$Modified.ID[treat_row]
}

# Rename 'N-number' column title with 'Treatment' column title; also convert volume column
# to a numeric column type
colnames(vol_tab)[2] = 'Treatment'
wb_data = data.frame(vol_tab)
wb_data[,3] = as.numeric(wb_data[,3])

# Conduct and report whole-brain ANOVA
wb_aov = anova(lm(Volume ~ Treatment, data=wb_data[,c('Treatment','Volume')]))

# Print TRUE if anova p-value < 0.05
paste('Whole-brain ANOVA result:',(wb_aov$`Pr(>F)` < 0.05)[1])


### BEGIN STATISTICAL ANALYSIS
# Read voxel volumes w/ treatments excel file
data=read.csv('/Users/nikhilgadiraju/Box Sync/Home Folder nvg6/Sharing/Bass Connections/Processed Data/Mean Intensity & Voxel Volumes/voxelvolumes.csv')

# Replace appended "X" to region names due to read.csv wrapper
old_colnames = colnames(data)[substr(colnames(data),1,1)=="X"][-1]
new_colnames = sub('.','',old_colnames)
colnames(data)[colnames(data) %in% old_colnames] <- new_colnames

# remove the first 3 columns (index, filename, ID) and leave treatment and data columns; omit NaNs
data=data[  , -c(1,2,3)] 
dim(data)
data=na.omit(data) 

# Convert voxel counts to proportions of total brain volume
data[,2:dim(data)[2]]=100*data[,2:dim(data)[2]] 

## Create blank matrix to represent p values (column) for each brain region (rows)
colnames_vec = c("Mean sedentary group", "Mean voluntary group", "Mean voluntary + enforced group",
                 "SD sedentary group", "SD voluntary group", "SD voluntary + enforced group", 
                 "Uncorrected Pvalue", "FDR corrected Pvalue", "F-value", "Cohen's F Effect Size", "Effect Size Eta^2", 
                 "CI lower bound", "CI upper bound", 
                 "Shapiro-Wilk Pvalue (norm)", "Levene Test Pvalue (homog)")
pvalsresults=matrix(NA,(dim(data)[2]-1), length(colnames_vec))
rownames(pvalsresults)=names(data)[2:dim(data)[2]]
colnames(pvalsresults)=colnames_vec

# Append whole-brain region to end of data
wb_data = wb_data[order(wb_data$Treatment, decreasing = TRUE), ]


# Populate created matrix with ANOVA-generated p-values
# Returning Partial Eta Squared (PES) effect size since we're using ANOVAs for each given brain region. Note
# that Cohen's d can be used however those are typically utilized when doing a t-test and evaluating a difference 
# between two groups
# Reference: https://www.youtube.com/watch?v=e5od9tH_QUo

# Comparison between ANOVA and linear model: 
# https://stats.libretexts.org/Bookshelves/Applied_Statistics/Book%3A_Learning_Statistics_with_R_-_A_tutorial_for_Psychology_Students_and_other_Beginners_(Navarro)/16%3A_Factorial_ANOVA/16.06%3A_ANOVA_As_a_Linear_Model
len = dim(data)[2]
for (i in 1:(len-1))  {
  # Set 'tempname' to the currently analyzed brain regions
  tempname=rownames(pvalsresults)[i]
  
  # Create a linear model object using a given brain region and the data associated with the three
  # categories: 'sedentary', 'treadmill', 'wheel_only'. Calculate eta-squared effect size for the overall model 
  mylm <- lm(get(tempname) ~ as.factor(Treatment), data=data) 
  eff=effectsize::eta_squared(mylm, partial = F)
  cf=effectsize::cohens_f(mylm)
  
  # Conduct ANOVA on linear model to evaluate whether there is a significant difference between exercise treatment groups
  # within a given brain region; expect the means to be different and reject null (u_sedentary = u_treadmill = u_wheel-only)
  aov_table = anova(mylm)
  
  # Ho: data come from a normal distribution, H1: data do not come from a normal distribution
  # If p > 0.05, do NOT reject null, and thus data is normal
  normality = shapiro.test(mylm$residuals)
  
  # Ho: variances are equal, H1: at least one variance is different
  # If p > 0.05, do NOT reject null, and thus data has equal variances (homogeneity == equality of variances)
  homogeneity = leveneTest(get(tempname) ~ as.factor(Treatment), data=data)

  # Calculate the mean and standard deviation for each treatment group within the current brain region dataset
  means=by(data[,i+1],as.factor(data$Treatment), mean)
  sds=by(data[,i+1],as.factor(data$Treatment), sd)
  
  # Output calculated values to the 'pvalresults' matrix in the user-defined order; the 'pvalresults' matrix will later
  # be used to create the 'posthoc_group_stats.csv' output
  val_list = c(means[1], means[2], means[3], sds[1], sds[2], sds[3], aov_table$`Pr(>F)`[1], aov_table$`Pr(>F)`[1], aov_table$'F value'[1], cf$Cohens_f, eff$Eta2, eff$CI_low, eff$CI_high, normality$p.value, homogeneity$`Pr(>F)`[1]) #normality$p.value>0.05, homogeneity$`Pr(>F)`[1]>0.05
  for (j in seq_along(val_list)){
    pvalsresults[i,j] <- val_list[j]
  }
}

# Create new 'pvalsresultadjusted' variable to eventually populated with FDR-corrected values
pvalsresultsadjusted=pvalsresults

## adjust pvalues Benjamini & Hochberg
# To understand what we're doing here, we are comparing all the p-values determined from the ANOVAs conducted on all
# 332 brain regions. Because we can treat each ANOVA as an individual hypothesis test, we need to account for the 
# multiple comparisons effect in the type 1 error rate. This can be done using various P-value (or significance) 
# correction methods, but we will be using the Benjamini & Hochberg Correction method (specified by either 'BH or 'fdr')
# References: https://www.youtube.com/watch?v=rZKa4tW2NKs&t=483s
#             https://www.youtube.com/watch?v=K8LQSvtjcEo&t=43s
pvalsresultsadjusted[,8] = p.adjust(pvalsresultsadjusted[,8], "fdr") #Benjamini & Hochberg


# Filter 'pvalsresultesadjusted' table to display brain regions that have significant p values (p<0.05)
#sig = pvalsresultsadjusted[pvalsresultsadjusted[,1]<=0.05,] 
sig = pvalsresultsadjusted[pvalsresultsadjusted[,8]<=0.05 & 
                           pvalsresultsadjusted[,ncol(pvalsresultsadjusted)]>0.05 & # For Levene's Test
                           pvalsresultsadjusted[,ncol(pvalsresultsadjusted)-1]>0.05,] # For Shapiro-Wilk Test
sig = as.data.frame(sig)
sig <- sig[order(sig$`Cohen's F Effect Size`, decreasing = TRUE),]

# Create empty 'posthoc' matrix to eventually populate with all gruop comparison data; 
# Set first 'posthoc' column equal to 8th column of 'sig' matrix (FDR-corrected P-value column)
posthoc = matrix(NA,dim(sig)[1],13)
posthoc[,1]=sig[,8]

# Loop through each significant (FDR-corrected p-value < 0.05) brain region (from 'sig' matrix) and conduct 
# a Tukey test and report p-values for each comparison group
for (i in 1:dim(sig)[1]) {
  # Set 'tempname' to the currently analyzed brain regions
  tempname=rownames(sig)[i]
  
  # Conduct ANOVA on currently analyzed brain region and use ANOVA output to conduct post-hoc Tukey test
  res.aov <- aov(get(tempname) ~ Treatment, data = data)
  tuk=tukey_hsd(res.aov)
  
  # Define Tukey comparison groups for calculating hedge's G effect size. Note that these groups are mannually
  # set so that effect sizes can correlate with changes in regional volumes (e.g. negative effect size correlates
  # with decrease in regional volumes)
  control = c('wheel_only', 'treadmill', 'treadmill')
  treatment = c('sedentary', 'sedentary', 'wheel_only')
  
  # For given brain region, loop through comparison groups and calculate hedge's g for each comparison groups
  for (j in 1:length(control)){
    hedges_out = hedges_g(get(tempname) ~ factor(Treatment, levels=c(control[j], treatment[j])), data=data)
    posthoc[i,j+4]=hedges_out$Hedges_g #Go through columns 5, 6, 7
  }
  
  # Add adjusted P-value and add low and high confidence intervals to 'posthoc' matrix
  posthoc[i,2:4]=tuk$p.adj
  posthoc[i,c(8,10,12)] = tuk$conf.low
  posthoc[i,c(9,11,13)] = tuk$conf.high
}

# Define output CSV column and row names
colnames(posthoc)=c(colnames(sig)[8],"ST Comparison Group Pvalue","SW Comparison Group Pvalue","TW Comparison Group Pvalue",
                    "ST Comparison Group Effect Size","SW Comparison Group Effect Size","TW Comparison Group Effect Size",
                    "ST Comparison Group Lower CI", "ST Comparison Group Higher CI", "SW Comparison Group Lower CI", "SW Comparison Group Higher CI", 
                    "TW Comparison Group Lower CI", "TW Comparison Group Higher CI")
rownames(posthoc)=rownames(sig)

# Save output CSV to specific folder string
write.csv(sig, '/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/User-generated Files/posthoc_group_stats.csv')
write.csv(posthoc, '/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/User-generated Files/posthoc_comparison_stats.csv')

# Activate Venv and run read_posthoc.py
# read_post.py serves to read above output CSV (posthoc_comparison_stats.csv specifically), and sort 
# brain regions by effect sizes and postivitiy/negativity
setwd("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22")
system("source BassVenv/bin/activate && python 'Statistical Analysis'/read_posthoc.py", wait=TRUE)

# Update Treatment names for plotting
data[,"Treatment"] = data[,"Treatment"] %>% str_replace_all(c("wheel_only" = "Voluntary", "treadmill" = "Voluntary + Enforced", "sedentary" = "Sedentary"))

# Define figure title hash/dict
dict <- hash()
dict[["st"]] = c("Sedentary", "Voluntary + Enforced")
dict[["sw"]] = c("Sedentary", "Voluntary")
dict[["tw"]] = c("Voluntary", "Voluntary + Enforced")

### PLOTTING
# Set comparisong group vector to loop through: 'st': sedentary vs. treadmill, 'sw': sedentary vs. wheel_only
# 'tw': treadmill vs. wheel_only
comp_groups = c('st','sw','tw')
plot_list = list()

# Nested for loop: main loop specifies whether positive or negative effect sizes are being analysed
# inner for loop specifies which comparison group is being analyzed
for (j in c('positive', 'negative')) {
  print(paste(j, 'effect size regions'))
  for (i in 1:length(comp_groups)) {
    comparison = comp_groups[i] # 'st', 'sw', or 'tw'
    
    # Based on whether we are looking at the positive/negative effect size group, or the 'st', 'sw', or 'tw' comparison group,
    # we access the appropriate CSV in the below line; CSV is the output of the read_posthoc.py script
    top_comp=read.csv(paste('/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Reference Files/User-generated Files/',comparison,'_regs/top_',substr(j,1,3),'_regions_',comparison,'.csv',sep=""))
    # data_comp = data[data$Treatment %in% dict[[i]],]
    
    # Assign values based on CSV columns to significant region names, abbreviations, and p-values and effect sizes
    sig_reg = top_comp$Abbreviation
    reg_struc = top_comp$Structure
    pvals_regs = formatC(top_comp$P.value, format = "e", digits = 2)
    eff_sizes = formatC(top_comp$Effect.Size, format = "e", digits = 2)
    
    graph_list = list()
    for (k in 1:3) {
      if (is.na(sig_reg[k]) == F) {
        res.aov <- aov(get(sig_reg[k]) ~ Treatment, data = data)
        tuk=tukey_hsd(res.aov)
        data_temp <- data[,c("Treatment",sig_reg[k])] %>% setNames(c("treatment","region"))
        tuk <- add_y_position(tuk, data=data_temp, formula=region ~ treatment)
        pbar_tab <- tuk[,c("group1", "group2", "p.adj", "y.position")]
  
        p <- ggplot(data, aes_string(x="Treatment", y=sig_reg[k])) + stat_pvalue_manual(pbar_tab, label = "p.adj", size = 3, tip.length = 0, hide.ns = TRUE) +
          geom_violin() + geom_boxplot(width=0.1) + geom_dotplot(binaxis= "y", stackdir = "center", dotsize=0.5, fill='red') + 
          labs(title=reg_struc[k], subtitle=paste("P-value of ",toString(pvals_regs[k])," | Effect size of ",toString(eff_sizes[k])), y="", x="") + theme_bw() +
          theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
        
        graph_list[[k]] = p
      }
      else {
        p <- ggplot(data, aes_string(x="Treatment", y=sig_reg[length(sig_reg)])) + geom_blank() + theme_bw() + labs(x="", y="")
        graph_list[[k]] = p
      }
    }
    
    p1 <- graph_list[[1]] #+ theme(axis.title.y = element_text(margin = margin(r = 20)), axis.title.x = element_blank())
    p2 <- graph_list[[2]] #+ theme(axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_blank())
    p3 <- graph_list[[3]] #+ theme(axis.title = element_blank())
    
    # Add total figure titles
    patchwork <- p1 + p2 + p3 + plot_annotation(
      title = paste(dict[[comparison]][1],'vs.',dict[[comparison]][2],'Exercise'),
      theme = theme(plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5, face = "bold"))
    )
    plot_list[[i]] = patchwork
    full_plot <- p1 + p2 + p3 + plot_annotation(
      title = paste('Regions of Significance following Post-Hoc Analysis (',str_to_title(j),' Effect Size)',sep=""),
      subtitle = paste(dict[[comparison]][1],'vs.',dict[[comparison]][2],'Exercise'),
      caption = 'DRAFT',
      theme = theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
                    plot.margin=unit(c(1,1,-0.5,0), 'cm'))
    )
    gt <- patchwork::patchworkGrob(full_plot)
    full_plot <- gridExtra::grid.arrange(gt, left = "Normalized Brain Proportion (%)", bottom = "Treatment")
    
    # Saving individual plots
    File <- paste("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Statistical Analysis/Output Figures/",j,"_eff/",comparison,'_',substr(j,1,3),'.png',sep="")
    ggsave(File, plot = full_plot, width=1213, height=514, dpi = 150, units='px', scale=2)
    
    # For tracking plotting and debugging
    print(paste(comparison, 'complete'))
  }
  
  comp_plot = plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], nrow=3, ncol=1)
  composite_figure <- comp_plot + plot_annotation(
    title = paste('Regions of Significance following Post-Hoc Analysis (',str_to_title(j),' Effect Size)',sep=""),
    caption = 'DRAFT',
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5),
                  plot.margin=unit(c(1,1,-0.5,0), 'cm'))
  )
  gt <- patchwork::patchworkGrob(composite_figure)
  composite_figure <- gridExtra::grid.arrange(gt, left = "Normalized Brain Proportion (%)", bottom = "Treatment")
  # Saving Plots
  File <- paste("/Volumes/GoogleDrive/My Drive/Education School/Duke University/Year 4 (2022-2023)/Courses/Semester 1/BME 493 (Badea Independent Study)/Bass-Connections-F22/Statistical Analysis/Output Figures/comparison_",substr(j,1,3),'.png',sep="")
  ggsave(File, plot = composite_figure, width=1322, height=1322, dpi = 150, units='px', scale=2)
}