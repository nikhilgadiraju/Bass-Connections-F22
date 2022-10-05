######## this says volume but in fact it reads and do analysis on fa

path_vol="/Users/ali/Desktop/aug/apoe2_paper_stats_volume/individual_label_statistics/"
file_list=list.files(path_vol)

temp=read.delim( paste0(path_vol,file_list[1]) )
len=length(temp$fa_mean)
volumes=matrix(NA,length(file_list),len)
 noreadcsf=c(148,152,161,314,318,327)

for (i in 1:length(file_list)) {
  #print(file)
  temp=read.delim( paste0(path_vol,file_list[i]) )
  # temp$volume_mm3
  #print(sum(temp$volume_mm3))
  volumes[i,2:len]=temp$fa_mean[2:len]
  # whole_volume[i,2]=sum(temp$volume_mm3[-noreadcsf])
  volumes[i,1]=substr( file_list[i] , 1, 6)
}
 volumes=as.data.frame(volumes);
 volumes$V1=as.numeric(substr(volumes$V1,2,6)) # make dwi numeric
 # volumes[,2:len]=as.numeric(volumes[,2:len])


library(readxl)
path_master='/Users/ali/Desktop/aug/MasterSheet_Experiments2021.xlsx'
data=read_xlsx(path_master, sheet = '18ABB11_readable02.22.22_BJ_Cor' )
datatemp=data%>%dplyr::select(DWI,Genotype,Weight, Sex, Diet, Age_Months)#subselect
#nchar(datatemp[111,1])
datatemp=na.omit(datatemp)
datatemp[nchar(datatemp$DWI)==1,]=matrix(NA,1,dim(datatemp)[2])
datatemp=na.omit(datatemp)
datatemp[substr(datatemp$DWI,1,1)!="N",]=matrix(NA,1,dim(datatemp)[2])
datatemp=na.omit(datatemp) ## ommit all na and zero character dwi and died durring
datatemp$DWI=as.numeric(substr(datatemp$DWI,2,6)) # make dwi numeric
datatemp=datatemp[datatemp$Genotype=="APOE22",]

datatemp$DWI%in%volumes$V1
datatemp$DWI[1]==volumes$V1
indeces_of_whole=match(datatemp$DWI, volumes$V1)
temp_bind=cbind(volumes[indeces_of_whole,],datatemp)
temp_bind=as.data.frame(temp_bind)
temp_bind=na.omit(temp_bind)

temp_bind=temp_bind%>%mutate( age_cat=case_when(  Age_Months<median(Age_Months)~1 ,
                                                  Age_Months>=median(Age_Months)~2            ) )
library(emmeans)
library(effectsize)

result=matrix(NA,6,(len-1))
for (i in 2:len) {
  var2=temp_bind[,i]
lm <- lm( var2 ~ as.factor(age_cat)*as.factor(Sex),data=temp_bind )
# lm <- lm( var2 ~ as.numeric(Age_Months)*as.factor(Sex),data=temp_bind )

an=anova(lm)
pvals=an$`Pr(>F)`
eff=eta_squared(lm, partial = TRUE)


result[1:3,(i-1)]=pvals[1:3]
 result[4:6,(i-1)]=eff$Eta2_partial


#emmeans(lm, ~ Age_Months, contr="tukey") 
# sd(temp_bind$V2[temp_bind$age_cat==1]) 
# 
# mean(temp_bind$V2[  temp_bind$age_cat==1  ])
# median(temp_bind$V2[  temp_bind$age_cat==1  ])
# mean(temp_bind$V2[  temp_bind$age_cat==2  ])
# median(temp_bind$V2[  temp_bind$age_cat==2 ])
# sd(temp_bind$Age_Months[temp_bind$age_cat==1])
# sd(temp_bind$Age_Months[temp_bind$age_cat==2])
}

rownames(result)=c("Age", "Sex", "Age*Sex", "Age ES", "Sex ES", "Age*Sex ES")

adjustresult=result
for (j in 1:3) {
  adjustresult[j,]=p.adjust(adjustresult[j,], method = "fdr")
  
}

sum(adjustresult[1,]<0.05) #age
sum(adjustresult[2,]<0.05) #sex
sum(adjustresult[3,]<0.05) #age*sex


pathnames='/Users/ali/Desktop/Jul/apoe/mouse_anatomy.csv'
ROI=read.csv(pathnames, header = TRUE, sep = ",", quote = "")
ROI=paste0(ROI$Bigpart, " ",ROI$ROI)


age_index_sig=which(adjustresult[1,] <0.05)
age_index_sig=setdiff(age_index_sig, noreadcsf)
sig_result_age=adjustresult[,age_index_sig]
colnames(sig_result_age)=age_index_sig
sig_result_age=sig_result_age[,order(sig_result_age[4,], decreasing = T)]


library(xlsx)
write.xlsx2(0, "fa.xlsx", sheetName = "0" )

table_age=matrix(NA, length(age_index_sig) ,12 )
colnames(table_age) = c("Number" , "Index", "Name of the region" , "FDR corrected Pvalue", "Effect Size Eta^2", 
                        "CI lower bound", "CI upper bound", "Mean group 1", "Mean group 2", 
                        "SD group 1", "SD group 2", "F-value")
for (i in 1:length(age_index_sig)) {
  table_age[i,1]=i
  index=as.numeric(colnames(sig_result_age)[i])
  table_age[i,2] = index
  table_age[i,3] = ROI[index]
 
  var2=as.numeric(temp_bind[,index+1])
  lm <- lm( as.numeric(var2) ~ as.factor(age_cat)*as.factor(Sex),data=temp_bind )
 
  # pvals=an$`Pr(>F)`
  temp=sig_result_age[,i]
  eff=eta_squared(lm, partial = TRUE)
  print=cbind(temp[1],t(unlist(eff[1,2:5]))  ) 
  table_age[i,4:7] = print[-c(3)]
  means=by(var2,as.factor(temp_bind$age_cat), mean )
  table_age[i,8:9]=c(means[1], means[2])
  sds=by(var2,as.factor(temp_bind$age_cat), sd )
  table_age[i,10:11]=c(sds[1], sds[2])
  an=anova(lm)
  table_age[i,12]=an$`F value`[1]
}
write.xlsx2(table_age , "fa.xlsx", sheetName =  paste0("age"), append=TRUE )





sex_index_sig=which(adjustresult[2,] <0.05)
sex_index_sig=setdiff(sex_index_sig, noreadcsf)
sig_result_sex=adjustresult[,sex_index_sig]
colnames(sig_result_sex)=sex_index_sig
sig_result_sex=sig_result_sex[,order(sig_result_sex[5,], decreasing = T)]

table_sex=matrix(NA, length(sex_index_sig) ,12 )
colnames(table_sex) = c("Number" , "Index", "Name of the region" , "FDR corrected Pvalue", "Effect Size Eta^2", 
                        "CI lower bound", "CI upper bound", "Mean group 1", "Mean group 2", 
                        "SD group 1", "SD group 2", "F-value")
for (i in 1:length(sex_index_sig)) {
  table_sex[i,1]=i
  index=as.numeric(colnames(sig_result_sex)[i])
  table_sex[i,2] = index
  table_sex[i,3] = ROI[index]
  
  var2=as.numeric(temp_bind[,index+1])
  lm <- lm( as.numeric(var2) ~ as.factor(age_cat)*as.factor(Sex),data=temp_bind )
  
  # pvals=an$`Pr(>F)`
  temp=sig_result_sex[,i]
  eff=eta_squared(lm, partial = TRUE)
  print=cbind(temp[2],t(unlist(eff[2,2:5]))  ) 
  table_sex[i,4:7] = print[-c(3)]
  means=by(var2,as.factor(temp_bind$Sex), mean )
  table_sex[i,8:9]=c(means[1], means[2])
  sds=by(var2,as.factor(temp_bind$Sex), sd )
  table_sex[i,10:11]=c(sds[1], sds[2])
  an=anova(lm)
  table_sex[i,12]=an$`F value`[2]
}
write.xlsx2(table_sex , "fa.xlsx", sheetName =  paste0("Sex"), append=TRUE )







agesex_index_sig=which(adjustresult[3,] <0.05)
agesex_index_sig=setdiff(agesex_index_sig, noreadcsf)
sig_result_agesex=adjustresult[,agesex_index_sig]
colnames(sig_result_agesex)=agesex_index_sig
sig_result_agesex=sig_result_agesex[,order(sig_result_agesex[6,], decreasing = T)]




table_agesex=matrix(NA, length(agesex_index_sig) ,8 )
colnames(table_agesex) = c("Number" , "Index", "Name of the region" , "FDR corrected Pvalue", "Effect Size Eta^2", 
                        "CI lower bound", "CI upper bound", "Mean group 1", "Mean group 2", 
                        "SD group 1", "SD group 2", "F-value")
for (i in 1:length(agesex_index_sig)) {
  table_agesex[i,1]=i
  index=as.numeric(colnames(sig_result_agesex)[i])
  table_agesex[i,2] = index
  table_agesex[i,3] = ROI[index]
  
  var2=as.numeric(temp_bind[,index+1])
  lm <- lm( as.numeric(var2) ~ as.factor(age_cat)*as.factor(Sex),data=temp_bind )
  
  # pvals=an$`Pr(>F)`
  temp=sig_result_agesex[,i]
  eff=eta_squared(lm, partial = TRUE)
  print=cbind(temp[3],t(unlist(eff[3,2:5]))  ) 
  table_agesex[i,4:7] = print[-c(3)]
  # means=by(var2,as.factor(temp_bind$Sex), mean )
  # table_sex[i,8:9]=c(means[1], means[2])
  # sds=by(var2,as.factor(temp_bind$Sex), sd )
  # table_sex[i,10:11]=c(sds[1], sds[2])
  an=anova(lm)
  table_sex[i,8]=an$`F value`[2]
}
write.xlsx2(table_sex , "fa.xlsx", sheetName =  paste0("Age*Sex"), append=TRUE )












# dodge <- position_dodge(width = 1)
# p= ggplot(data=temp_bind, aes(x=as.factor(age_cat), y=V2, fill =Sex, color=Sex, alpha=0.3 ))+
#   geom_violin(inherit.aes=TRUE, position=dodge, alpha=0.3) +
#   scale_color_manual(values=c('blueviolet', 'chartreuse1', 'red'))+
#   scale_fill_manual(values=c('blueviolet', 'chartreuse1', 'red'))+
#   
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=2, alpha=0.95, position=dodge)+
#   geom_boxplot(color="black", outlier.color="black", width=0.2, alpha=.8, position=dodge) +
#   theme_bw()
# plot(p) 
# #dev.off()
# outpath='/Users/ali/Desktop/aug/apoe2_paper_stats_volume'
# ggsave(paste(outpath,'apoe22brainvol_nocsf.png',sep=''), plot = last_plot(), device='png', scale=1, width=4, height=4, unit=c("in"), dpi=200)
