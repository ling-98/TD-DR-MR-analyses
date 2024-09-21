library(MendelianRandomization)
library(TwoSampleMR)
library(MVMR)
#install.packages("data.table")
library(data.table)


#reading the exposure:
setwd("F:/欧阳欧阳喜欢太阳/甲状腺-DR(MR)/TOS-DR/exposure")
a<-fread("finngen_R9_THYROTOXICOSIS.txt")#TOS
setwd("E:/MR/Greas 病")
b<-fread("finngen_R9_E4_GRAVES_STRICT.txt")#GD
setwd("F:/欧阳欧阳喜欢太阳/甲状腺-DR(MR)/中介孟德尔")
d<-fread("finngen_R9_M13_RHEUMA.txt")#RA



#reading the outcome:

setwd("E:/MR/iron-GWAS data")
c<-fread("finngen_R9_DM_RETINA_PROLIF.txt")

#extracting the snps associated with the exposure:
head(a)
head(b)
head(c)
exp_a <- a[a$pva < 5e-8,]
exp_b <- b[b$pva < 5e-8,] 
exp_d <- d[d$pva < 5e-8,] 

colnames(exp_a)[colnames(exp_a)=="rsids"] <- "SNP" 
colnames(exp_b)[colnames(exp_b)=="rsids"] <- "SNP" 
colnames(exp_d)[colnames(exp_d)=="rsids"] <- "SNP" 

##########第二步，分别去LD############
exp_a=clump_data(exp_a,clump_r2 = 0.01,clump_kb = 5000)
exp_b=clump_data(exp_b,clump_r2 = 0.01,clump_kb = 5000)
exp_d=clump_data(exp_d,clump_r2 = 0.01,clump_kb = 5000)

head(exp_a)
head(exp_b)

#########第三步，去掉重复的SNP######
exp_a<-exp_a[,"SNP"]
exp_b<-exp_b[,"SNP"]
exp_d<-exp_d[,"SNP"]

exp<-rbind(exp_a,exp_b,exp_d)
exp<-unique(exp)   ###去重
exp<-clump_data(exp,clump_r2 = 0.01,clump_kb = 5000)
########第四步，提取所exp的snp信息在暴露1和2，以及结局中############
exp_a1<-merge(exp,a,by.x = "SNP",by.y = "rsids")
exp_b1<-merge(exp,b,by.x = "SNP",by.y = "rsids")
exp_d1<-merge(exp,d,by.x = "SNP",by.y = "rsids")

exp_c<-merge(exp,c,by.x = "SNP",by.y = "rsids")

#rs185466530 (22) rs36229731(36)
exp_c <- subset(exp_c, row.names(exp_c) != 56) 


#读取实例数据
#准备数据
MRMVInputObject <- mr_mvinput(bx = cbind(exp_a1$beta,exp_b1$beta,exp_d1$beta),
                              bxse = cbind(exp_a1$sebeta,exp_b1$sebeta,exp_d1$sebeta),
                              by = exp_c$beta, 
                              byse = exp_c$sebeta)

MRMVInputObject

#IVW方法
MRMVObject <- mr_mvivw(MRMVInputObject, 
                       model = "default",
                       correl = FALSE,
                       distribution = "normal",
                       alpha = 0.05)

MRMVObject


#egger方法
MRMVObject<-mr_mvegger(
  MRMVInputObject,
  orientate = 1,
  correl = FALSE,
  distribution = "normal",
  alpha = 0.05)

MRMVObject

