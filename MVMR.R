#设置工作目录
setwd("D:")
#读入数据包
library(TwoSampleMR)
########################IEU在线提取数据MVMR###########################
id_exposure <- c("ieu-a-299","ieu-a-300","ieu-a-302") # 三个暴露分别是HDL cholesterol，LDL cholesterol和Triglycerides
id_outcome <- "ieu-a-7"  #定义结局变量ID
exposure_dat <- mv_extract_exposures(id_exposure)  #三个显著相关的SNP提出
outcome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome)  #结局数据
mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
res <- mv_multiple(mvdat)
res_OR<-generate_odds_ratios(res$result)
res_OR
#write.table(res_OR, file="MV孟德尔随机化结果.xls",sep="\t",quote=F)

##################MVMR所用的SNP是根据什么规则来的######################
exp1 <- extract_instruments(outcomes="ieu-a-299")  
exp2 <- extract_instruments(outcomes="ieu-a-300")  
exp3 <- extract_instruments(outcomes="ieu-a-302")

sum=dim(exp1)[1]+dim(exp2)[1]+dim(exp3)[1]

sum

exp<-rbind(exp1,exp2,exp3)  #按列合并
exp<-exp["SNP"]  #提出素有SNP
exp<-unique(exp)  #去非重复SNP
dim(exp)  #计算最后SNP

####################LASSO去除高度共线MR###############
mv_lasso_feature_selection(mvdat)   #LASSO

mv_lasso<-mv_subset(
  mvdat,
  features = mv_lasso_feature_selection(mvdat),   #LASSO后进行MR
  intercept = FALSE,
  instrument_specific = FALSE,
  pval_threshold = 5e-08,
  plots = FALSE
)
mv_lasso
mv_lasso_OR<-generate_odds_ratios(mv_lasso$result)
mv_lasso_OR
####################单个SNP的MVMR################
mv_residual(
  mvdat,
  intercept = FALSE,
  instrument_specific = FALSE,
  pval_threshold = 5e-08,
  plots = FALSE
)
#######################PRESSO####################
# mr_presso可以完成多变量
#if (!require("devtools")) { install.packages("devtools") } else {}
#devtools::install_github("rondolab/MR-PRESSO")

SummaryStats<-cbind(mvdat[["outcome_beta"]],
                    mvdat[["exposure_beta"]][,1],
                    mvdat[["exposure_beta"]][,2],
                    mvdat[["exposure_beta"]][,3],
                    mvdat[["exposure_se"]][,1],
                    mvdat[["exposure_se"]][,2],
                    mvdat[["exposure_se"]][,3],
                    mvdat[["outcome_se"]])
SummaryStats<-data.frame(SummaryStats)

library(MRPRESSO)
mr_presso(BetaOutcome = "X1",
          BetaExposure = c("X2", "X3","X4"), 
          SdOutcome = "X8", 
          SdExposure = c("X5", "X6","X7"),
          OUTLIERtest = TRUE, 
          DISTORTIONtest = TRUE, 
          data = SummaryStats,
          NbDistribution = 1000, 
          SignifThreshold = 0.05)
#################MendelianRandomization包的几种方法############

#if (!requireNamespace("MendelianRandomization"))
#  install.packages("MendelianRandomization")
library(MendelianRandomization)
#读取实例数据
#准备数据
MRMVInputObject <- mr_mvinput(bx = cbind(ldlc, hdlc, trig),
                              bxse = cbind(ldlcse, hdlcse, trigse),
                              by = chdlodds, 
                              byse = chdloddsse)

MRMVInputObject

MRMVInputObject_1<- mr_mvinput(bx = cbind(SummaryStats$X2,SummaryStats$X3,SummaryStats$X4),
                              bxse = cbind(SummaryStats$X5,SummaryStats$X6,SummaryStats$X7),
                              by = SummaryStats$X1, 
                              byse = SummaryStats$X8)
MRMVInputObject_1
#IVW方法
MRMVObject <- mr_mvivw(MRMVInputObject, 
                       model = "default",
                       correl = FALSE,
                       distribution = "normal",
                       alpha = 0.05)

MRMVObject

MRMVObject <- mr_mvivw(MRMVInputObject_1, 
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

MRMVObject<-mr_mvegger(
                      MRMVInputObject_1,
                      orientate = 1,
                      correl = FALSE,
                      distribution = "normal",
                      alpha = 0.05)
MRMVObject

#LASSO
MRMVObject<-mr_mvlasso(
  MRMVInputObject,
  orientate = 1,
  distribution = "normal",
  alpha = 0.05,
  lambda = numeric(0)
)
MRMVObject

MRMVObject<-mr_mvlasso(
  MRMVInputObject_1,
  orientate = 1,
  distribution = "normal",
  alpha = 0.05,
  lambda = numeric(0)
)
MRMVObject

#median
MRMVObject<-mr_mvmedian(
  MRMVInputObject,
  distribution = "normal",
  alpha = 0.05,
  iterations = 10000,
  seed = 314159265
)
MRMVObject

MRMVObject<-mr_mvmedian(
  MRMVInputObject_1,
  distribution = "normal",
  alpha = 0.05,
  iterations = 10000,
  seed = 314159265
)
MRMVObject

####################RMVMR###########################
#remotes::install_github("WSpiller/RMVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
library(MVMR)
#rawdat_mvmr<-rawdat_rmvmr
head(rawdat_mvmr)

F.data <- format_mvmr(BXGs = rawdat_mvmr[,c(1,2,3)],
                      BYG = rawdat_mvmr[,7],
                      seBXGs = rawdat_mvmr[,c(4,5,6)],
                      seBYG = rawdat_mvmr[,8],
                      RSID = rawdat_mvmr[,9])
head(F.data)

sres <- strength_mvmr(r_input = F.data, gencov = 0)

pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)

res <- ivw_mvmr(r_input = F.data)












#下边给大家放了TwoSampleMR本地的代码，个人感觉很难用，不灵活，有兴趣可以研究
###################本地暴露读取#################
exp_l<-mv_extract_exposures_local(
 c("exp_a.txt","exp_b.txt"),
  sep = " ",
  phenotype_col = "exposure",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  id_col = "id.exposure"
)

outcome_dat1<-format_data(
  a,
  type = "outcome",
  snps = exp$"SNP",
  header = TRUE,
  snp_col = "SNP",
  beta_col = "BETA_INSOMNIA",
  se_col = "SE_INSOMNIA",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_INSOMNIA",
)



