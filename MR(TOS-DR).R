library(TwoSampleMR)
library(ggplot2)
library(RadialMR)

#读取暴露
TOS_exp <- read_exposure_data(
  filename = "finngen_R9_THYROTOXICOSIS.txt",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "ref",
  other_allele_col = "alt",
  eaf_col = "af_alt",
  pval_col = "pval"
)
TOS_exp <- TOS_exp[TOS_exp$pval.exposure < 5e-8,]
TOS_exp <- clump_data(TOS_exp,clump_r2 = 0.01,clump_kb = 5000)



#计算F，TOS_exp1是F大于10
TOS_exp$R2.exposure<- 2*(1-TOS_exp$eaf.exposure)*TOS_exp$eaf.exposure*(TOS_exp$beta.exposure)^2
TOS_exp$F.exposure<- (TOS_exp$R2.exposure)/(1-TOS_exp$R2.exposure)*(375580-2)
TOS_exp1<- subset(TOS_exp,F.exposure>10)
dim(TOS_exp1)
write.csv(TOS_exp1, file="TOS_exp1.csv")

#读取结局
DR_out <- read_outcome_data(
  snps = TOS_exp1$SNP,
  filename = "finngen_R9_DM_RETINOPATHY_EXMORE.txt",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "ref",
  other_allele_col = "alt",
  eaf_col = "af_alt",
  pval_col = "pval")


TOS_DR <- harmonise_data(
  exposure_dat =  TOS_exp1, 
  outcome_dat = DR_out
)

#排除和结局相关的SNP。
TOS_DR1<- subset(TOS_DR,pval.outcome>=5e-08)
dim(TOS_DR1)
write.csv(TOS_DR1, file="TOS_DR1.csv")


library("MRPRESSO")
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data =TOS_DR2, NbDistribution = 1200,  
          SignifThreshold = 0.05, seed=1234)
TOS_DR1的mrpresso的结果：a=0.05（1个离群值）,p<0.000833333333333333，Distortion Pvalue 0.6108333
TOS_DR2<- subset(TOS_DR1,SNP != "rs2160215")

TOS_DR2的mrpresso的结果：a=0.05（0个离群值）,p<0.000833333333333333，Distortion Pvalue NA


ivw.object<-ivw_radial(r_input = TOS_DR2,alpha = 0.05,weights = 1,tol=0.0001,summary = TRUE)
plotly_radial(ivw.object)（直线） 12个离群值
plot_radial(ivw.object)（曲线）
TOS_DR3<- subset(TOS_DR2,SNP != "rs146542010"&SNP != "rs11571297"&SNP != "rs149465829"&SNP != "rs185466530"
                 &SNP != "rs12198492"&SNP != "rs2495994"&SNP != "rs1893592"&SNP != "rs4804413"
                 &SNP != "rs2466074"&SNP != "rs17767904"&SNP != "rs12284404"&SNP != "rs4338740")
write.csv(TOS_DR1, file="TOS_DR1.csv")
write.csv(TOS_DR2, file="TOS_DR2.csv")
write.csv(TOS_DR3, file="TOS_DR3.csv")


#剔除混杂：
TOS_DR4<- subset(TOS_DR3,SNP != "rs6679677"&SNP != "rs11571297"&SNP != "rs12190030"&SNP != "rs2524106"
                  &SNP != "rs28633917"&SNP != "rs2647004"&SNP != "rs3819714"
                  &SNP != "rs707958"&SNP != "rs61734579"&SNP != "rs6065926"&SNP != "rs1893592")
#MR analysis
#After first checking the outliers
TOS_DRresult<-mr(TOS_DR2)
TOS_DRresult
OR<-generate_odds_ratios(TOS_DRresult)
write.csv(OR, file="OR.csv")
write.csv(TOS_DRresult, file="TOS_DRresult.csv")

#After checking the outliers by radialMR
TOS_DRresult1<-mr(TOS_DR3)
TOS_DRresult1
OR1<-generate_odds_ratios(TOS_DRresult1)
write.csv(OR1, file="OR1.csv")
write.csv(TOS_DRresult1, file="TOS_DRresult1.csv")

#剔除混杂后分析：
TOS_DRresult2<-mr(TOS_DR4)
TOS_DRresult2
OR2<-generate_odds_ratios(TOS_DRresult2)
write.csv(OR2, file="OR2.csv")


#directional test
首先计算结局的R2
TOS_DR3$R2.outcome <- 2*(1-TOS_DR3$eaf.outcome)*TOS_DR3$eaf.outcome*(TOS_DR3$beta.outcome)^2
TOS_DR3$F.outcome<- (TOS_DR3$R2.outcome)/(1-TOS_DR3$R2.outcome)*(319046-2)
write.csv(TOS_DR3, file="TOS_DR3.csv")

out <- directionality_test(TOS_DR3) 
knitr::kable(out) TRUE

多效性
TOS_DR_pleio <- mr_pleiotropy_test(TOS_DR2)
TOS_DR_pleio
write.csv(TOS_DR_pleio, file="TOS_DR_pleio.csv")

TOS_DR_pleio1 <- mr_pleiotropy_test(TOS_DR3)
TOS_DR_pleio1
write.csv(TOS_DR_pleio1, file="TOS_DR_pleio1.csv")

TOS_DR_pleio2 <- mr_pleiotropy_test(TOS_DR4)
TOS_DR_pleio2
write.csv(TOS_DR_pleio2, file="TOS_DR_pleio2.csv")


异质性检验
TOS_DR_het <- mr_heterogeneity(TOS_DR2)
TOS_DR_het
write.csv(TOS_DR_het, file="TOS_DR_het.csv")

TOS_DR_het1 <- mr_heterogeneity(TOS_DR3)
TOS_DR_het1
write.csv(TOS_DR_het1, file="TOS_DR_het1.csv")

TOS_DR_het2 <- mr_heterogeneity(TOS_DR4)
TOS_DR_het2
write.csv(TOS_DR_het2, file="TOS_DR_het2.csv")

留一法
TOS_DR_leave <- mr_leaveoneout(TOS_DR3)
mr_leaveoneout_plot(TOS_DR_leave)
TOS_DR_leave1 <- mr_leaveoneout(TOS_DR4)
mr_leaveoneout_plot(TOS_DR_leave1)

森林图
TOS_DR_sin <- mr_singlesnp(TOS_DR3)
p1 <- mr_forest_plot(TOS_DR_sin)
p1

散点图
p2 <-mr_scatter_plot(TOS_DRresult1, TOS_DR3)
p2
漏斗图
p3<- mr_funnel_plot(TOS_DR_sin)
p3


#searching confounders
library(TwoSampleMR)
library(MendelianRandomization)
confounder=MendelianRandomization::phenoscanner(
  snpquery = unique(TOS_exp1$SNP),
  catalogue = "GWAS",pvalue = 5e-08,
  proxies = "None",r2=0.8,build = 37
)

a=confounder$results[,c("snp","trait")];a
a1=dplyr::count(a,trait,sort = TRUE)
write.csv(a1, file="a1.csv")
write.csv(a, file="a.csv")

使用网站计算Power：0.97
