library(TwoSampleMR)
library(ggplot2)
library(RadialMR)

#reading exposure
TOS_exp <- read_exposure_data(
  filename = "finngen_R9_THYROTOXICOSIS.txt",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
   effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval"
)
TOS_exp <- TOS_exp[TOS_exp$pval.exposure < 5e-8,]
TOS_exp <- clump_data(TOS_exp,clump_r2 = 0.01,clump_kb = 5000)



#countering the F statistic，TOS_exp1 is greater than 10.
TOS_exp$R2.exposure<- 2*(1-TOS_exp$eaf.exposure)*TOS_exp$eaf.exposure*(TOS_exp$beta.exposure)^2
TOS_exp$F.exposure<- (TOS_exp$R2.exposure)/(1-TOS_exp$R2.exposure)*(375580-2)
TOS_exp1<- subset(TOS_exp,F.exposure>10)
dim(TOS_exp1)
write.csv(TOS_exp1, file="TOS_exp1.csv")

#reading outcome
NPDR_out <- read_outcome_data(
  snps = TOS_exp1$SNP,
  filename = "finngen_R7_DM_BCKGRND_RETINA_NONPROLIF.txt",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval")


TOS_NPDR <- harmonise_data(
  exposure_dat =  TOS_exp1, 
  outcome_dat = NPDR_out
)

#excluding the snp associated with the outcome.
TOS_NPDR1<- subset(TOS_NPDR,pval.outcome>=5e-08)
dim(TOS_NPDR1)


#checking the outliers:
library("MRPRESSO")
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data =TOS_NPDR2, NbDistribution = 1200,  
          SignifThreshold = 0.05, seed=1234)
TOS_NPDR1的mrpresso的结果：a=0.05（2个离群值）,p<0.000833333333333333，Distortion Pvalue 0.19
TOS_NPDR2<- subset(TOS_NPDR1,SNP != "rs3819714"&SNP != "rs707958")

TOS_NPDR2的mrpresso的结果：a=0.05（0个离群值）,p 0.0025，Distortion Pvalue NA


ivw.object<-ivw_radial(r_input = TOS_NPDR2,alpha = 0.05,weights = 1,tol=0.0001,summary = TRUE)
plotly_radial(ivw.object)（直线）5个离群值
plot_radial(ivw.object)（曲线）
TOS_NPDR3<- subset(TOS_NPDR2,SNP != "rs146542010"&SNP != "rs28360634"&SNP != "rs12198492"&SNP != "rs937033"
                  &SNP != "rs28633917")
write.csv(TOS_NPDR1, file="TOS_NPDR1.csv")
write.csv(TOS_NPDR2, file="TOS_NPDR2.csv")
write.csv(TOS_NPDR3, file="TOS_NPDR3.csv")

#剔除混杂：
TOS_NPDR4<- subset(TOS_NPDR3,SNP != "rs6679677"&SNP != "rs11571297"&SNP != "rs12190030"&SNP != "rs2524106"
                  &SNP != "rs28633917"&SNP != "rs2647004"&SNP != "rs3819714"
                  &SNP != "rs707958"&SNP != "rs61734579"&SNP != "rs6065926"&SNP != "rs1893592")

#MR analysis
#After first checking the outliers
TOS_NPDRresult<-mr(TOS_NPDR2)
TOS_NPDRresult
OR<-generate_odds_ratios(TOS_NPDRresult)
write.csv(OR, file="OR.csv")
write.csv(TOS_NPDRresult, file="TOS_NPDRresult.csv")

#After checking the outliers by radialMR
TOS_NPDRresult1<-mr(TOS_NPDR3)
TOS_NPDRresult1
OR1<-generate_odds_ratios(TOS_NPDRresult1)
write.csv(OR1, file="OR1.csv")
write.csv(TOS_NPDRresult1, file="TOS_NPDRresult1.csv")


#剔除混杂后分析：
result2<-mr(TOS_NPDR4)
result2
OR2<-generate_odds_ratios(result2)
write.csv(OR2, file="OR2.csv")

#directional test
首先计算结局的R2
TOS_NPDR3$R2.outcome <- 2*(1-TOS_NPDR3$eaf.outcome)*TOS_NPDR3$eaf.outcome*(TOS_NPDR3$beta.outcome)^2
TOS_NPDR3$F.outcome<- (TOS_NPDR3$R2.outcome)/(1-TOS_NPDR3$R2.outcome)*(297584-2)
write.csv(TOS_NPDR3, file="TOS_NPDR3.csv")

out <- directionality_test(TOS_NPDR3) 
knitr::kable(out) TRUE

#pleiotropy
TOS_NPDR_pleio <- mr_pleiotropy_test(TOS_NPDR2)
TOS_NPDR_pleio
write.csv(TOS_NPDR_pleio, file="TOS_NPDR_pleio.csv")

TOS_NPDR_pleio1 <- mr_pleiotropy_test(TOS_NPDR3)
TOS_NPDR_pleio1
write.csv(TOS_NPDR_pleio1, file="TOS_NPDR_pleio1.csv")


#heterogeneity test
TOS_NPDR_het <- mr_heterogeneity(TOS_NPDR2)
TOS_NPDR_het
write.csv(TOS_NPDR_het, file="TOS_NPDR_het.csv")

TOS_NPDR_het1 <- mr_heterogeneity(TOS_NPDR3)
TOS_NPDR_het1
write.csv(TOS_NPDR_het1, file="TOS_NPDR_het1.csv")

#leave-out-one 
TOS_NPDR_leave <- mr_leaveoneout(TOS_NPDR3)
mr_leaveoneout_plot(TOS_NPDR_leave)

#forest plot
TOS_NPDR_sin <- mr_singlesnp(TOS_NPDR3)
p1 <- mr_forest_plot(TOS_NPDR_sin)
p1

#scatter plot
p2 <-mr_scatter_plot(TOS_NPDRresult1, TOS_NPDR3)
p2
#funnel plot
p3<- mr_funnel_plot(TOS_NPDR_sin)
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

#counter the power by the website:
statistic Power：0.99
