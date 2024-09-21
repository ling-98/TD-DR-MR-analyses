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
GD_out <- read_outcome_data(
  snps = TOS_exp1$SNP,
  filename = "finngen_R9_E4_GRAVES_STRICT.txt",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
   effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval")


TOS_GD <- harmonise_data(
  exposure_dat =  TOS_exp1, 
  outcome_dat = GD_out
)

#excluding the snp associated with the outcome.
TOS_GD1<- subset(TOS_GD,pval.outcome>=5e-08)
dim(TOS_GD1)


#checking the outliers:
library("MRPRESSO")
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data =TOS_GD1, NbDistribution = 1200,  
          SignifThreshold = 0.05, seed=1234)
TOS_GD1的mrpresso的结果：a=0.05（0个离群值）,p 0.05916667，Distortion Pvalue NA



ivw.object<-ivw_radial(r_input = TOS_GD1,alpha = 0.05,weights = 1,tol=0.0001,summary = TRUE)
plotly_radial(ivw.object)（直线）8个离群值
plot_radial(ivw.object)（曲线）
TOS_GD2<- subset(TOS_GD1,SNP != "rs1893592"&SNP != "2046045"&SNP != "2928167")
write.csv(TOS_GD1, file="TOS_GD1.csv")
write.csv(TOS_GD2, file="TOS_GD2.csv")


#MR analysis
#After first checking the outliers
TOS_GDresult<-mr(TOS_GD1)
TOS_GDresult
OR<-generate_odds_ratios(TOS_GDresult)
write.csv(OR, file="OR.csv")
write.csv(TOS_GDresult, file="TOS_GDresult.csv")

#After checking the outliers by radialMR
TOS_GDresult1<-mr(TOS_GD2)
TOS_GDresult1
OR1<-generate_odds_ratios(TOS_GDresult1)
write.csv(OR1, file="OR1.csv")
write.csv(TOS_GDresult1, file="TOS_GDresult1.csv")

#directional test
#First, countering the R2 of outcome:
TOS_GD2$R2.outcome <- 2*(1-TOS_GD2$eaf.outcome)*TOS_GD2$eaf.outcome*(TOS_GD2$beta.outcome)^2
TOS_GD2$F.outcome<- (TOS_GD2$R2.outcome)/(1-TOS_GD2$R2.outcome)*(377277-2)
write.csv(TOS_GD2, file="TOS_GD2.csv")

out <- directionality_test(TOS_GD2) 
knitr::kable(out) TRUE

#pleiotropy
TOS_GD_pleio <- mr_pleiotropy_test(TOS_GD1)
TOS_GD_pleio
write.csv(TOS_GD_pleio, file="TOS_GD_pleio.csv")

TOS_GD_pleio1 <- mr_pleiotropy_test(TOS_GD2)
TOS_GD_pleio1
write.csv(TOS_GD_pleio1, file="TOS_GD_pleio1.csv")


#heterogeneity test
TOS_GD_het <- mr_heterogeneity(TOS_GD1)
TOS_GD_het
write.csv(TOS_GD_het, file="TOS_GD_het.csv")

TOS_GD_het1 <- mr_heterogeneity(TOS_GD2)
TOS_GD_het1
write.csv(TOS_GD_het1, file="TOS_GD_het1.csv")

#leave-out-one 
TOS_GD_leave <- mr_leaveoneout(TOS_GD2)
mr_leaveoneout_plot(TOS_GD_leave)

#forest plot
TOS_GD_sin <- mr_singlesnp(TOS_GD2)
p1 <- mr_forest_plot(TOS_GD_sin)
p1

#scatter plot
p2 <-mr_scatter_plot(TOS_GDresult1, TOS_GD2)
p2
#funnel plot
p3<- mr_funnel_plot(TOS_GD_sin)
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
statistic Power：1
