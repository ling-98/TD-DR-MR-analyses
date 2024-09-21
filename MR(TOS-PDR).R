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
  filename = "finngen_R9_DM_RETINA_PROLIF.txt",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
   effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval")
PDR_out <-(NPDR_out)

TOS_PDR <- harmonise_data(
  exposure_dat =  TOS_exp1, 
  outcome_dat = PDR_out
)

#excluding the snp associated with the outcome.
TOS_PDR1<- subset(TOS_PDR,pval.outcome>=5e-08)
dim(TOS_PDR1)


#checking the outliers:
library("MRPRESSO")
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data =TOS_PDR2, NbDistribution = 1200,  
          SignifThreshold = 0.05, seed=1234)
TOS_PDR1的mrpresso的结果：a=0.05（2个离群值）,p<0.000833333333333333，Distortion Pvalue 0.9541667
TOS_PDR2<- subset(TOS_PDR1,SNP != "rs28360634"&SNP != "rs28633917")

TOS_PDR2的mrpresso的结果：a=0.05（0个离群值）,p 0.0008333333，Distortion Pvalue NA


ivw.object<-ivw_radial(r_input = TOS_PDR2,alpha = 0.05,weights = 1,tol=0.0001,summary = TRUE)
plotly_radial(ivw.object)（直线）8个离群值
plot_radial(ivw.object)（曲线）
TOS_PDR3<- subset(TOS_PDR2,SNP != "rs146542010"&SNP != "rs149465829"&SNP != "rs12198492"&SNP != "rs937033"
                   &SNP != "rs2160215"&SNP != "rs2524106"&SNP != "rs1861628"&SNP != "rs4804413")
write.csv(TOS_PDR1, file="TOS_PDR1.csv")
write.csv(TOS_PDR2, file="TOS_PDR2.csv")
write.csv(TOS_PDR3, file="TOS_PDR3.csv")

#MR analysis
#After first checking the outliers
TOS_PDRresult<-mr(TOS_PDR2)
TOS_PDRresult
OR<-generate_odds_ratios(TOS_PDRresult)
write.csv(OR, file="OR.csv")
write.csv(TOS_PDRresult, file="TOS_PDRresult.csv")

#After checking the outliers by radialMR
TOS_PDRresult1<-mr(TOS_PDR3)
TOS_PDRresult1
OR1<-generate_odds_ratios(TOS_PDRresult1)
write.csv(OR1, file="OR1.csv")
write.csv(TOS_PDRresult1, file="TOS_PDRresult1.csv")

#directional test
#First, countering the R2 of outcome:
TOS_PDR3$R2.outcome <- 2*(1-TOS_PDR3$eaf.outcome)*TOS_PDR3$eaf.outcome*(TOS_PDR3$beta.outcome)^2
TOS_PDR3$F.outcome<- (TOS_PDR3$R2.outcome)/(1-TOS_PDR3$R2.outcome)*(372092-2)
write.csv(TOS_PDR3, file="TOS_PDR3.csv")

out <- directionality_test(TOS_PDR3) 
knitr::kable(out) TRUE

#pleiotropy
TOS_PDR_pleio <- mr_pleiotropy_test(TOS_PDR2)
TOS_PDR_pleio
write.csv(TOS_PDR_pleio, file="TOS_PDR_pleio.csv")

TOS_PDR_pleio1 <- mr_pleiotropy_test(TOS_PDR3)
TOS_PDR_pleio1
write.csv(TOS_PDR_pleio1, file="TOS_PDR_pleio1.csv")


#heterogeneity test
TOS_PDR_het <- mr_heterogeneity(TOS_PDR2)
TOS_PDR_het
write.csv(TOS_PDR_het, file="TOS_PDR_het.csv")

TOS_PDR_het1 <- mr_heterogeneity(TOS_PDR3)
TOS_PDR_het1
write.csv(TOS_PDR_het1, file="TOS_PDR_het1.csv")

#leave-out-one 
TOS_PDR_leave <- mr_leaveoneout(TOS_PDR3)
mr_leaveoneout_plot(TOS_PDR_leave)

#forest plot
TOS_PDR_sin <- mr_singlesnp(TOS_PDR3)
p1 <- mr_forest_plot(TOS_PDR_sin)
p1

#scatter plot
p2 <-mr_scatter_plot(TOS_PDRresult1, TOS_PDR3)
p2
#funnel plot
p3<- mr_funnel_plot(TOS_PDR_sin)
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
