library(TwoSampleMR)
library(ggplot2)
library(RadialMR)



#reading outcome
PDR_out <- read_outcome_data(
  snps = HPT_exp1$SNP,
  filename = "finngen_R9_DM_RETINA_PROLIF.txt",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval")

#rs1000263, rs11694732, rs1464085, rs7144089, rs7754251
HPT_PDR <- harmonise_data(
  exposure_dat =  HPT_exp1, 
  outcome_dat = PDR_out
)

#excluding the snp associated with the outcome.
HPT_PDR1<- subset(HPT_PDR,pval.outcome>=5e-08)
dim(HPT_PDR1)


#checking the outliers:
library("MRPRESSO")
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data =HPT_PDR2, NbDistribution = 4500,  
          SignifThreshold = 0.05, seed=1234)
HPT_PDR1的mrpresso的结果：a=0.05（2个离群值）,p<0.000222222222222222，Distortion Pvalue 0.3297778
HPT_PDR2<- subset(HPT_PDR1,SNP != "rs10818050"&SNP != "rs9348859")

HPT_PDR2的mrpresso的结果：a=0.05（0个离群值）,p <0.000222222222222222，Distortion Pvalue NA


ivw.object<-ivw_radial(r_input = HPT_PDR2,alpha = 0.05,weights = 1,tol=0.0001,summary = TRUE)
plotly_radial(ivw.object)（直线）19个离群值
plot_radial(ivw.object)（曲线）
HPT_PDR3<- subset(HPT_PDR2,SNP != "rs3184504"&SNP != "rs56175143"&SNP != "rs61938962"
                  &SNP != "rs9368744"&SNP != "rs2774290"&SNP != "rs11799926"
                  &SNP != "rs212388"&SNP != "rs1808192"&SNP != "rs724078"
                  &SNP != "rs71641308"&SNP != "rs2046045"&SNP != "rs3130186"
                  &SNP != "rs753760"&SNP != "rs12897126"&SNP != "rs61916675"
                  &SNP != "rs66760320"&SNP != "rs927984"&SNP != "rs7863943"
                  &SNP != "rs146527047")
write.csv(HPT_PDR1, file="HPT_PDR1.csv")
write.csv(HPT_PDR2, file="HPT_PDR2.csv")
write.csv(HPT_PDR3, file="HPT_PDR3.csv")


#MR analysis
#After first checking the outliers
result<-mr(HPT_PDR2)
result
OR<-generate_odds_ratios(result)
write.csv(OR, file="OR.csv")
write.csv(result, file="result.csv")

#After checking the outliers by radialMR
result1<-mr(HPT_PDR3)
result1
OR1<-generate_odds_ratios(result1)
write.csv(OR1, file="OR1.csv")
write.csv(result1, file="result1.csv")

#directional test
首先计算结局的R2
HPT_PDR3$R2.outcome <- 2*(1-HPT_PDR3$eaf.outcome)*HPT_PDR3$eaf.outcome*(HPT_PDR3$beta.outcome)^2
HPT_PDR3$F.outcome<- (HPT_PDR3$R2.outcome)/(1-HPT_PDR3$R2.outcome)*(372092-2)
write.csv(HPT_PDR3, file="HPT_PDR3.csv")

out <- directionality_test(HPT_PDR3) 
knitr::kable(out) TRUE

#pleiotropy
pleio <- mr_pleiotropy_test(HPT_PDR2)
pleio
write.csv(pleio, file="pleio.csv")

pleio1 <- mr_pleiotropy_test(HPT_PDR3)
pleio1
write.csv(pleio1, file="pleio1.csv")


#heterogeneity test
het <- mr_heterogeneity(HPT_PDR2)
het
write.csv(het, file="het.csv")

het1 <- mr_heterogeneity(HPT_PDR3)
het1
write.csv(het1, file="het1.csv")

#leave-out-one 
leave <- mr_leaveoneout(HPT_PDR3)
mr_leaveoneout_plot(leave)

#forest plot
sin <- mr_singlesnp(HPT_PDR3)
p1 <- mr_forest_plot(sin)
p1

#scatter plot
p2 <-mr_scatter_plot(result1, HPT_PDR3)
p2
#funnel plot
p3<- mr_funnel_plot(sin)
p3


#searching confounders
library(TwoSampleMR)
library(MendelianRandomization)
confounder=MendelianRandomization::phenoscanner(
  snpquery = unique(HPT_exp$SNP),
  catalogue = "GWAS",pvalue = 5e-08,
  proxies = "None",r2=0.8,build = 37
)

a=confounder$results[,c("snp","trait")];a
a1=dplyr::count(a,trait,sort = TRUE)
write.csv(a1, file="a1.csv")
write.csv(a, file="a.csv")

#counter the power by the website:
statistic Power：1
