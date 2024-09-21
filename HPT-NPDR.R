library(TwoSampleMR)
library(ggplot2)
library(RadialMR)


#reading outcome
NPDR_out <- read_outcome_data(
  snps = HPT_exp1$SNP,
  filename = "finngen_R7_DM_BCKGRND_RETINA_NONPROLIF.txt",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval")


#rs1000263, rs11694732, rs1464085, rs7144089, rs7754251
HPT_NPDR <- harmonise_data(
  exposure_dat =  HPT_exp1, 
  outcome_dat = NPDR_out
)

#excluding the snp associated with the outcome.
HPT_NPDR1<- subset(HPT_NPDR,pval.outcome>=5e-08)
dim(HPT_NPDR1)


#checking the outliers:
library("MRPRESSO")
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data =HPT_NPDR3, NbDistribution = 4500,  
          SignifThreshold = 0.05, seed=1234)
HPT_NPDR1的mrpresso的结果：a=0.05（1个离群值）,p<0.000833333333333333，Distortion Pvalue 0.77
HPT_NPDR2<- subset(HPT_NPDR1,SNP != "rs62404122")

HPT_NPDR2的mrpresso的结果：a=0.05（1个离群值）,p <0.000222222222222222，Distortion Pvalue 0.8433333
HPT_NPDR3<- subset(HPT_NPDR2,SNP != "rs385251")
HPT_NPDR3的mrpresso的结果：a=0.05（0个离群值）,p 0.001333333，Distortion Pvalue NA

ivw.object<-ivw_radial(r_input = HPT_NPDR3,alpha = 0.05,weights = 1,tol=0.0001,summary = TRUE)
plotly_radial(ivw.object)（直线）16个离群值
plot_radial(ivw.object)（曲线）
HPT_NPDR4<- subset(HPT_NPDR3,SNP != "rs3184504"&SNP != "rs56175143"&SNP != "rs193477"
                   &SNP != "rs661891"&SNP != "rs2893970"&SNP != "rs12540388"
                   &SNP != "rs724078"&SNP != "rs4912068"&SNP != "rs147791594"
                   &SNP != "rs10818050"&SNP != "rs2046045"&SNP != "rs2275710"
                   &SNP != "rs13137589"&SNP != "rs11969311"&SNP != "rs10055404"
                   &SNP != "rs58451984")

write.csv(HPT_NPDR, file="HPT_NPDR.csv")
write.csv(HPT_NPDR1, file="HPT_NPDR1.csv")
write.csv(HPT_NPDR2, file="HPT_NPDR2.csv")
write.csv(HPT_NPDR3, file="HPT_NPDR3.csv")
write.csv(HPT_NPDR4, file="HPT_NPDR4.csv")

#MR analysis
#After first checking the outliers
result<-mr(HPT_NPDR2)
result
OR<-generate_odds_ratios(result)
write.csv(OR, file="OR.csv")

#Second
result1<-mr(HPT_NPDR3)
result1
OR1<-generate_odds_ratios(result1)
write.csv(OR1, file="OR1.csv")

#After checking the outliers by radialMR
result2<-mr(HPT_NPDR4)
result2
OR2<-generate_odds_ratios(result2)
write.csv(OR2, file="OR2.csv")


#directional test
首先计算结局的R2
HPT_NPDR4$R2.outcome <- 2*(1-HPT_NPDR4$eaf.outcome)*HPT_NPDR4$eaf.outcome*(HPT_NPDR4$beta.outcome)^2
HPT_NPDR4$F.outcome<- (HPT_NPDR4$R2.outcome)/(1-HPT_NPDR4$R2.outcome)*(297584-2)
write.csv(HPT_NPDR3, file="HPT_NPDR3.csv")

out <- directionality_test(HPT_NPDR4) 
knitr::kable(out) TRUE

#pleiotropy
pleio <- mr_pleiotropy_test(HPT_NPDR2)
pleio
write.csv(pleio, file="pleio.csv")

pleio1 <- mr_pleiotropy_test(HPT_NPDR3)
pleio1
write.csv(pleio1, file="pleio1.csv")

pleio2 <- mr_pleiotropy_test(HPT_NPDR4)
pleio2
write.csv(pleio2, file="pleio2.csv")


#heterogeneity test
het <- mr_heterogeneity(HPT_NPDR2)
het
write.csv(het, file="het.csv")

het1 <- mr_heterogeneity(HPT_NPDR3)
het1
write.csv(het1, file="het1.csv")

het2 <- mr_heterogeneity(HPT_NPDR4)
het2
write.csv(het2, file="het2.csv")

#leave-out-one 
leave <- mr_leaveoneout(HPT_NPDR4)
mr_leaveoneout_plot(leave)

#forest plot
sin <- mr_singlesnp(HPT_NPDR4)
p1 <- mr_forest_plot(sin)
p1

#scatter plot
p2 <-mr_scatter_plot(result1, HPT_NPDR4)
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
statistic Power：0.92
