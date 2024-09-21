library(TwoSampleMR)
library(ggplot2)
library(RadialMR)



#读取结局
DR_out <- read_outcome_data(
  snps = HPT_exp1$SNP,
  filename = "finngen_R9_DM_RETINOPATHY_EXMORE.txt",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval")


#rs1000263, rs11694732, rs1464085, rs7144089, rs7754251
HPT_DR <- harmonise_data(
  exposure_dat =  HPT_exp1, 
  outcome_dat = DR_out
)

#排除和结局相关的SNP。
HPT_DR1<- subset(HPT_DR,pval.outcome>=5e-08)
dim(HPT_DR1)


library("MRPRESSO")
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data =HPT_DR2, NbDistribution = 4500,  
          SignifThreshold = 0.05, seed=1234)
HPT_DR1的mrpresso的结果：a=0.05,p<0.000222222222222222，Distortion Pvalue 0.3982222
HPT_DR2<- subset(HPT_DR1,SNP != "rs10818050"&SNP != "rs1317983"&SNP != "rs147791594"
                 &SNP != "rs2893970"&SNP != "rs61938962"&SNP != "rs706779"
                 &SNP != "rs9348859")

HPT_DR2的mrpresso的结果：a=0.05（0个离群值）,p <0.000222222222222222，Distortion Pvalue NA


ivw.object<-ivw_radial(r_input = HPT_DR2,alpha = 0.05,weights = 1,tol=0.0001,summary = TRUE)
plotly_radial(ivw.object)（直线） 43个离群值
plot_radial(ivw.object)（曲线）
HPT_DR3<- subset(HPT_DR2,SNP != "rs11571297"&SNP != "rs9267873"&SNP != "rs10814915"
                 &SNP != "rs193477"&SNP != "rs9368744"&SNP != "rs74745605"
                 &SNP != "rs4409785"&SNP != "rs1990760"&SNP != "rs6448432"
                 &SNP != "rs1057373"&SNP != "rs9497965"&SNP != "rs111937080"
                 &SNP != "rs75329808"&SNP != "rs724078"&SNP != "rs11969311"
                 &SNP != "rs1885013"&SNP != "rs12736474"&SNP != "rs72729322"
                 &SNP != "rs2387397"&SNP != "rs140367581"&SNP != "rs9378805"
                 &SNP != "rs1788232"&SNP != "rs212388"&SNP != "rs12697352"&SNP != "rs3103991"&SNP != "rs6831973"
                 &SNP != "rs2046045"&SNP != "rs3130186"&SNP != "rs7746336"
                 &SNP != "rs753760"&SNP != "rs17364832"&SNP != "rs71430783"
                 &SNP != "rs2275710"&SNP != "rs244687"&SNP != "rs7902146"
                 &SNP != "rs12897126"&SNP != "rs13137589"&SNP != "rs76428106"
                 &SNP != "rs8006310"&SNP != "rs568999"&SNP != "rs927984"
                 &SNP != "rs61916675"&SNP != "rs35937663")
write.csv(HPT_DR, file="HPT_DR.csv")
write.csv(HPT_DR1, file="HPT_DR1.csv")
write.csv(HPT_DR2, file="HPT_DR2.csv")
write.csv(HPT_DR3, file="HPT_DR3.csv")

#MR analysis
#After first checking the outliers
result<-mr(HPT_DR2)
result
OR<-generate_odds_ratios(result)
write.csv(OR, file="OR.csv")
write.csv(result, file="result.csv")

#After checking the outliers by radialMR
result1<-mr(HPT_DR3)
result1
OR1<-generate_odds_ratios(result1)
write.csv(OR1, file="OR1.csv")
write.csv(result1, file="result1.csv")

#directional test
首先计算结局的R2
HPT_DR3$R2.outcome <- 2*(1-HPT_DR3$eaf.outcome)*HPT_DR3$eaf.outcome*(HPT_DR3$beta.outcome)^2
HPT_DR3$F.outcome<- (HPT_DR3$R2.outcome)/(1-HPT_DR3$R2.outcome)*(319046-2)
write.csv(HPT_DR3, file="HPT_DR3.csv")

out <- directionality_test(HPT_DR3) 
knitr::kable(out) TRUE

#pleiotropy
pleio <- mr_pleiotropy_test(HPT_DR2)
pleio
write.csv(pleio, file="pleio.csv")

pleio1 <- mr_pleiotropy_test(HPT_DR3)
pleio1
write.csv(pleio1, file="pleio1.csv")


#heterogeneity test
het <- mr_heterogeneity(HPT_DR2)
het
write.csv(het, file="het.csv")

het1 <- mr_heterogeneity(HPT_DR3)
het1
write.csv(het1, file="het1.csv")

#leave-out-one 
leave <- mr_leaveoneout(HPT_DR3)
mr_leaveoneout_plot(leave)

#forest plot
sin <- mr_singlesnp(HPT_DR3)
p1 <- mr_forest_plot(sin)
p1

#scatter plot
p2 <-mr_scatter_plot(result1, HPT_DR3)
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
