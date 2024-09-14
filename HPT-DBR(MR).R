library(TwoSampleMR)
library(ggplot2)
library(RadialMR)

#reading exposure
HPT_exp <- read_exposure_data(
  filename = "finngen_R9_E4_HYTHY_AI_STRICT.txt",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "ref",
  other_allele_col = "alt",
  eaf_col = "af_alt",
  pval_col = "pval"
)
HPT_exp <- HPT_exp[HPT_exp$pval.exposure < 5e-8,]
HPT_exp <- clump_data(HPT_exp,clump_r2 = 0.01,clump_kb = 5000)



#countering the F statistic，TOS_exp1 is greater than 10.
HPT_exp$R2.exposure<- 2*(1-HPT_exp$eaf.exposure)*HPT_exp$eaf.exposure*(HPT_exp$beta.exposure)^2
HPT_exp$F.exposure<- (HPT_exp$R2.exposure)/(1-HPT_exp$R2.exposure)*(314995-2)
HPT_exp1<- subset(HPT_exp,F.exposure>10)
dim(HPT_exp1)
write.csv(HPT_exp1, file="HPT_exp1.csv")



#reading outcome
DBR_out <- read_outcome_data(
  snps = HPT_exp1$SNP,
  filename = "finngen_R7_DM_BCKGRND_RETINA.txt",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "ref",
  other_allele_col = "alt",
  eaf_col = "af_alt",
  pval_col = "pval")


#s1000263, rs11694732, rs1464085, rs7144089, rs7754251
HPT_DBR <- harmonise_data(
  exposure_dat =  HPT_exp1, 
  outcome_dat = DBR_out
)

#excluding the snp associated with the outcome.
HPT_DBR1<- subset(HPT_DBR,pval.outcome>=5e-08)
dim(HPT_DBR1)


#checking the outliers:
library("MRPRESSO")
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data =HPT_DBR2, NbDistribution = 4500,  
          SignifThreshold = 0.05, seed=1234)
HPT_DBR1的mrpresso的结果：a=0.05（1个离群值）,p<0.000222222222222222，Distortion Pvalue 0.4226667
HPT_DBR2<- subset(HPT_DBR1,SNP != "rs10818050"&SNP != "rs147791594"&SNP != "rs706779"&SNP != "rs9348859")

HPT_DBR2的mrpresso的结果：a=0.05（0个离群值）,p<0.000222222222222222，Distortion Pvalue NA


ivw.object<-ivw_radial(r_input = HPT_DBR2,alpha = 0.05,weights = 1,tol=0.0001,summary = TRUE)
plotly_radial(ivw.object)（直线）38个离群值
plot_radial(ivw.object)（曲线）
HPT_DBR3<- subset(HPT_DBR2,SNP != "rs3184504"&SNP != "rs61938962"&SNP != "rs1990760"
                  &SNP != "rs9368744"&SNP != "rs10814915"&SNP != "rs193477"
                  &SNP != "rs12967678"&SNP != "rs724078"&SNP != "rs75329808"
                  &SNP != "rs9981704"&SNP != "rs12736474"&SNP != "rs2387397"
                  &SNP != "rs9393698"&SNP != "rs2893970"&SNP != "rs9378805"
                  &SNP != "rs2136600"&SNP != "rs28391281"&SNP != "rs142647938"
                  &SNP != "rs72729322"&SNP != "rs2046045"&SNP != "rs3130186"
                  &SNP != "rs1317983"&SNP != "rs11935941"&SNP != "rs10917469"
                  &SNP != "rs17364832"&SNP != "rs10748781"&SNP != "rs753760"
                  &SNP != "rs71430783"&SNP != "rs2275710"&SNP != "rs13137589"
                  &SNP != "rs76428106"&SNP != "rs66760320"&SNP != "rs61916675"
                  &SNP != "rs414755"&SNP != "rs10118880"&SNP != "rs3780446"
                  &SNP != "rs12550413"&SNP != "rs35937663")
write.csv(HPT_DBR, file="HPT_DBR.csv")
write.csv(HPT_DBR1, file="HPT_DBR1.csv")
write.csv(HPT_DBR2, file="HPT_DBR2.csv")
write.csv(HPT_DBR3, file="HPT_DBR3.csv")

#MR analysis
#After first checking the outliers
result<-mr(HPT_DBR2)
result
OR<-generate_odds_ratios(result)
write.csv(OR, file="OR.csv")
write.csv(result, file="result.csv")

#After checking the outliers by radialMR
result1<-mr(HPT_DBR3)
result1
OR1<-generate_odds_ratios(result1)
write.csv(OR1, file="OR1.csv")
write.csv(result1, file="result1.csv")

#directional test
首先计算结局的R2
HPT_DBR3$R2.outcome <- 2*(1-HPT_DBR3$eaf.outcome)*HPT_DBR3$eaf.outcome*(HPT_DBR3$beta.outcome)^2
HPT_DBR3$F.outcome<- (HPT_DBR3$R2.outcome)/(1-HPT_DBR3$R2.outcome)*(300010-2)
write.csv(HPT_DBR3, file="HPT_DBR3.csv")

out <- directionality_test(HPT_DBR3) 
knitr::kable(out) TRUE

#pleiotropy
pleio <- mr_pleiotropy_test(HPT_DBR2)
pleio
write.csv(pleio, file="pleio.csv")

pleio1 <- mr_pleiotropy_test(HPT_DBR3)
pleio1
write.csv(pleio1, file="pleio1.csv")


#heterogeneity test
het <- mr_heterogeneity(HPT_DBR2)
het
write.csv(het, file="het.csv")

het1 <- mr_heterogeneity(HPT_DBR3)
het1
write.csv(het1, file="het1.csv")

#leave-out-one 
leave <- mr_leaveoneout(HPT_DBR3)
mr_leaveoneout_plot(leave)

#forest plot
sin <- mr_singlesnp(HPT_DBR3)
p1 <- mr_forest_plot(sin)
p1

#scatter plot
p2 <-mr_scatter_plot(result1, HPT_DBR3)
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
