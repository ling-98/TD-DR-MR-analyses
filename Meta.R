

data<-read.csv("data.csv") #改成自己的路径
View(data)
library(meta)





#数据正态性检验
shapiro.test(data$OR)
#其具体命令为：
lnor<-  log(data[,"OR"])
lnuci<-  log(data[,"UCI"])
lnlci<-  log(data[,"LCI"])
selnor<-  (lnuci-lnlci)/(2*1.96)

data <- apply(data[,1:5],2,as.numeric)  
#随后进行meta分析
MetaOR=metagen(lnor, #效应量的对数
               selnor, #效应量的标准误
               sm="OR",   #可选合并效应量类别"RD", "RR", "OR"
               data=data,
               studlab=paste(data$Items,sep="-"),#添加标签
               subgroup = subgroups,
               fixed = F,
               random = T,
               backtransf=T
)
MetaOR
head(data)
#森林图
forest(MetaOR)
settings.meta('JAMA')#JAMA"RevMan5""IQWiG5""IQWiG6", "geneexpr", or "meta4"(灰色)
settings.meta('RevMan5')

settings.meta("reset") #恢复默认格式
forest(MetaOR,
       label.left = "pre-EST+LC", 
       label.right = "LCBDE+LC",
       colgap.studlab = "0.8cm", 
       colgap.forest.left = "3cm",
       layout = "RevMan5",
       leftcols="Items"
)
MetaOR

P1<-plot(MetaOR, "OR")
