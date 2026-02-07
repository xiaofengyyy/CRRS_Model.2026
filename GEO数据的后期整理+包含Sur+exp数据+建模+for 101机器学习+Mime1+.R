EGAS00001004656

EGAS00001005300

GSE138794

GSE141383

GSE162631

GSE173278

GSE182109

GSE223063

GSE235676


https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42127
load("GSE42127_data.RData") #symbol
GPL=GSE42127[[1]]
exp=GSE42127[[2]]
pd =GSE42127[[3]]
library(tidyverse)
colnames(pd)
surdata=select(pd,2,41:47)

table(pd$source_name_ch1)
table(pd$`histology:ch1`)
surdata$OS=ifelse(surdata$`survival status:ch1`=="D",1,0)

#surdata=filter(surdata, OS.time>0)
#surdata=filter(surdata,DFS.time>0)
#surdata=filter(surdata,surdata$`histology:ch1` =="Squamous")
surdata=filter(surdata,surdata$`histology:ch1` =="Adenocarcinoma")

surdata$`histology:ch1` %>% table()

# # 使用strsplit()函数按冒号分割
# #split_data <- strsplit(surdata$characteristics_ch1, ":")
# split_data <- strsplit(surdata$characteristics_ch1, split = "[:;]")
# # 将列表转换为数据框
# split_df <- do.call(rbind, lapply(split_data, function(x) {
#   as.data.frame(t(x))
# }))

# 给新生成的数据框命名列
#colnames(split_df) <- c("Column1", "Column2", "Column3")

# 将分割后的数据框与原始数据框合并
#surdata <- cbind(surdata, split_df)
surdata=select(surdata,1,7,9)
colnames(surdata)=c("ID", "OS.time","OS" )

surdata$OS.time=as.numeric(surdata$OS.time)
surdata$OS.time=surdata$OS.time/30

exp=t(exp) %>% as.data.frame()
surdata=filter(surdata, OS.time>0)

com=intersect(rownames(exp),rownames(surdata))
surdata=surdata[com,]
exp=exp[com,]
identical(rownames(exp),rownames(surdata))


surdata=filter(surdata,OS.time>0)


sur_exp=cbind(surdata,exp)
sur_exp=na.omit(sur_exp)
#sur_exp$OS=ifelse(sur_exp$OS=="Dead",1,0)

#save(sur_exp,file = "GSE42127_sur_exp.Rdata")

save(sur_exp,file = "GSE42127_sur_exp.Squamous.Rdata")
load("GSE42127_sur_exp.Squamous.Rdata")

#Adenocarcinoma
save(sur_exp,file = "GSE42127_sur_exp.Adenocarcinoma.Rdata")
