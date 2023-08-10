data=read.table("clipboard",header=T)
View(data)
 str(data)


library("psych")
b<- corr.test(data[1:2])
b$r
b$t 
b$p 
b$se

sink("correlation.doc")
print(b)
sink()

pairs.panels(data[,-5],pch = 35, stars = T)
detach("package:psych",unload = TRUE)

library("PerformanceAnalytics")
require("PerfomarmanceAnalytics")
chart.Correlation(data[1:4], histrogam = TRUE, pch=100)
detach("package:PerfomarmanceAnalytics", unload = TRUE)

