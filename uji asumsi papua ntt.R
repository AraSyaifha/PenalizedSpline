#membaca data
ara=read.csv("papuabaratntt.csv", sep=";", header = TRUE)
ara

#persiapan data 
y= ara$y
x1= ara$x1

#membentuk model 
library(lmtest)
model1 = lm(y~x1, ara)
summary(model1)

#ujiasumsi
error = model1$residuals
error

#ujinormalitas dari residual 
library(nortest)
shapiro.test(error)  #shapirowilk
ad.test(error)
lillie.test(error)  #kolmogrof

#Uji autokorelasi
library(car)
durbinWatsonTest(model1)

#uji heterokedastisitas
abs_error = abs(error)
abs_error
uji_identik=lm(abs_error~x1)
summary(uji_identik)

library(lmtest)
bptest(model1)
