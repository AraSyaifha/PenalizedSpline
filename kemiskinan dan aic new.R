library(splines)
MPL<-function(x,eps=1e-009)  #mengatasi matriks singular menjadi nonsingular
{
  x<-as.matrix(x)  
  xsvd<-svd(x) #fungsi buat mengatasi
  diago<-xsvd$d[xsvd$d>eps] 
  if(length(diago)==1)
  {
    xplus<-as.matrix(xsvd$v[,1]) %*% t(as.matrix(xsvd$u[,1])/diago)
  }
  else
  {
    xplus<-xsvd$v[,1:length(diago)] %*% diag(1/diago) %*% t(xsvd$u[,1:length(diago)])
  }
  return(xplus)
}

data=read.csv("papuabaratntt.csv",sep=";",header=TRUE)
data
x = data$x1[1:64]
x
y = data$y[1:64]
y
plot(x, type = "l", lwd=2, col="#29648A", ylim = c(13300,16000))

# ynaive = x

aic1<-function(x,y,m,k,lambda)
{
  n=length(y) #panjang data
  print(n)
  aic=10^10 #10 = mendefinisikan nilainya
  t1=seq(min(x),max(x),length.out=100) #seq = membuat vektor
  p=rep(0,(n-2))  #perulangan dari 0 sampai n-2
  for(z in 1:(n-2))
  {
    p[z]=t1[z+1] 
  }
  comb1=combn(p,1,FUN = NULL) #combinasi
  c1=t(comb1) #tranpose
  #for(m in 2)
    for (m in 2:4) #orde 2-4
  {
    k1=c1[,1] #kolom = kanan 
    K1=length(k1) 
    for (v1 in 1:K1) { #letak titik knot 1 - K titik knot
      bs1=bs(x, df=NULL, knot=k1[v1], degre=m-1, intercept=TRUE, Boundary.knots=range(x))
      B = cbind(bs1)
      D = diag(ncol(bs1)) #membuat matriks D
      BtB = t(B) %*% B
      deter=det(BtB) 
      z = MPL(BtB)
      beta= MPL (BtB + lambda * t(D) %*% D) %*% (t(B) %*% y)
      S= B %*% MPL (BtB + lambda * t(D) %*% D) %*% t(B)
      yhat= B %*% beta
      MSE = t(y- yhat) %*% (y-yhat)/n
      I <- matrix(0, ncol = n, nrow = n)
      for (v in 1:n) 
        I[v,v]=1
      k = sum(diag(I-S))
      AIC = (n+((n*log(2*pi))+n*(log(MSE))+(2*(1+m))))
      if(aic>=AIC)
      {
        aic=AIC
        knot1=k1[v1]
        det=deter
        orde1=m
      }
    }
    cat("orde = ", orde1, "knot 1 = ", knot1, "AIC = ", AIC, "determinan = ", deter,"\n")
    
    aic=10^10
    }
}
  return(c("orde"=orde1, "knot"=knot1,  "AIC"=AIC, "determinan"=deter))
}
listLambda= seq(from = 0.001, to = 0.99, length.out = 100)
listAIC=c()
for (lambdaVal in seq(from = 0.001, to = 0.99, length.out = 100)) {
  aicResult = aic1(x, y, m, k, lambdaVal) # TODO: MASUKAN x, y, m, k
  listAIC = c(listAIC, aicResult["AIC"])
}

dfResult=data.frame(Lambda=listLambda, AIC=listAIC)
View(dfResult)

library(ggplot2)
ggplotAIC<-ggplot(dfResult, aes(listLambda, listAIC)) +
  geom_line(color = "red") +
  geom_point(color = "red", size = 1) +
  xlab("Lambda") + ylab("AIC")+ ggtitle("Hubungan Antara Lambda dan Nilai AIC")
ggplotAIC+theme(plot.title=element_text(hjust = 0.5))


aic2<-function(x,y,m,k,lambda)
{
  n=length(y) #panjang data
  aic=10^10 #10 = mendefinisikan nilainya
  t1=seq(min(x),max(x),length.out=100) #seq = membuat vektor
  p=rep(0,(n-2))  #perulangan dari 0 sampai n-2
  for(z in 1:(n-2))
  {
    p[z]=t1[z+1] 
  }
  comb1=combn(p,2,FUN = NULL) #combinasi
  c1=t(comb1) #tranpose
  #for(m in 3)
  for (m in 2:4) #orde 2-4
  {
    k1=c1[,1] #kolom = kanan 
    k2=c1[,2] 
    K1=length(k1) 
    for (v1 in 1:K1) { #letak titik knot 1 - K titik knot
      bs1=bs(x, df=NULL, knot=c(k1[v1],k2[v1]), degre=m-1, intercept=TRUE, Boundary.knots=range(x))
      B = cbind(bs1)
      D = diag(ncol(bs1))
      BtB = t(B) %*% B
      deter=det(BtB) 
      z = MPL(BtB)
      beta= MPL (BtB + lambda * t(D) %*% D) %*% (t(B) %*% y)
      S= B %*% MPL (BtB + lambda * t(D) %*% D) %*% t(B)
      yhat= B %*% beta
      MSE = t(y- yhat) %*% (y-yhat)/n
      I <- matrix(0, ncol = n, nrow = n)
      for (v in 1:n) 
        I[v,v]=1
      k = sum(diag(I-S))
      AIC = (n+((n*log(2*pi))+n*(log(MSE))+(2*(2+m))))
      if(aic>=AIC)
      {
        aic=AIC
        knot1=k1[v1]
        knot2=k2[v1]
        orde1=m
      }
    }
    cat("orde = ", orde1, "knot 1 = ", knot1, "knot 2 =",knot2,  "AIC = ", AIC, "determinan = ", deter,"\n")
    
    aic=10^10
  }
}
  return(c("orde"=orde1, "knot 1"=knot1, "knot 2"=knot2, "AIC"= AIC, "determinan"=deter))
} 
listLambda= seq(from = 0.001, to = 0.99, length.out = 100)
listAIC = c()
for (lambdaVal in seq(from = 0.001, to = 0.99, length.out = 100)) {
  aicResult = aic2(x, y, m, k, lambdaVal) # TODO: MASUKAN x, y, m, k
  listAIC = c(listAIC, aicResult["AIC"])
}

dfResult=data.frame(Lambda=listLambda, AIC=listAIC)
View(dfResult)

library(ggplot2)
ggplotAIC<-ggplot(dfResult, aes(listLambda, listAIC)) +
  geom_point(color = "red", size = 1) +
  xlab("Lambda") + ylab("AIC")+ ggtitle("Hubungan Antara Lambda dan Nilai AIC")
ggplotAIC+theme(plot.title=element_text(hjust = 0.5))


aic3<-function(x,y,m,k,lambda)
{
  n=length(y) #panjang data
  print(n)
  aic=10^10 #10 = mendefinisikan nilainya
  t1=seq(min(x),max(x),length.out=100) #seq = membuat vektor
  p=rep(0,(n-2))  #perulangan dari 0 sampai n-2
  for(z in 1:(n-2))
  {
    p[z]=t1[z+1] 
  }
  comb1=combn(p,3,FUN = NULL) #combinasi
  c1=t(comb1) #tranpose
  #for (m in 2)
  for (m in 2:4) #orde 2-4
  {
    k1=c1[,1] #kolom = kanan 
    k2=c1[,2] 
    k3=c1[,3]
    K1=length(k1) 
    for (v1 in 1:K1) { #letak titik knot 1 - K titik knot
      bs1=bs(x, df=NULL, knot=c(k1[v1],k2[v1],k3[v1]), degre = m-1, intercept=TRUE, Boundary.knots=range(x))
      B = cbind(bs1)
      D = diag(ncol(bs1))
      BtB = t(B) %*% B
      deter=det(BtB) 
      z = MPL(BtB)
      beta= MPL (BtB + lambda * t(D) %*% D) %*% (t(B) %*% y)
      S= B %*% MPL (BtB + lambda * t(D) %*% D) %*% t(B)
      yhat= B %*% beta
      MSE = t(y- yhat) %*% (y-yhat)/n
      I <- matrix(0, ncol = n, nrow = n)
      for (v in 1:n) 
        I[v,v]=1
      k = sum(diag(I-S))
      AIC = (n+((n*log(2*pi))+n*(log(MSE))+(2*(3+m))))
      if(aic>=AIC)
      {
        aic=AIC
        knot1=k1[v1]
        knot2=k2[v1]
        knot3=k3[v1]
        orde1=m
      }
    }
    cat("orde = ", orde1, "knot 1 = ", knot1, "knot 2 =",knot2, "knot 3 = ", knot3 ,  "AIC = ", AIC, "determinan = ", deter,"\n")
    
    aic=10^10
  }
}
  return(c("orde"=orde1, "knot"=knot1, "knot 2 =",knot2, "knot 3 = ", knot3 , "AIC"=AIC, "determinan"=deter))
}
listLambda= seq(from = 0.001, to = 0.99, length.out = 100)
listAIC = c()
for (lambdaVal in seq(from = 0.001, to = 0.99, length.out = 100)) {
  aicResult = aic3(x, y, m, k, lambdaVal) # TODO: MASUKAN x, y, m, k
  listAIC = c(listAIC, aicResult["AIC"])
}

dfResult=data.frame(Lambda=listLambda, AIC=listAIC)
View(dfResult)

library(ggplot2)
ggplotAIC<-ggplot(dfResult, aes(listLambda, listAIC)) +
  geom_point(color = "red", size = 1) +
  xlab("Lambda") + ylab("AIC")+ ggtitle("Hubungan Antara Lambda dan Nilai AIC")
ggplotAIC+theme(plot.title=element_text(hjust = 0.5))



aic4<-function(x,y,m,k,lambda)
{
  n=length(y) #panjang data
  print(n)
  aic=10^10 #10 = mendefinisikan nilainya
  t1=seq(min(x),max(x),length.out=100) #seq = membuat vektor
  p=rep(0,(n-2))  #perulangan dari 0 sampai n-2
  for(z in 1:(n-2))
  {
    p[z]=t1[z+1] 
  }
  comb1=combn(p,4,FUN = NULL) #combinasi
  c1=t(comb1) #tranpose
  #for (m in 2)
  for (m in 2:4) #orde 2-4
  {
    k1=c1[,1] #kolom = kanan 
    k2=c1[,2] 
    k3=c1[,3]
    k4=c1[,4]
    K1=length(k1) 
    for (v1 in 1:K1) { #letak titik knot 1 - K titik knot
      bs1=bs(x, df=NULL, knot=c(k1[v1],k2[v1],k3[v1],k4[v1]), degre=m-1, intercept=TRUE, Boundary.knots=range(x))
      B = cbind(bs1)
      D = diag(ncol(bs1))
      BtB = t(B) %*% B
      deter=det(BtB) 
      z = MPL(BtB)
      beta= MPL (BtB + lambda * t(D) %*% D) %*% (t(B) %*% y)
      S= B %*% MPL (BtB + lambda * t(D) %*% D) %*% t(B)
      yhat= B %*% beta
      MSE = t(y- yhat) %*% (y-yhat)/n
      I <- matrix(0, ncol = n, nrow = n)
      for (v in 1:n) 
        I[v,v]=1
      k = sum(diag(I-S))
      AIC = (n+((n*log(2*pi))+n*(log(MSE))+(2*(4+m))))
      if(aic>=AIC)
      {
        aic=AIC
        knot1=k1[v1]
        knot2=k2[v1]
        knot3=k3[v1]
        knot4=k4[v1]
        orde1=m
      }
    }
    cat("orde = ", orde1, "knot 1 = ", knot1, "knot 2 =",knot2, "knot 3 =", knot3, "knot 4 =", knot4 ,  "AIC = ", AIC, "determinan = ", deter,"\n")
    
    aic=10^10
  }
  return(c("orde"=orde1, "knot"=knot1, "knot 2 =",knot2, "knot 3 = ", knot3, "knot 4 =", knot4  , "AIC"=AIC, "determinan"=deter))
}
listLambda= seq(from = 0.001, to = 0.99, length.out = 100)
listAIC = c()
for (lambdaVal in seq(from = 0.001, to = 0.99, length.out = 100)) {
  aicResult = aic4(x, y, m, k, lambdaVal) # TODO: MASUKAN x, y, m, k
  listAIC = c(listAIC, aicResult["AIC"])
}

dfResult=data.frame(Lambda=listLambda, AIC=listAIC)
View(dfResult)

library(ggplot2)
ggplotAIC<-ggplot(dfResult, aes(listLambda, listAIC)) +
  geom_point(color = "red", size = 1) +
  xlab("Lambda") + ylab("AIC")+ ggtitle("Hubungan Antara Lambda dan Nilai AIC")
ggplotAIC+theme(plot.title=element_text(hjust = 0.5))







#estimasi parameter
pspline<-function(x,y,m,k,lambda){
  n<-length(y)
  print(n)
  knot<-c(k)
  knot<-as.matrix(knot)
  k1<-length(knot)
  bs1 = bs(x,df=NULL, knots=k, degree=m-1, intercept=TRUE, Boundary.knots=range(x))
  B = cbind(bs1)
  D = diag(ncol(bs1))
  BtB = t(B) %*% B
  z = MPL(BtB)
  beta= MPL (BtB + lambda * t(D) %*% D) %*% (t(B) %*% y)
  S= B %*% MPL (BtB + lambda * t(D) %*% D) %*% t(B)
  cat("Nilai parameter adalah", "\n", beta, "\n")
  yhat = (B %*% beta)
  MSE = ((t(y-yhat)) %*% (y-yhat))/n
  I <- matrix(0, ncol=n, nrow = n)
  for (v in 1:n) 
    I[v,v] = 1
  l = sum(diag(I-S))
  AIC = (n+((n*log(2*pi))+n*(log(MSE))+(2*(3+m))))
  cat("NIlai AIC adalah ","\n",AIC,"\n")
}

pspline(x,y,2, k=c(69.39747 , 69.81707), 0.001)

#pengujian parameter
Uji_Parameter<-function(x,y,m,k,lambda){
  n<-length(y)
  knot<-c(k)
  knot<-as.matrix(knot)
  k1<-length(knot)
  orde = (m-1)
  bs1=bs(x,df=NULL, knots=k, degree=m-1, intercept=TRUE, Boundary.knots=range(x))
  B = cbind(bs1)
  D = diag(ncol(bs1))
  BtB = t(B) %*% B
  z = MPL(BtB)
  beta= MPL (BtB + lambda * t(D) %*% D) %*% (t(B) %*% y)
  Beta = as.matrix(beta)
  cat("Nilai parameter adalah", "\n", "\n")
  print(Beta)
  S= B %*% MPL (BtB + lambda * t(D) %*% D) %*% t(B)
  yhat = (B %*% beta)
  ybar <- mean(y)
  MSE = ((t(y-yhat)) %*% (y-yhat))/n
  MSE <- as.numeric(MSE)
  SSR <- sum((yhat-ybar)^2)
  MSR <- SSR/((orde+k1)-1)
  SSE <- sum((y-yhat)^2)
  MSER <- SSE/(n-(orde+k1))
  JKT = sum((y-ybar)^2)
  Fhit_1 = MSR/MSER
  cat("Nilai F hitung pengujian serentak adalah","\n", Fhit_1,"\n")
  cat("Kesimpulan hasil uli serentak","\n")
  cat("-----------------------------------","\n")
  cat("Analysis Of Variance (ANOVA)","\n")
  cat("===================================================","\n")
  cat("Sumber df   SS    MS      Fhitung","\n")
  cat("Regresi",((orde+k1)-1), "",SSR,"",MSR,"",Fhit_1,"\n")
  cat("Error", (n-(orde+k1)), "", SSE, "", MSER, "", "\n")
  cat("Total ", (n-1), "", JKT, "", "\n")
  print(MSE)
  cat("===================================================", "\n")
}

Uji_Parameter(x,y,2, k=c(69.39747 , 69.81707), 0.001)

#mencari nilai ftabel
qf(p = 0.05, df1=2, df2=61, lower.tail = FALSE)

#Menghitung nilai prediksi
prediksi<-function(x,y,m,k,lambda){
  knot<-c(k)
  knot<-as.matrix(knot)
  k<-length(knot)
  bs1=bs(x,df=NULL, knots=k, degree=m-1, intercept=TRUE, Boundary.knots=range(x))
  B = cbind(bs1)
  D = diag(ncol(bs1))
  BtB = t(B) %*% B
  z = MPL(BtB)
  beta= MPL(BtB + lambda*t(D) %*% D) %*% (t(B) %*% y)
  Beta = as.matrix(beta)
  yhat = (B %*% beta)
  prediksi = cbind(y,yhat)
  prediksi_B = as.matrix(prediksi)
  cat("Nilai prediksi adalah \n \n")
  print(prediksi_B)
}
prediksi(x,y,2, k=c(69.39747 , 69.81707), 0.001)
yhat = prediksi(x,y,2, k=c(69.39747 , 69.81707), 0.001)
MAPE = mean(abs(y-yhat)/y)*100
MAPE
