.libPaths(new = "/work/statsgeneral/vcdim/Code/packages")
.libPaths() #Check to see it is #1 in the search path
install.packages(c('ncvreg', 'doParallel', 'polynom', 'parallel'), repos="http://cran.r-project.org")

#library(polynom)
library(MASS)
library(doParallel)
library(parallel)

#library(MASS)
#library(doParallel)
#library(parallel)

Model_cv = function(data, n_folds){
  set.seed(10)
  data1 = data.frame(scale(data, scale = TRUE, center = TRUE))
  df = 2:ncol(data1) 
  folds_i <- sample(rep(1:n_folds, length.out = dim(data)[1]))
  cv_tmp <- matrix(NA, nrow = n_folds, ncol = length(df))
  for (k in 1:n_folds) {
    test_i <- which(folds_i == k)
    train_xy <- data1[-test_i, ]
    test_x <- data1[test_i, ]
    y = data1[test_i, ][,"YIELD"]
    
    
    fitted_models <- apply(t(df), 2, function(degf) lm(YIELD ~ ., data = train_xy[,1:degf]))
    
    pred <- mapply(function(obj, degf) predict(obj, test_x[, 1:degf]), 
                   fitted_models, df)
    cv_tmp[k, ] <- sapply(as.list(data.frame(pred)), function(y_hat) mean((y - y_hat)^2, na.rm = TRUE))
  }
  return(cv_tmp)
}

ERM = function(Loss, eta, n, h, m){
  coef1 = (m^2)/(2*n)*log((2*m/eta)*((2*n*exp(1)/h)^h))
  coef2 = (1 + sqrt(1 + (4*n*Loss)/((m^2)*log((2*m/eta)*((2*n*exp(1)/h)^h)))))
  Bound = Loss + coef1*coef2
  return(Bound)
}
Gaby = function(Loss, eta, n, h, m){
  coef = (m)*sqrt((1/n)*log((2*m/eta)*((2*n*exp(1)/h)^h)))
  Bound = Loss + coef
  return(Bound)
}
phiTheo5 = function(n,x, c1, c2){
  c2 = 0
  c1*sqrt((x/n)*log(2*n*exp(1)/x)) + c2*(x/n)*log(2*n*exp(1)/x)
  #0.5*((m^2)/n)*log(2*n*exp(1)/(x))*(1 + sqrt(1 + (x/n)*log((2*n*exp(1))/x)))
}

phiTheo = function(x){
  0.2*sqrt((x/250)*log(2*250*exp(1)/x)) #+ 0.2*(x/250)*log(2*250*exp(1)/x)
  #0.5*((m^2)/n)*log(2*n*exp(1)/(x))*(1 + sqrt(1 + (x/n)*log((2*n*exp(1))/x)))
}

phiTheo51 = function(x){
  0.33*sqrt((10/x)*log(2*x*exp(1)/10)) + 0.01*(10/x)*log(2*x*exp(1)/10)
  #0.5*((m^2)/n)*log(2*n*exp(1)/(x))*(1 + sqrt(1 + (x/n)*log((2*n*exp(1))/x)))
}

C1C2 = function(MatChixi, h, NL, c1, c2){
  x1 = c1*sqrt((h/NL)*log(2*NL*exp(1)/h)) 
  x2 = c2*(h/NL)*log(2*NL*exp(1)/h)
  x2 = 0
  out = (1/length(NL))*sum((MatChixi - x1 - x2)^2)
}

C1C2ratio = function(MatChixi, h, NL, c1, c2){
  x1 = c1*sqrt((h/NL)*log(2*NL*exp(1)/h)) 
  x2 = c2*(h/NL)*log(2*NL*exp(1)/h)
  x = x1 + x2
  out = (1/length(NL))*sum((MatChixi/x - 1)^2)
}



vcfunctratio = function(MatChixi,NL,x,m,c1,c2){
  row = 1
  Sum = 0
  while(row <= length(MatChixi)){
    Sum = Sum + (MatChixi[row]/phiTheo5(n=NL[row],x,c1,c2) - 1 )^2
    row = row + 1
  }
  Sum = (1/length(NL))*Sum
  return(Sum)
}


vcfunct = function(MatChixi,NL,x,m,c1,c2){
  row = 1
  Sum = 0
  while(row <= length(MatChixi)){
    Sum = Sum + (MatChixi[row] - phiTheo5(n=NL[row],x,c1,c2))^2
    row = row + 1
  }
  Sum = (1/length(NL))*Sum
  return(Sum)
}

Vapbound = function(x,NL){
  0.16*((log(2*(NL/x))+1)/(NL/x-0.15))*(1+ sqrt(1+ 1.2*(NL/x-0.15)/(log(2*NL/x)+1)))
}

vapfunct = function(MatChixi,NL,x){
  row = 1
  Sum = 0
  while(row <= length(MatChixi)){
    Sum = Sum + (MatChixi[row] - Vapbound(x,NL[row]))^2
    row = row + 1
  }
  Sum = (1/length(MatChixi))*Sum
}

#data = read.csv("C:/Users/merli/OneDrive/Documents/DataSet/SNPWheatData.csv", header = T)
data = read.csv("C:/Users/merli/OneDrive/Documents/DataSet/FullWheatData.csv", header = T)

str(data)
names(data)
#data = read.csv(file = "/work/statsgeneral/vcdim/Code/FullWheatData.csv", header = T)
head(data)

B=50

#data = read.csv(file = "/work/statsgeneral/vcdim/Code/WheatData.csv", header = T)
#data = as.data.frame(data, center = TRUE, scale = TRUE)
NL = c(450, 500, 550, 600, 650, 700, 750)
Loc = levels(data$LOCATION)
k = 2
data237 = subset(data, data$LOCATION == Loc[k])
var = c('YIELD','HT', 'TSTWT', 'TKWT', 'SPSM', 'KPS', 'KPSM', 
        'barc67', 'cmwg680bcd366', 'bcd141', 'barc86', 'gwm155',
        'barc12','IBLK')
data2 = data237[,var]


BigModel = lm(YIELD ~ TKWT + TSTWT + SPSM + KPS + KPSM + HT + 
                I(TKWT^2)+ I(TKWT*TSTWT) + I(TKWT*SPSM) + I(TKWT*KPS) + I(TKWT*KPSM) + I(TKWT*SPSM) + I(TKWT*HT) +
                I(TSTWT^2) + I(TSTWT*SPSM) + I(TSTWT*KPS) + I(TSTWT*KPSM) + I(TSTWT*HT) +
                I(SPSM^2) + I(SPSM*KPS) + I(SPSM*KPSM) + I(SPSM*HT) + 
                I(KPS^2) + I(KPS*KPSM) + I(KPS*HT) + 
                I(KPSM^2) + I(KPSM*HT) +
                I(HT^2) + barc67 + cmwg680bcd366 + bcd141 + barc86 + gwm155 +
                barc12, data = data2, x=TRUE, y=TRUE)
Xdat = BigModel$x[,-1]
Ydat = BigModel$y
cor(Ydat, Xdat)
ddd = as.matrix(cbind(Ydat,Xdat))
cor(ddd)[1,]
Name = c('TKWT', 'TSTWT', 'SPSM', 'KPS', 'KPSM', 'HT', 
         'I(TKWT^2)', 'I(TKWT*TSTWT)', 'I(TKWT*SPSM)', 'I(TKWT*KPS)', 'I(TKWT*KPSM)', 'I(TKWT*HT)', 
         'I(TSTWT^2)', 'I(TSTWT*SPSM)',  'I(TSTWT*KPS)', 'I(TSTWT*KPSM)', 'I(TSTWT*HT)',
         'I(SPSM^2)', 'I(SPSM*KPS)', 'I(SPSM*KPSM)', 'I(SPSM*HT)',
         'I(KPS^2)', 'I(KPS*KPSM)', 'I(KPS*HT)', 
         'I(KPSM^2)', 'I(KPSM*HT)', 
         'I(HT^2)',
         'barc67', 'cmwg680bcd366', 'bcd141', 'barc86', 'gwm155', 'barc12')
Cor = as.matrix(round(abs(cor(Ydat, Xdat)),4))
Name[order(cor(ddd)[1,][-1], decreasing = TRUE)]
######################################################################################################
# Order of inclusion of covariates using SNP data in Licoln 01
######################################################################################################


AA = Name[order(cor(ddd)[1,][-1], decreasing = TRUE)]
Big_order = lm(YIELD ~ I(TKWT*KPSM) + (TSTWT*KPS)  + KPSM + I(SPSM*KPS) + I(KPSM^2) + I(SPSM*KPSM) +    
                 I(TKWT*SPSM) + I(TSTWT*SPSM) + SPSM + I(KPSM*HT) + I(SPSM^2) + I(KPS*KPSM) +        
                 I(SPSM*HT) + TSTWT + I(TSTWT^2) +         
                 barc67 + I(TKWT*TSTWT)+ barc86 +              
                 TKWT+ I(TKWT^2) + cmwg680bcd366 +       
                 bcd141 + I(TKWT*KPS) + gwm155 +     
                 I(TSTWT*KPS)  + barc12 + I(KPS^2) +          
                 KPS + I(TKWT*HT) + I(KPS*HT) +           
                 I(TSTWT*HT) + HT + I(HT^2),
               data = data2, x = TRUE, y = TRUE)

X_big = Big_order$x[,-1]
Y_big = Big_order$y

X_data = data.frame(YIELD = Y_big, X_big)
X_data = data.frame(scale(X_data, center = TRUE, scale = TRUE))
#######################################################################################################

Chxi = function(data, NL, B, m){
  MatChixi = matrix(NA, nrow = B, ncol = length(NL))
  
  Qstar = function(j,B,m){
    (2*j+1)*B/(2*m)
  }
  Lower = function(j,B,m){
    j*B/m
  }
  Upper = function(j,B,m){
    (j+1)*B/m
  }
  
  for(i in 1:length(NL)){
    # step one: we need to generate 2n data points
    n = NL[i]
    bprim = 1
    while(bprim<(B+1)){
      b = 1
      #sumdiff = 0
      Matxi = matrix(NA, nrow = B, ncol = m)
      while(b < (B+1)){
        cat('Bootstrap #', bprim, 'second boot', b, '\n')
        set.seed(i*bprim*b+1)
        Index1 = sample(nrow(data),size = 2*n,replace = TRUE)
        Mydata = data[Index1,]
        # Step two: split the data into two groups
        Index = sample(nrow(Mydata),size = n,replace = FALSE)
        SampleData = Mydata
        G1 = SampleData[Index,]
        G2 = SampleData[-Index,]
        
        # Now lets fit a model using the modify dataset.
        Model1 = lm(YIELD ~ .,  x = TRUE, y = TRUE, data = G1)
        
        Model2 = lm(YIELD ~ .,  x = TRUE, y = TRUE, data = G2)
        
        FirstHalfN = matrix(NA, nrow =length(Index), ncol =m)
        FirstHalfQ = matrix(NA, nrow =length(Index), ncol =m)
        colnames(FirstHalfN) =paste("N",1:m,sep="")
        colnames(FirstHalfQ) = paste("Q",1:m,sep="")
        
        X1 = data.frame(Model1$x)
        Y1 = data.frame(Model1$y)
        X2 = data.frame(Model2$x)
        Y2 = data.frame(Model2$y)
        
        #Pred1 = PredSqrt[Index]
        qstarj = matrix(NA, nrow = 1, ncol = m)
        Pred1 = (predict(Model1, X2)-Y2)^2
        for(ind in 0:(m-1)){
          qstarj[1,(1+ind)]= Qstar(ind,max(Pred1),m)
          for(k in 1:length(Pred1)){
            if(Pred1[k]>=Lower(ind,max(Pred1),m) & Pred1[k]<Upper(ind,max(Pred1),m)){
              FirstHalfN[k,(ind+1)] = 1
              FirstHalfQ[k,(ind+1)] = Qstar(ind,max(Pred1),m)
            }
          }
        }
        N1 = apply(FirstHalfN, 2, sum, na.rm = TRUE)
        Q1 = apply(FirstHalfQ, 2, sum, na.rm = TRUE)
        nustar1 = N1*qstarj/(n)
        #nustar1 = N1*Q1/(n)
        
        
        SecondHalfN = matrix(NA, nrow =length(Index), ncol =m)
        SecondHalfQ = matrix(NA, nrow =length(Index), ncol =m)
        
        colnames(SecondHalfN) =paste("N",1:m,sep="")
        colnames(SecondHalfQ) =paste("Q",1:m,sep="")
        
        qstarj = matrix(NA, nrow = 1, ncol = m)
        
        Pred2 = (predict(Model2 ,X1)-Y1)^2
        for(indx in 0:(m-1)){
          qstarj[1,(1+indx)]= Qstar(indx,max(Pred1),m)
          for(kk in 1:length(Pred2)){
            if(Pred2[kk]>=Lower(indx,max(Pred2),m) & Pred2[kk]<Upper(indx,max(Pred2),m)){
              SecondHalfN[k,(indx+1)] = 1
              SecondHalfQ[k,(indx+1)] = Qstar(indx,max(Pred2),m)
            }
          }
        }
        
        N2 = apply(SecondHalfN, 2, sum, na.rm = TRUE)
        Q2 = apply(SecondHalfQ, 2, sum, na.rm = TRUE)
        nustar2 = N2*qstarj/(n)
        #nustar2 = N2*Q2/(n)
        
        diff = abs(nustar2 - nustar1)
        #sumdiff = sumdiff + diff
        Matxi[b,] = diff
        
        b = b + 1
      }
      #MatChixi[bprim,i] = sum (apply(Matxi, 2, max, na.rm = T))
      MatChixi[bprim,i] = sum (apply(Matxi, 2, mean, na.rm = T))
      
      
      bprim = bprim + 1
    }
    
  }
  
  MeanChi = apply(MatChixi, MARGIN =2, FUN = mean)
}


Output = matrix(NA, ncol = 6, nrow = length(AA))
colnames(Output) = c('h', 'ERM1', 'ERM2', 'AIC', 'BIC', 'C')

for(l in 2:ncol(X_data)){
  dat = X_data[, 1:l]
  X = data.frame(cbind(scale(dat, center = TRUE, scale = TRUE)))
  #est = Chxi(data = X, NL = NL, B=5, m=10)
  Model1 = lm(YIELD ~ ., data = X, x = TRUE, y = TRUE)
  
  # estimate Chixi
  est = Chxi(data = X, NL = NL, B=5, m=10)
  
  #####################################################################
  # Estimate c1 and c2
  #####################################################################
  c1 =  seq(from = 0.01, to = 10, by = 0.01)
  c2 = 0
  Mertt = numeric()
  # for (t in 1:length(hh)) {
  ourestMat = matrix(NA, nrow = length(c1), ncol = length(c2))
  ourestMat2 = matrix(NA, nrow = length(c1), ncol = length(c2))
  for (j in 1:length(c1)) {
    for (g in 1:length(c2)) {
      ourestMat[j,g] = C1C2(est,ncol(data)-1, NL, c1 = c1[j], c2 = c2[g])
      ourestMat2[j,g] = C1C2ratio(est,ncol(data)-1, NL, c1 = c1[j], c2 = c2[g])
    }
  }
  Indexx  = which(ourestMat == min(ourestMat), arr.ind = TRUE)
  Indexx2 = which(ourestMat2 == min(ourestMat2), arr.ind = TRUE)
  c111 = c1[Indexx[1,1]]
  c222 = c2[Indexx[1,2]]
  
  c11r = c1[Indexx2[1,1]]
  c22r = c2[Indexx2[1,2]]
  cat("C1 is ", c111, "c2 is ", c222, 'the number of col in data is  is', ncol(data)-1, "\n")
  # estimate vc dimension using grid search
  range2 = seq(from = 1, to = 100, by = 1)
  MerlinMat1 = numeric(length(range2))
  MerlinMat1r = numeric(length(range2))
  
  for (kk in 1:length(range2)) {
    MerlinMat1[kk] = vcfunct(est,NL=NL,x=range2[kk], m=10, c1=c111, c2=c222)
    MerlinMat1r[kk] = vcfunctratio(est,NL=NL,x=range2[kk], m=10, c1=c11r, c2=c22r)
  }
  cat('The estimate vcdim is: ', range2[which.min(MerlinMat1)], " for Loc ", Loc[k], "\n")
  cat('The estimate vcdim is: ', range2[which.min(MerlinMat1r)], " for Loc ", Loc[k], "\n")
  Risk2 = sum(Model1$residuals^2)
  BIC = BIC(Model1)
  AIC = AIC(Model1)
  cat("The BIC is:", round(BIC), '\n')
  ERM1 = Gaby(Loss = Risk2,eta = 0.05,
              n = nrow(data), h = range2[which.min(MerlinMat1)], m=10)
  round(ERM1)
  ERM2 = ERM(Loss = Risk2,eta = 0.05,
             n = nrow(data), h = range2[which.min(MerlinMat1)], m=10)
  round(ERM2)
  dataa = data.frame(h = range2[which.min(MerlinMat1)], ERM1 = round(ERM1), ERM2 = round(ERM2), BIC = round(BIC), AIC = round(AIC))
  Output[l-1,] = c(range2[which.min(MerlinMat1)], round(ERM1), round(ERM2),  round(AIC), round(BIC), c111
  )
  #dataa
  
}

Output

cv_tmp = Model_cv(dat, n_folds <- 10)

cv = colMeans(cv_tmp)

data.frame(Output[,-6], cv)


X_big = Big_order$x[,-1]
Y_big = Big_order$y
X1Scaled = scale(X_big, center = TRUE, scale = TRUE)
Anal_data = data.frame(X1Scaled, IBLK = data2$IBLK)
#Anal_data = data.frame(X1Scaled)
library(ncvreg)
Model_SCAD = ncvreg(Anal_data, Y_big, penalty = 'SCAD')
plot(Model_SCAD)
cvfit2 <- cv.ncvreg(Anal_data, Y_big, penalty = "SCAD")
round(coef(cvfit2),2)
plot(cvfit2)
cvfit2$lambda.min
Param = t(Model_SCAD$beta)
source("C:/Users/merli/OneDrive/Documents/Code/LSA.r")
source("C:/Users/merli/OneDrive/Documents/Code/Lasso.r")
Ada_lasso = lasso.adapt.bic2(Anal_data, Y_big)
round(Ada_lasso$coeff,2)
round(Ada_lasso$intercept,2)
ada = lsa(Model1)

#######################################################################################
# model dev, impl , valida, risk model, forcasting credit handle 
# risk model enterprisewise data


