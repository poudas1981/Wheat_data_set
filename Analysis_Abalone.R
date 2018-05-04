
.libPaths(new = "/work/statsgeneral/vcdim/Code/packages")
.libPaths() #Check to see it is #1 in the search path
install.packages('ncvreg', repos="http://cran.r-project.org")
install.packages('doParallel', repos="http://cran.r-project.org")
install.packages('doRNG', repos="http://cran.r-project.org")


library(ncvreg)
library(MASS)
library(doParallel)
#library(doRNG)



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


Func_Stat = function(df, data, VCD, m = 10){
  data = data[,1:df]
  Model1 = lm(SPEED ~., data = data)
  Risk2 = sum(Model1$residuals^2)
  BIC = round(BIC(Model1))
  AIC = round(AIC(Model1))
  ERM1 = round(Gaby(Loss = Risk2,eta = 0.05,
                    n = nrow(data), h = VCD[df-1], m))
  ERM2 = round(ERM(Loss = Risk2,eta = 0.05,
                   n = nrow(data), h = VCD[df-1],  m))
  output = c(VCD[df-1],ERM1, ERM2, AIC, BIC)
  return(output)
  
}

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
    y = data1[test_i, ][,1]
    
    
    fitted_models <- apply(t(df), 2, function(degf) lm(SPEED ~ ., data = train_xy[,1:degf]))
    
    pred <- mapply(function(obj, degf) predict(obj, test_x[, 1:degf]), 
                   fitted_models, df)
    cv_tmp[k, ] <- sapply(as.list(data.frame(pred)), function(y_hat) mean((y - y_hat)^2, na.rm = TRUE))
  }
  return(cv_tmp)
}

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
        set.seed(i*bprim*b+1)
        Index1 = sample(nrow(data),size = 2*n,replace = TRUE)
        Mydata = data[Index1,]
        # Step two: split the data into two groups
        Index = sample(nrow(Mydata),size = n,replace = FALSE)
        SampleData = Mydata
        G1 = SampleData[Index,]
        G2 = SampleData[-Index,]
        
        # Now lets fit a model using the modify dataset.
        Model1 = lm(Rings ~ .,  x = TRUE, y = TRUE, data = G1)
        
        Model2 = lm(Rings ~ .,  x = TRUE, y = TRUE, data = G2)
        
        
        FirstHalfN = matrix(NA, nrow =length(Index), ncol =m)
        FirstHalfQ = matrix(NA, nrow =length(Index), ncol =m)
        colnames(FirstHalfN) =paste("N",1:m,sep="")
        colnames(FirstHalfQ) = paste("Q",1:m,sep="")
        
        
        X1 = as.matrix(Model1$x)
        Y1 = data.frame(Model1$y)
        X2 = as.matrix(Model2$x)
        Y2 = data.frame(Model2$y)
        Para1 = Model1$coefficients
        Para1 = matrix(Para1, ncol = 1)
        
        Para2 = matrix(Model2$coefficients, ncol = 1)
        #Para2 = (Para2[2:length(Para2)])
        
        #Pred1 = PredSqrt[Index]
        qstarj = matrix(NA, nrow = 1, ncol = m)
        Pred1 = (X2%*%Para1 - Y2)^2
        
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
        Pred2 = (X1%*%Para2 - Y1)^2
        
        #Pred2 = (predict(Model2, X1)-Y1)^2
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


B = 50
m = 10
NL = seq(from = 1000, to = 4000, by = 500)

data = read.csv(file = "/work/statsgeneral/vcdim/Code/MerlinAbalone.csv", header = T)

#data = read.csv("C:/Users/merli/OneDrive/Documents/Code/MerlinAbalone.csv", header = TRUE)
SexData = data[,-1]
datascale = as.data.frame(scale(SexData, center = TRUE, scale = TRUE))
corr = cor(datascale)
Name = c("Length",  "Diameter", "Height", "Wholeweight", "Shuckedweight", "Visceraweight", "Shellweight")
cor_name = corr[,"Rings"][-8]
name_order = Name[order(cor_name, decreasing = TRUE)]
datascale = datascale[, name_order]
set.seed(10)






# estimate Chixi
est <- sapply(2:ncol(datascale), function(x, data, NL, B, m) Chxi(data[,1:x], NL, B, m), data = datascale, NL = NL, B = B, m = m)


VCD = sapply(1:ncol(est), function(x, data, NL, B, m) vcd_funct(est[,x], data, NL, B, m), data = datascale, NL = NL, B=B, m = m)
cat("The estimated VCD are:", "\n")
VCD


cv_tmp = Model_cv(Tour, n_folds <- 10)

cv = colMeans(cv_tmp)

st = t(sapply(df, function(x, data, VCD, m) Func_Stat(x, data, VCD, m), data = datascale,VCD = VCD, m = 10))
colnames(st) = c("VCD","ERM1", "ERM2", "AIC", "BIC")
st
data.frame(st, cv = cv)