
.libPaths(new = "/work/statsgeneral/vcdim/Code/packages")
.libPaths() #Check to see it is #1 in the search path

library(polynom)
library(MASS)
library(nloptr)
library(MASS)
require(Hmisc)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

datasim = function(beta, mu, sigma, size,seed){
  Mean = rnorm(nrow(beta), mean = 5, sd = 2)
  Sigma = diag(x = 1, nrow=nrow(beta), ncol=nrow(beta))
  set.seed(seed)
  e.sim<-matrix(data = rnorm(n = size, mean = mu, sd = sigma), nrow = size, ncol = 1) 
  #x.sim = matrix(data = runif(n = N, min = -5, max = 5), nrow = nrow, ncol = ncol)
  x.sim = mvrnorm(n = size, Mean, Sigma, empirical = TRUE)
  y.sim = x.sim%*%beta + e.sim
  
  output = data.frame(y=y.sim, x.sim)
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

B=50

data = read.csv(file = "/work/statsgeneral/vcdim/Code/Tour1.csv", header = T)
#file_tour1 = 'C:/Users/merli/OneDrive/Documents/Code/Tour1.csv'
file_tour = 'C:/Users/merli/OneDrive/Documents/Thesis/Code/MerlinTourDeFrance.csv'
Tour = read.csv(file_tour, header = T)
#Tour1 = read.csv(file_tour1, header = T)
#Tour = Tour[-c(13:20),]
head(Tour)
pairs(Tour)

plot(Tour$YEAR, Tour$SPEED, ylab = paste("Average Speed of the winner"),
     xlab = "Year", main = "Tour De France Data")
library(ggplot2)
P1 = ggplot(data = Tour) + 
  geom_point(mapping = aes(x = YEAR, y = SPEED)) +
  labs(x = "YEAR", y = "")

P1
P2 = ggplot(data=Tour, aes(AGE)) + 
  geom_histogram(breaks=seq(15, 40, by =5), 
                 col="red", 
                 aes(fill=..count..)) +
  scale_fill_gradient("Count", low = "green", high = "red")+
  #ggtitle("Histogram for Age") +
  labs(x="Age", y="") 

P3 = ggplot(data=Tour, aes(STAGESWIN)) + 
  geom_histogram(breaks=seq(1, 30, by =5), 
                 col="red", 
                 aes(fill=..count..)) +
  scale_fill_gradient("Count", low = "green", high = "red")+
  #ggtitle("Histogram for the number \nof of Stages won") +
  labs(x="Stages Won", y="") 

multiplot(P2, P3, cols=2)
plot(Tour$YEAR, Tour$SPEED)
identify(Tour$YEAR, Tour$SPEED)

plot(Tour$DIST, Tour$SPEED, ylab = paste("Average Speed of the winner"),
     xlab = "Distance Travel", main = "Tour De France Data")

identify(Tour$DIST, Tour$SPEED)
plot(Tour$SPEED~Tour$DIST, ylab = paste("Average Speed of the winner"),
     xlab = "Distance Travel", main = "Tour De France Data")

plot(Tour$STAGES~Tour$YEAR, ylab = paste("Average Speed of the winner"),
     xlab = "STAGES", main = "Year")

points(Tour$STAGESWIN~Tour$YEAR, ylab = paste("Average Speed of the winner"),
     xlab = "STAGES", main = "Year", col = "red")

data = data.frame(SPEED = Tour$SPEED, YEAR = Tour$YEAR, 
                  DIST=Tour$DIST, STAGESWIN = Tour$STAGES, 
                  AGE = Tour$AGE)
dataYD = data.frame(data, YD = data$YEAR*data$DIST)
cor(dataYD)
data1 = data.frame(SPEED = Tour$SPEED, YEAR = Tour$YEAR, 
                  DIST=Tour$DIST)
data$YEAR = seq(from = 1, to = nrow(data), by = 1)
data1$YEAR = seq(from = 1, to = nrow(data1), by = 1)
library(monreg)

df  = monreg(x=data$YEAR, y = data$SPEED, degree = 0, hd = .5, hr = .5)
data$SPEED = df$estimation
data = as.data.frame(scale(data, center = TRUE, scale = TRUE))
cor(data)
dataYD = data.frame(data, YD = data$YEAR*data$DIST)
cor(dataYD)
pairs(data)
data1 = as.data.frame(scale(data1, center = TRUE, scale = TRUE))
cor(data1)
cor(data, method = "kendall")
#Tour1 = as.data.frame(scale(data, center = TRUE, scale = TRUE))
#Tour1$YEAR = sign(Tour1$YEAR)*(abs(Tour1$YEAR)^(1.5))
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
        #G1 = as.data.frame(scale(SampleData[Index,], center = TRUE, scale = TRUE))
        #G2 = as.data.frame(scale(SampleData[-Index,], center = TRUE, scale = TRUE))
        
        # Now lets fit a model using the modify dataset.
        Model1 = lm(G1$SPEED ~  I(G1$YEAR^2)+ G1$YEAR  + G1$DIST + I(G1$DIST^2)  + I(G1$YEAR*G1$DIST), x = TRUE, y = TRUE, data = G1)
        Model2 = lm(G2$SPEED ~  I(G2$YEAR^2)+ G2$YEAR  + G2$DIST + I(G2$DIST^2)  + I(G2$YEAR*G2$DIST), x = TRUE, y = TRUE, data = G2)
        #PredSqrt = (Model$residuals)^2
        
        FirstHalfN = matrix(NA, nrow =length(Index), ncol =m)
        FirstHalfQ = matrix(NA, nrow =length(Index), ncol =m)
        colnames(FirstHalfN) =paste("N",1:m,sep="")
        colnames(FirstHalfQ) = paste("Q",1:m,sep="")
        
        X1 = data.frame(Model1$x)
        Y1 = Model1$y
        X2 = data.frame(Model2$x)
        Y2 = Model2$y
        #Pred1 = PredSqrt[Index]
        qstarj = matrix(NA, nrow = 1, ncol = m)
        Pred1 = (predict(Model1, X2)- Y2)^2
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
      MatChixi[bprim,i] = sum (apply(Matxi, 2, max, na.rm = T))
      #MatChixi[bprim,i] = sum (apply(Matxi, 2, mean, na.rm = T))
      
      
      bprim = bprim + 1
    }
    
  }
  
  MeanChi = apply(MatChixi, MARGIN =2, FUN = mean)
}
set.seed(13000000)


Size = 400
NL = seq(from = 20, to = 100, by = 10)
B = 50

# estimate Chixi
    est = Chxi(data = data, NL = NL, B=B, m=10)
    #####################################################################
    # Estimate c1 and c2
    #####################################################################
    c1 =  seq(from = 0.01, to = 100, by = 0.01)
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
    range2 = seq(from = 1, to = 30, by = 1)
    MerlinMat1 = numeric(length(range2))
    MerlinMat1r = numeric(length(range2))
	
    for (kk in 1:length(range2)) {
      MerlinMat1[kk] = vcfunct(est,NL=NL,x=range2[kk], m=10, c1=c111, c2=c222)
      MerlinMat1r[kk] = vcfunctratio(est,NL=NL,x=range2[kk], m=10, c1=c11r, c2=c22r)
    }
    cat('The estimate vcdim is: ', range2[which.min(MerlinMat1)], "\n")
    cat('The estimate vcdim is: ', range2[which.min(MerlinMat1r)], "\n")

    Name1 = paste('TourDeFrancePlot', 'pdf', sep='.')
    FileName1 = "/work/statsgeneral/vcdim/Code"
    pdf(paste(FileName1,Name1, sep='/'))
    plot(MerlinMat1,ylab =paste('Expected Max diff') , xlab = paste('The estimated VC Dim is: ',range2[which.min(MerlinMat1)]),
         main = paste('Analysis of Tour De France Data'))
    dev.off()
    library(monreg)
    df  = monreg(x=data$YEAR, y = data$SPEED, degree = 0, hd = .5, hr = .5)
    df2 = monreg(x=(data$YEAR)^1.5,  y = data$SPEED, degree = 0, hd = .5, hr = .5)
    df3 = monreg(x=(data$YEAR)^1.25, y = data$SPEED, degree = 0, hd = .5, hr = .5)
    df4 = monreg(x=log(data$YEAR),   y = data$SPEED, degree = 0, hd = .5, hr = .5)
    df5 = monreg(x=(data$YEAR)^1.1,   y = data$SPEED, degree = 0, hd = .5, hr = .5)
    par(mfrow = c(2,1))
      plot(data$YEAR, data$SPEED, lty = 1,   col = 1)
    points(df$t, df$estimation,   lty = 2,   col = 2)
    points(df2$t, df2$estimation, lty = 3,   col = 3)
    points(df3$t, df3$estimation, lty = 4,   col = 4)
    points(df4$t, df4$estimation, lty = 5,   col = 5)
    points(df5$t, df5$estimation, lty = 6,   col = 6)
    legend(5, 24, legend=c("alpha=1", "alpha=NW", "alpha=1.5", 
                           "alpha=1.25","log","alpha=1.1"),
           col=c(1,2,3,4,5,6), lty=1:6, cex=0.8, title = "diff alpha",
           text.font=4, bg='lightblue')
    plot(data$YEAR, data$SPEED,         type = 'l', col = 1, lty=1)
    points((data$YEAR)^1.1, data$SPEED, type = 'l', col = 2, lty=2)
    points((data$YEAR)^1.2, data$SPEED, type = 'l', col = 3, lty=3)
    points((data$YEAR)^1.3, data$SPEED, type = 'l', col = 4, lty=4)
    points((data$YEAR)^1.32, data$SPEED,type = 'l', col = 5, lty=5)
    points(sqrt(data$YEAR), data$SPEED, type = 'l', col = 6, lty=6)
    points(df$t, df$estimation,         type = 'l', col = 7, lty=7)
    legend(5, 24, legend=c("alpha=1", "alpha=1.1", "alpha=1.2", 
                  "alpha=1.3","alpha=1.32", "alpha=sqrt","alpha=NW"),
           col=1:7, lty=1:7, cex=0.8, title = "diff alpha",
           text.font=4, bg='lightblue')
    par(mfrow = c(1,1))
    library(ncvreg)
    Model = lm(SPEED ~ YEAR + DIST + I(YEAR^2) + I(DIST^2) 
               + I(YEAR*DIST), data = data, x = TRUE, y = TRUE)
    Risk = sum(Model$residuals^2)
    Model_step = step(Model, direction = "both", k = log(nrow(data)))
    Model_data = Model$x
    head(Model_data)
    Model_data = Model_data[,-1]
    head(Model_data)
    yy = Model$y
    fit1 <- ncvreg(Model_data,yy,penalty="SCAD")
    fit1
    plot(fit1$loss)
    ##############################################################
    Model2 = lm(SPEED ~ YEAR + DIST + I(YEAR^2) + I(DIST^2) +
                  I(YEAR*DIST), data = data, x = TRUE, y = TRUE)
    #VIF(Model2)
    Risk2 = sum(Model2$residuals^2)
    BIC = BIC(Model2)
    BIC
    Model_data2 = Model2$x
    Model_data2 = Model_data2[,-1]
    yyy = Model2$y
    cor_data = data.frame(speed = yyy, Model_data2)
    head(Model_data2)
    cor(cor_data)
    head(Model_data2)
    yy2 = data$SPEED
    fit2 <- ncvreg(Model_data2,yyy,penalty="SCAD")
    plot(fit2)
    
    Param = fit2$beta
    Lamb = fit2$lambda
    Data_Param = data.frame(Lamb = Lamb, YEAR = Param[2,], DIST = Param[3,],
                            YEARSQR = Param[4,], DISTSQR = Param[5,], 
                            YEAR_DIST = Param[6,])
    #, YEAR_AGE = Param[7,], DIST_AGE = Param[8,],
                            YEAR_STAGESWIN = Param[9,], DIST_STAGESWIN = Param[10,],
                            AGE_STAGESWIN = Param[11,], YEAR_DIST_AGE = Param[12,],
                            YEAR_DIST_STAGESWIN = Param[13,], 
                            YEAR_AGE_STAGESWIN = Param[14,],
                            DIST_AGE_STAGESWIN = Param[15,],
                            YEAR_DIST_AGE_STAGESWIN = Param[16,]
    )
    Data_Param  = Data_Param [order(Data_Param$Lamb, decreasing=TRUE), ]
    plot(Data_Param$Lamb, Data_Param[,2], col = 1, type = "l", xlab = expression(lambda), 
         ylab = expression(beta), lty = 1, ylim = c(-0.5,1),
         main = "Order of inclusion of parameters using SCAD")
    points(Data_Param$Lamb, Data_Param[,3], col = 2, lty = 1, type = "l")
    points(Data_Param$Lamb, Data_Param[,4], col = 3, lty = 1, type = "l")
    points(Data_Param$Lamb, Data_Param[,5], col = 4, lty = 1, type = "l")
    points(Data_Param$Lamb, Data_Param[,6], col = 5, lty = 1, type = "l")
    points(Data_Param$Lamb, Data_Param[,7], col = 6, lty = 1, type = "l")
    points(Data_Param$Lamb, Data_Param[,8], col = 7, lty = 1, type = "l")
    points(Data_Param$Lamb, Data_Param[,9], col = 8, lty = 1, type = "l")
    points(Data_Param$Lamb, Data_Param[,10], col = 9, lty = 1, type = "l")
    points(Data_Param$Lamb, Data_Param[,11], col = 10, lty = 1, type = "l")
    points(Data_Param$Lamb, Data_Param[,12], col = 11, lty = 1, type = "l")
    points(Data_Param$Lamb, Data_Param[,13], col = 12, lty = 1, type = "l")
    points(Data_Param$Lamb, Data_Param[,14], col = 13, lty = 1, type = "l")
    points(Data_Param$Lamb, Data_Param[,15], col = 14, lty = 1, type = "l")
    points(Data_Param$Lamb, Data_Param[,16], col = 15, lty = 1, type = "l")
    
    #points(Lamb, Param[10,], col = 9, lty = 1, type = "l")
    legend("topright", row.names(Param)[2:6], col = 1:5, lty = 1
    )
    
    ###############################################################
    # Empirical risk minimization
    ##############################################################
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
    
    
    ERM1 = Gaby(Loss = Risk2,eta = 0.05,
               n = nrow(data), h =4, m=10)
    ERM1
    ERM2 = ERM(Loss = Risk2,eta = 0.05,
               n = nrow(data), h =4, m=10)
    ERM2
    Bound = numeric()
    hhat = seq(from = 1, to = 100, by = 1)
    hhat1 = seq(from = 1, to = 100, by = 1)
    for (k in 1:length(hhat)) {
      hhat[k] = ERM(Loss = Risk2,eta = 0.05,
                            n = nrow(data), h =k, m=10)
      hhat1[k] = Gaby(Loss = Risk2,eta = 0.05,
                            n = nrow(data), h =k, m=10)
    }
    plot(seq(from = 1, to = 100, by = 1), hhat, col = 1)
    plot(seq(from = 1, to = 100, by = 1), hhat1, col = 2)
    Bound = ERM(fit1$loss,eta = 0.05,
                n = nrow(data), h =2, m=10)
    plot(1:length(fit1$loss),Bound, ylab = 'ERM', xlab = 'Index',
         main = paste("Upper Bound of ERM for each", " Model class",'\n'," using SCAD",sep = ""))
######################################################################
    require(grDevices)
    f = function(n,x){
      sqrt((x/n)*log(2*n*exp(1)/x)) 
    }
    
    x <- seq(100, 2000, length= 30)
    y <- x
    z <- outer(x, y, f)
    z[is.na(z)] <- 1
    op <- par(bg = "white")
    persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue")
    persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
          ltheta = 120, shade = 0.75, ticktype = "detailed",
          xlab = "h", ylab = "N", zlab = "RHS-2.37"
    ) -> res
    round(res, 3)
    ############################################################################
    # Cross validation
    ###########################################################################
    

    Model_cv = function(nparam, Size, n_folds, sigma, n_zero){
      set.seed(10)
      Param_beta = rnorm(nparam, mean = 5, sd = 3)
      Beta = matrix((c(Param_beta, rep(0,n_zero))), ncol = 1)
      data1 = datasim(Beta, mu=0, sigma, size = Size,seed=345)
      data1 = data.frame(scale(data1, scale = TRUE, center = TRUE))
      df = 2:ncol(data1) 
      folds_i <- sample(rep(1:n_folds, length.out = Size))
      cv_tmp <- matrix(NA, nrow = n_folds, ncol = length(df))
      
      for (k in 1:n_folds) {
        test_i <- which(folds_i == k)
        train_xy <- data1[-test_i, ]
        test_x <- data1[test_i, ]
        y = data1[test_i, ][,1]
        
        
        fitted_models <- apply(t(df), 2, function(degf) lm(y ~ ., data = train_xy[,1:degf]))
        
        pred <- mapply(function(obj, degf) predict(obj, test_x[, 1:degf]), 
                       fitted_models, df)
        cv_tmp[k, ] <- sapply(as.list(data.frame(pred)), function(y_hat) mean((y - y_hat)^2, na.rm = TRUE))
      }
      return(cv_tmp)
    }
    
    
    ###########################################################################
    Size_60 = 50:68
    dat_60_1 = matrix(c(
    100, 300983.8012, 297002.9979, 6562.445,
    100, 291446.0768, 287528.9701, 6546.303,
    100, 291409.6099, 287492.7487, 6552.765,
    100, 269741.8654, 265973.7111, 6504.855,
    100, 189473.1256, 186316.1848, 6262.231,
    100, 187152.9284, 184015.4210, 6260.084,
    100, 142518.5487, 139781.5457, 6074.169,
    100, 142062.4166, 139329.8086, 6078.454,
    100,  64558.1199,  62718.3890, 5526.228,
    100,   4835.8495,   4337.6109, 3661.755,
    57,    184.4617,    109.9579, 1058.243
    ), ncol = 4, byrow = TRUE)
    
    dat_60_2 = matrix(c(
    66, 189.2694, 109.5122, 1059.554,
    68, 190.4764, 109.5611, 1065.968,
    63, 186.8104, 108.9616, 1069.638,
    66, 188.6882, 109.0631, 1076.156,
    67, 188.9589, 108.8298, 1080.879,
    65, 187.7004, 108.7460, 1087.335,
    62, 185.7993, 108.6358, 1093.868,
    66, 188.2975, 108.7614, 1100.303,
    63, 186.3635, 108.6152, 1106.582,
    61, 184.7189, 108.2593, 1111.198), ncol = 4, byrow = TRUE
    )
    dat_60 = rbind(dat_60_1, dat_60_2)
    plot(dat_60[,2],xlab = "VCD", ylab = "ERM2")
    points(dat_60[,3],xlab = "VCD", ylab = "ERM2", col = 2)
    title(main = paste("ERM2~ VCD", " N = " ,600 , "\n"," NL = ","c(75,150,225,300,375,450,525,600)"  ))
    
    ######################################################################################################################
    # p = 15, N=400, NL = seq(from = 50, to = 400, by = 50)
    ##################################################################
    
    dat_15 = matrix(c(50, 102727.09446, 104980.1709, 3426.3679,
                      50,   74945.40604,  76872.4018, 3306.2229,
                      50,  53641.16581,  55274.1218, 3178.4223,
                      50,  46353.31935,  47872.5317, 3125.9961,
                      50,  21047.24401,  22076.6787, 2816.1058,
                      16,     66.86905,    110.6398,  493.0155,
                      25, 68.61495, 123.4442, 504.8823,
                      25, 68.30051, 123.0282, 508.8784,
                      26, 68.34595, 124.1460, 514.6572,
                      27, 68.36705, 125.2136, 520.2935,
                      25, 67.97364, 122.5954, 524.7679
    ), ncol = 4, byrow = TRUE)
    
    par(mfrow=c(2,1),oma=c(0,0,3,0),pch=42,font.sub=3,adj=0.5)
    
    par(mar = c(5,5,5,5))
    plot(10:20,dat_15[,2], pch=16, col = 1, xlab=NA, ylab= "ERM1 and ERM2", cex=1.2, lty = 1)
    points(10:20,dat_15[,3], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(10:20,dat_15[,1],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'VC Dimension')
    
    legend("topright", c("ERM1", "ERM2", "VCD"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),lwd = c(1,2,3),
           merge = TRUE, bg = "gray90")
    title(main = "Complexity Methods")
    ####################################
    # Plot for AIC, BIC and CV
    ###################################
    nparam = 15
    n_zero = 10
    Size = 400
    set.seed(10)
    Param_beta = rnorm(nparam, mean = 5, sd = 3)
    Beta = matrix((c(Param_beta, rep(0,n_zero))), ncol = 1)
    data1 = datasim(Beta, mu=0, sigma = 0.4, size = Size,seed=345)
    data1 = data.frame(scale(data1, scale = TRUE, center = TRUE))
    mat_aic_bic = matrix(NA, nrow = ncol(data1)-1, ncol = 2)
    colnames(mat_aic_bic) = c("AIC", "BIC")
    
    for(i in 2:ncol(data1)){
      Model = lm(y ~., data = data1[,1:i])
      mat_aic_bic[i-1,] = c(round(AIC(Model)), round(BIC(Model)))
    }
    
    cv_tmp = Model_cv(nparam = 15, Size = 400, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    out = data.frame(mat_aic_bic, cv = colMeans(cv_tmp))[10:20,]
    
    par(mar = c(5,5,5,5))
    plot(10:20,out[,1], pch=16, col = 1, xlab=NA, ylab= "AIC and BIC", cex=1.2, lty = 1)
    points(10:20,out[,2], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(10:20,out[,3],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'Cross Validation')
    
    legend("topright", c("AIC", "BIC", "CV"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),
           merge = TRUE, bg = "gray90")
    title(main = "traditional Methods")
     
    mtext(side=3,outer=T,line=0,  paste("p = ", 15, ","," n = ",400 , ",", "\n", "n_l", " = ","50,100,150,200,250,300,350,400"),cex=1.3) 
    
  
    #par(mfrow = c(3,2))
    #plot(10:20, cv, lwd = 2, col = gray(0.4), ylab = "Prediction error", ylim = c(0,0.7),
     #    xlab = "Size of the conjectured Model ", main = paste0(n_folds, 
      #                                                          "-fold Cross-Validation", "\n",
       #                                                         "p = ", 15,","," N = " ,400))
    #lines(10:20, cv, lwd = 2, col = "darkred", lty = 2)
    #cv_sd <- apply(cv_tmp, 2, sd)/sqrt(n_folds)
    #errbar(10:20, cv, cv + cv_sd, cv - cv_sd, add = TRUE, col = "steelblue2", pch = 19, 
     #      lwd = 0.5)
    ##################################################################
    # p = 30, N = 400, NL = 50, 100, 150, 200, 250, 300, 350, 400
    #################################################################
    dat30 = matrix(c(50, 116360.87, 118757.68, 3572.08,
                     50, 102187.23, 104434.42, 3526.1150,
                     50, 90273.38, 92386.57, 3482.52,
                     50, 6583.36, 7166.91, 2440.7975,
                     30, 66.30, 125.44, 571.34,
                     42, 82.38, 150.10, 684.47,
                     39, 82.21, 147.48, 690.63, 
                     40, 82.09, 148.14, 695.71,
                     41, 81.84, 148.5, 699.91,
                     42, 81.58, 149.03, 704.01,
                     38, 80.61, 144.55, 705.24,
                     38, 79.60, 143.22, 704.67,
                     40, 79.61, 144.85 , 710.20), ncol = 4, byrow = TRUE)
    par( mfrow = c(1,1) )
    
    par(mfrow=c(2,1),oma=c(0,0,3,0),pch=42,font.sub=3,adj=0.5)
    
    par(mar = c(5,5,5,5))
    plot(26:38,dat30[,2], pch=16, col = 1, xlab=NA, ylab= "ERM1 and ERM2", cex=1.2, lty = 1)
    points(26:38,dat30[,3], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(26:38,dat30[,1],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'VC Dimension')
    
    legend("topright", c("ERM1", "ERM2", "VCD"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),lwd = c(1,2,3),
           merge = TRUE, bg = "gray90")
    title(main = "Complexity Methods")
    ####################################
    # Plot for AIC, BIC and CV
    ###################################
    nparam = 30
    n_zero = 10
    Size = 400
    set.seed(10)
    Param_beta = rnorm(nparam, mean = 5, sd = 3)
    Beta = matrix((c(Param_beta, rep(0,n_zero))), ncol = 1)
    data1 = datasim(Beta, mu=0, sigma = 0.4, size = Size,seed=345)
    data1 = data.frame(scale(data1, scale = TRUE, center = TRUE))
    mat_aic_bic = matrix(NA, nrow = ncol(data1)-1, ncol = 2)
    colnames(mat_aic_bic) = c("AIC", "BIC")
    
    for(i in 2:ncol(data1)){
      Model = lm(y ~., data = data1[,1:i])
      mat_aic_bic[i-1,] = c(round(AIC(Model)), round(BIC(Model)))
    }
    
    cv_tmp = Model_cv(nparam = 30, Size = 400, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    out = data.frame(mat_aic_bic, cv = colMeans(cv_tmp))[26:38,]
    
    par(mar = c(5,5,5,5))
    plot(26:38,out[,1], pch=16, col = 1, xlab=NA, ylab= "AIC and BIC", cex=1.2, lty = 1)
    points(26:38,out[,2], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(26:38,out[,3],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'Cross Validation')
    
    legend("topright", c("AIC", "BIC", "CV"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),
           merge = TRUE, bg = "gray90")
    title(main = "traditional Methods")
    
    mtext(side=3,outer=T,line=0,paste("p = ", 30,","," n = " ,400 , ",", "\n","n_l = ","50,100,150,200,250,300,350,400"),cex=1.3) 
    
    
    #plot(26:38,dat301[,2],xlab = "Size of the nested model", ylab = "ERM2", col = 1, lty = 1, pch = 1, lwd = 1, ylim = c(-1,2))
    #points(26:38,dat301[,3],xlab = "VCD", ylab = "ERM2", col = 2, lty = 2, pch = 2, lwd = 2)
    #points(26:38,dat301[,4], col = 3, lty = 3, pch = 3, lwd = 3)
    #points(30:50,1000*dat_40[,1], col = 4, lty = 4, pch = 4, lwd = 4)
    #legend("topright", c("ERM2", "ERM1", "BIC"), col = c(1,2,3),
    #       text.col = "green4", lty = c(1,2,3), pch = c(1,2,3),lwd = c(1,2,3),
    #       merge = TRUE, bg = "gray90")
    #title(main = paste("ERM1, ERM2, BIC~ VCD,", "p = ", 30,",N = " ,400 , "\n",",NL = ","50,100,150,200,250,300,350,400"  ))
    #plot(26:38,dat30[,1], xlab = "Size of the nested model", ylab = "Estimated Hhat", col = 4, lty = 4, pch = 4, lwd = 4)
    #title(main = paste("p = ", 30,",","N = " ,400 , ",", "\n","NL = ","50,100,150,200,250,300,350,400"  ))
    
    #cv_tmp = Model_cv(nparam = 30, Size = 400, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    #cv_tmp = cv_tmp[,26:38]
    #cv = colMeans(cv_tmp)
    #plot(26:38, cv, lwd = 2, col = gray(0.4), ylab = "Prediction error", ylim = c(0,0.11),
    #     xlab = "Size of the conjectured Model ", main = paste0(n_folds, 
    #                                                            "-fold Cross-Validation", "\n",
    #                                                            "p = ", 30,","," N = " ,400 ))
    #lines(26:38, cv, lwd = 2, col = "darkred", lty = 2)
    #cv_sd <- apply(cv_tmp, 2, sd)/sqrt(n_folds)
    #errbar(26:38, cv, cv + cv_sd, cv - cv_sd, add = TRUE, col = "steelblue2", pch = 19, 
    #       lwd = 0.5)
    ##########################################################################
    # p = 40, N = 600, NL = 75-600 by 75
    #########################################################################
    
    dat_40 = matrix(c( 75, 110058.2246, 107763.00054, 5021.8454,
                       75, 109865.7575, 107572.54731, 5027.1809,
                       75,  96154.4606,  94009.56177, 4952.7108,
                       75,  58241.7034,  56573.92385, 4654.3695,
                       75,  39611.8789,  38237.67853, 4425.6950,
                       75,  39158.7359,  37792.45811, 4425.0636,
                       75,  38860.0620,  37499.03115, 4426.7830,
                       75,  16189.8521,  15313.82068, 3895.6834,
                       75,  16144.8663,  15270.06251, 3900.3626,
                       75,   6119.0625,   5583.16753, 3302.6033,
                       41,    163.4084,     97.90365,  848.5549,
                       48,  167.7513, 97.83740, 852.3346,
                       50,  169.2108, 97.88200, 858.4368,
                       46,  166.1366, 97.69479, 864.8002,
                       48, 167.6096, 97.72774, 870.8093,
                       50, 169.1004, 97.79682, 877.0713,
                       47, 166.8030, 97.65402, 883.4210,
                       45, 165.2236, 97.54999, 889.7503,
                       47, 166.4983, 97.41791, 894.6697,
                       45, 164.9074, 97.30405, 900.9345,
                       48, 167.0091, 97.26315, 906.1478
                       
    ), ncol = 4, byrow = TRUE)
    par( mfrow = c(1,1) )
    #dat_40 = scale(dat_40)
    par(mfrow=c(2,1),oma=c(0,0,3,0),pch=42,font.sub=3,adj=0.5)
    
    par(mar = c(5,5,5,5))
    plot(30:50,dat_40[,2], pch=16, col = 1, xlab=NA, ylab= "ERM1 and ERM2", cex=1.2, lty = 1)
    points(30:50,dat_40[,3], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(30:50,dat_40[,1],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'VC Dimension')
    
    legend("topright", c("ERM1", "ERM2", "VCD"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),lwd = c(1,2,3),
           merge = TRUE, bg = "gray90")
    title(main = "Complexity Methods")
    ####################################
    # Plot for AIC, BIC and CV
    ###################################
    nparam = 40
    n_zero = 10
    Size = 600
    set.seed(10)
    Param_beta = rnorm(nparam, mean = 5, sd = 3)
    Beta = matrix((c(Param_beta, rep(0,n_zero))), ncol = 1)
    data1 = datasim(Beta, mu=0, sigma = 0.4, size = Size,seed=345)
    data1 = data.frame(scale(data1, scale = TRUE, center = TRUE))
    mat_aic_bic = matrix(NA, nrow = ncol(data1)-1, ncol = 2)
    colnames(mat_aic_bic) = c("AIC", "BIC")
    
    for(i in 2:ncol(data1)){
      Model = lm(y ~., data = data1[,1:i])
      mat_aic_bic[i-1,] = c(round(AIC(Model)), round(BIC(Model)))
    }
    
    cv_tmp = Model_cv(nparam = 40, Size = 600, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    out = data.frame(mat_aic_bic, cv = colMeans(cv_tmp))[30:50,]
    
    par(mar = c(5,5,5,5))
    plot(30:50,out[,1], pch=16, col = 1, xlab=NA, ylab= "AIC and BIC", cex=1.2, lty = 1)
    points(30:50,out[,2], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(30:50,out[,3],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'Cross Validation')
    
    legend("topright", c("AIC", "BIC", "CV"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),
           merge = TRUE, bg = "gray90")
    title(main = "Traditional Methods")
    
    mtext(side=3,outer=T,line=0,paste("p = ", 40,","," n = " ,600 , ",", "n_l = ","75-600 by 75"),cex=1.3) 
    
    
    #plot(30:50,dat_40[,2],xlab = "Size of the nested model", ylab = "ERM2", col = 1, lty = 1, pch = 1, lwd = 1, ylim = c(-1,2))
    #points(30:50,dat_40[,3],xlab = "VCD", ylab = "ERM2", col = 2, lty = 2, pch = 2, lwd = 2)
    #points(30:50,dat_40[,4], col = 3, lty = 3, pch = 3, lwd = 3)
    #points(30:50,1000*dat_40[,1], col = 4, lty = 4, pch = 4, lwd = 4)
    #legend("topright", c("ERM2", "ERM1", "BIC"), col = c(1,2,3),
    #       text.col = "green4", lty = c(1,2,3), pch = c(1,2,3),lwd = c(1,2,3),
    #       merge = TRUE, bg = "gray90")
    #title(main = paste("ERM1, ERM2, BIC~ VCD,", "p = ", 40,",N = " ,600 , "\n",",NL = ","c(75,150,225,300,375,450,525,600)"  ))
    #plot(30:50,dat_40[,1], xlab = "Size of the nested model", ylab = "Estimated Hhat", col = 4, lty = 4, pch = 4, lwd = 4)
    #title(main = paste("p = ", 40,",","N = " ,600 , ",", "\n","NL = ","c(75,150,225, 300, 375, 450, 525, 600)"  ))
    
    #cv_tmp = Model_cv(nparam = 40, Size = 600, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    #cv_tmp = cv_tmp[,30:50]
    #cv = colMeans(cv_tmp)
    #plot(30:50, cv, lwd = 2, col = gray(0.4), ylab = "Prediction error", ylim = c(0,0.25),
    #     xlab = "Size of the conjectured Model ", main = paste0(n_folds, 
    #                                                            "-fold Cross-Validation", "\n",
    #                                                            "p = ", 40,","," N = " ,600))
    #lines(30:50, cv, lwd = 2, col = "darkred", lty = 2)
    #require(Hmisc)
    #cv_sd <- apply(cv_tmp, 2, sd)/sqrt(n_folds)
    #errbar(30:50, cv, cv + cv_sd, cv - cv_sd, add = TRUE, col = "steelblue2", pch = 19, 
    #       lwd = 0.5)
    ######################################################################
    # p = 50, N = 600, NL = 75 to 600 by 75
    ####################################################################
    dat50 = matrix(c(75, 63643.52, 65411.39, 4808.19,
                     75, 44226.92, 45703.52, 4596.18,
                     75, 34270.58, 35572.47, 4449.52,
                     75, 29017.58, 30216.93,  4356.06,
                     48,  97.26, 167.71,  906.15,
                     55, 97.28, 172.60, 911.16,
                     53, 97.40, 171.16, 917.57,
                     52, 97.36,  170.43,  923.96,
                     52, 96.94, 170.16,  929.03,
                     52, 96.86, 169.99,   934.57,
                     50, 96.65, 168.49,  940.81,
                     51, 96.65, 169.22,  947.18), ncol = 4, byrow = TRUE)
    par( mfrow = c(1,1) )
    #dat501 = scale(dat50)
    
    
    par(mfrow=c(2,1),oma=c(0,0,3,0),pch=42,font.sub=3,adj=0.5)
    
    par(mar = c(5,5,5,5))
    plot(46:57,dat50[,2], pch=16, col = 1, xlab=NA, ylab= "ERM1 and ERM2", cex=1.2, lty = 1)
    points(46:57,dat50[,3], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(46:57,dat50[,1],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'VC Dimension')
    
    legend("topright", c("ERM1", "ERM2", "VCD"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),lwd = c(1,2,3),
           merge = TRUE, bg = "gray90")
    title(main = "Complexity Methods")
    ####################################
    # Plot for AIC, BIC and CV
    ###################################
    nparam = 50
    n_zero = 10
    Size = 600
    set.seed(10)
    Param_beta = rnorm(nparam, mean = 5, sd = 3)
    Beta = matrix((c(Param_beta, rep(0,n_zero))), ncol = 1)
    data1 = datasim(Beta, mu=0, sigma = 0.4, size = Size,seed=345)
    data1 = data.frame(scale(data1, scale = TRUE, center = TRUE))
    mat_aic_bic = matrix(NA, nrow = ncol(data1)-1, ncol = 2)
    colnames(mat_aic_bic) = c("AIC", "BIC")
    
    for(i in 2:ncol(data1)){
      Model = lm(y ~., data = data1[,1:i])
      mat_aic_bic[i-1,] = c(round(AIC(Model)), round(BIC(Model)))
    }
    
    cv_tmp = Model_cv(nparam = 50, Size = 600, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    out = data.frame(mat_aic_bic, cv = colMeans(cv_tmp))[46:57,]
    
    par(mar = c(5,5,5,5))
    plot(46:57,out[,1], pch=16, col = 1, xlab=NA, ylab= "AIC and BIC", cex=1.2, lty = 1)
    points(46:57,out[,2], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(46:57,out[,3],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'Cross Validation')
    
    legend("topright", c("AIC", "BIC", "CV"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),
           merge = TRUE, bg = "gray90")
    title(main = "Traditional Methods")
    
    mtext(side=3,outer=T,line=0,paste("p = ", 50,","," n = " ,600 , ",", "n_l = ","75-600 by 75"),cex=1.3) 
    
    
    #plot(46:57,dat501[,2],xlab = "Size of the nested model", ylab = "ERM2", col = 1, lty = 1, pch = 1, lwd = 1, ylim = c(-1,2))
    #points(46:57,dat501[,3],xlab = "VCD", ylab = "ERM2", col = 2, lty = 2, pch = 2, lwd = 2)
    #points(46:57,dat501[,4], col = 3, lty = 3, pch = 3, lwd = 3)
    #points(30:50,1000*dat_40[,1], col = 4, lty = 4, pch = 4, lwd = 4)
    #legend("topright", c("ERM2", "ERM1", "BIC"), col = c(1,2,3),
    #       text.col = "green4", lty = c(1,2,3), pch = c(1,2,3),lwd = c(1,2,3),
    #       merge = TRUE, bg = "gray90")
    #title(main = paste("ERM1, ERM2, BIC~ VCD,", "p = ", 50,",N = " ,600 , "\n",",NL = ","75,150,225,300,375,450,525,600"  ))
    #plot(46:57,dat50[,1], xlab = "Size of the nested model", ylab = "Estimated Hhat", col = 4, lty = 4, pch = 4, lwd = 4, ylim = c(40,90))
    #title(main = paste("p = ", 50,",","N = " ,600 , ",", "\n","NL = ","75,150,225,300,375,450,525,600"  ))
    
    #cv_tmp = Model_cv(nparam = 50, Size = 600, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    #cv_tmp = cv_tmp[,46:57]
    #cv = colMeans(cv_tmp)
    #plot(46:57, cv, lwd = 2, col = gray(0.4), ylab = "Prediction error", ylim = c(0,0.12),
    #     xlab = "Size of the conjectured Model ", main = paste0(n_folds, 
    #                                                            "-fold Cross-Validation", "\n",
    #                                                            "p = ", 50,","," N = " ,600 ))
    #lines(46:57, cv, lwd = 2, col = "darkred", lty = 2)
    #cv_sd <- apply(cv_tmp, 2, sd)/sqrt(n_folds)
    #errbar(46:57, cv, cv + cv_sd, cv - cv_sd, add = TRUE, col = "steelblue2", pch = 19, 
    #       lwd = 0.5)
    ###################################################################
    # p = 70, N = 2000, NL = 500, 700, 1000, 1500, 2000
    ##############################################################
    dat_70_2000 = matrix(c(500,       50000,       50000, 17620.27,
                           500,       50000,       50000, 17617.12,
                           500,       50000,       50000, 17536.70,
                           500,       50000,       50000, 17501.47,
                           500,       50000,       50000, 17240.56,
                           500,       50000,       50000, 16629.22,
                           500,       50000,       50000, 15606.13,
                           500,       50000,       50000, 14349.52,
                           500,       50000,       50000, 13875.20,
                           71, 3739.5619, 3483.2204,  7322.59,
                           78, 389.9254, 307.0259, 2453.663,
                           79, 390.4185, 307.0479, 2461.263,
                           79, 390.3845, 307.0177, 2468.665,
                           79, 390.3838, 307.0171, 2476.262,
                           80, 390.8566, 307.0241, 2483.763,
                           80, 390.7139, 306.8975, 2490.527,
                           80, 390.6939, 306.8798, 2498.011,
                           80, 389.9303, 306.2024, 2501.127,
                           80,  389.7455, 306.0386, 2507.642,
                           80, 389.7455, 306.0386, 2515.243,
                           71,  385.8239,  306.3781,  2442.84
    ), ncol = 4, byrow = TRUE)
    par( mfrow = c(1,1) )
    #dat_70_20001 = scale(dat_70_2000)
    
    par(mfrow=c(2,1),oma=c(0,0,3,0),pch=42,font.sub=3,adj=0.5)
    
    par(mar = c(5,5,5,5))
    plot(60:80,dat_70_2000[,2], pch=16, col = 1, xlab=NA, ylab= "ERM1 and ERM2", cex=1.2, lty = 1)
    points(60:80,dat_70_2000[,3], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(60:80,dat_70_2000[,1],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'VC Dimension')
    
    legend("topright", c("ERM1", "ERM2", "VCD"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),lwd = c(1,2,3),
           merge = TRUE, bg = "gray90")
    title(main = "Complexity Methods")
    ####################################
    # Plot for AIC, BIC and CV
    ###################################
    nparam = 70
    n_zero = 10
    Size = 2000
    set.seed(10)
    Param_beta = rnorm(nparam, mean = 5, sd = 3)
    Beta = matrix((c(Param_beta, rep(0,n_zero))), ncol = 1)
    data1 = datasim(Beta, mu=0, sigma = 0.4, size = Size,seed=345)
    data1 = data.frame(scale(data1, scale = TRUE, center = TRUE))
    mat_aic_bic = matrix(NA, nrow = ncol(data1)-1, ncol = 2)
    colnames(mat_aic_bic) = c("AIC", "BIC")
    
    for(i in 2:ncol(data1)){
      Model = lm(y ~., data = data1[,1:i])
      mat_aic_bic[i-1,] = c(round(AIC(Model)), round(BIC(Model)))
    }
    
    cv_tmp = Model_cv(nparam = 70, Size = 2000, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    out = data.frame(mat_aic_bic, cv = colMeans(cv_tmp))[60:80,]
    
    par(mar = c(5,5,5,5))
    plot(60:80,out[,1], pch=16, col = 1, xlab=NA, ylab= "AIC and BIC", cex=1.2, lty = 1)
    points(60:80,out[,2], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(60:80,out[,3],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'Cross Validation')
    
    legend("topright", c("AIC", "BIC", "CV"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),
           merge = TRUE, bg = "gray90")
    title(main = "Traditional Methods")
    
    mtext(side=3,outer=T,line=0,paste("p = ", 70,","," n = " ,2000 , ",", "n_l = ","500, 700, 1000, 1500, 2000"),cex=1.3) 
    
    
    #plot(60:80,dat_70_20001[,2],xlab = "Size of the nested model", ylab = "ERM2", col = 1, lty = 1, pch = 1, lwd = 1, ylim = c(-1,2))
    #points(60:80,dat_70_20001[,3],xlab = "VCD", ylab = "ERM2", col = 2, lty = 2, pch = 2, lwd = 2)
    #points(60:80,dat_70_20001[,4], col = 3, lty = 3, pch = 3, lwd = 3)
    #points(30:50,1000*dat_40[,1], col = 4, lty = 4, pch = 4, lwd = 4)
    #legend("topright", c("ERM2", "ERM1", "BIC"), col = c(1,2,3),
    #       text.col = "green4", lty = c(1,2,3), pch = c(1,2,3),lwd = c(1,2,3),
    #       merge = TRUE, bg = "gray90")
    #title(main = paste("ERM1, ERM2, BIC~ VCD,", "p = ", 70,",N = " ,2000 , "\n",",NL = ","500,700,1000,1500,2000"  ))
    #plot(60:80,dat_70_2000[,1], xlab = "Size of the nested model", ylab = "Estimated Hhat", col = 4, lty = 4, pch = 4, lwd = 4, ylim = c(30,500))
    #title(main = paste("p = ", 70,",","N = " ,2000 , ",", "\n","NL = ","500,700,1000,1500,2000"  ))
    
    #cv_tmp = Model_cv(nparam = 70, Size = 2000, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    #cv_tmp = cv_tmp[,60:80]
    #cv = colMeans(cv_tmp)
    #plot(60:80, cv, lwd = 2, col = gray(0.4), ylab = "Prediction error", ylim = c(0,0.20),
    #     xlab = "Size of the conjectured Model ", main = paste0(n_folds, 
    #                                                            "-fold Cross-Validation", "\n",
    #                                                            "p = ", 70,","," N = " ,2000))
    #lines(60:80, cv, lwd = 2, col = "darkred", lty = 2)
    #cv_sd <- apply(cv_tmp, 2, sd)/sqrt(n_folds)
    #errbar(60:80, cv, cv + cv_sd, cv - cv_sd, add = TRUE, col = "steelblue2", pch = 19, 
    #       lwd = 0.5)
    ####################################################################
    # p = 70 , N = 700, NL = 100 to 700 by 100
    ####################################################################
    dat70 = matrix(c(100, 220239.0913, 216834.9130, 6407.725,
                     100, 219093.8419, 215698.5450, 6410.598,
                     100, 209775.2781, 206453.1268, 6386.482,
                     100, 205426.1843, 202138.7269, 6378.249,
                     100, 179857.2140, 176781.6118, 6290.970,
                     100, 132655.3773, 130015.0371, 6082.425,
                     100,  79834.3988,  77787.7338, 5729.383,
                     100,  42960.1970,  41460.7738, 5295.412,
                     100,  33927.0741,  32595.3923, 5133.526,
                      63,   1486.4357,   1255.6499, 2857.236,
                      61,    184.7189,    108.2593, 1111.198,
                     68, 187.7582, 107.4676, 1110.587,
                     67, 187.1443, 107.4311, 1117.123,
                     66, 186.3614, 107.2675, 1122.782,
                     70, 188.7221, 107.3475, 1128.940,
                     66, 186.1702, 107.1200, 1134.862,
                     68, 187.3888, 107.1835, 1141.374,
                     67, 186.7367, 107.1172, 1147.703,
                     65, 185.3436, 106.9256, 1153.408,
                     64, 184.6342, 106.8247, 1159.504,
                     64, 184.3821, 106.6300, 1164.700
                     
    ), ncol = 4, byrow = TRUE)
    par( mfrow = c(1,1) )
    #dat701 = scale(dat70)
    
    par(mfrow=c(2,1),oma=c(0,0,3,0),pch=42,font.sub=3,adj=0.5)
    
    par(mar = c(5,5,5,5))
    plot(60:80,dat70[,2], pch=16, col = 1, xlab=NA, ylab= "ERM1 and ERM2", cex=1.2, lty = 1)
    points(60:80,dat70[,3], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(60:80,dat70[,1],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'VC Dimension')
    
    legend("topright", c("ERM1", "ERM2", "VCD"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),lwd = c(1,2,3),
           merge = TRUE, bg = "gray90")
    title(main = "Complexity Methods")
    ####################################
    # Plot for AIC, BIC and CV
    ###################################
    nparam = 70
    n_zero = 10
    Size = 700
    set.seed(10)
    Param_beta = rnorm(nparam, mean = 5, sd = 3)
    Beta = matrix((c(Param_beta, rep(0,n_zero))), ncol = 1)
    data1 = datasim(Beta, mu=0, sigma = 0.4, size = Size,seed=345)
    data1 = data.frame(scale(data1, scale = TRUE, center = TRUE))
    mat_aic_bic = matrix(NA, nrow = ncol(data1)-1, ncol = 2)
    colnames(mat_aic_bic) = c("AIC", "BIC")
    
    for(i in 2:ncol(data1)){
      Model = lm(y ~., data = data1[,1:i])
      mat_aic_bic[i-1,] = c(round(AIC(Model)), round(BIC(Model)))
    }
    
    cv_tmp = Model_cv(nparam = 70, Size = 700, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    out = data.frame(mat_aic_bic, cv = colMeans(cv_tmp))[60:80,]
    
    par(mar = c(5,5,5,5))
    plot(60:80,out[,1], pch=16, col = 1, xlab=NA, ylab= "AIC and BIC", cex=1.2, lty = 1)
    points(60:80,out[,2], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(60:80,out[,3],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'Cross Validation')
    
    legend("topright", c("AIC", "BIC", "CV"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),
           merge = TRUE, bg = "gray90")
    title(main = "Traditional Methods")
    
    mtext(side=3,outer=T,line=0,paste("p = ", 70,","," n = " ,700 , ",", "n_l = ","100 to 700 by 100"),cex=1.3) 
    
    #plot(60:80,dat701[,2],xlab = "Size of the nested model", ylab = "ERM2", col = 1, lty = 1, pch = 1, lwd = 1, ylim = c(-1,2))
    #points(60:80,dat701[,3],xlab = "VCD", ylab = "ERM2", col = 2, lty = 2, pch = 2, lwd = 2)
    #points(60:80,dat701[,4], col = 3, lty = 3, pch = 3, lwd = 3)
    #points(30:50,1000*dat_40[,1], col = 4, lty = 4, pch = 4, lwd = 4)
    #legend("topright", c("ERM2", "ERM1", "BIC"), col = c(1,2,3),
    #       text.col = "green4", lty = c(1,2,3), pch = c(1,2,3),lwd = c(1,2,3),
    #       merge = TRUE, bg = "gray90")
    #title(main = paste("ERM1, ERM2, BIC~ VCD,", "p = ", 70,",N = " ,700 , "\n",",NL = ","100,200,300,400,500,600,700"  ))
    #plot(60:80,dat70[,1], xlab = "Size of the nested model", ylab = "Estimated Hhat", col = 4, lty = 4, pch = 4, lwd = 4, ylim = c(50,100))
    #title(main = paste("p = ", 70,",","N = " ,700 , ",", "\n","NL = ","100,200,300,400,500,600,700"  ))
    
    #cv_tmp = Model_cv(nparam = 70, Size = 700, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    #cv_tmp = cv_tmp[,60:80]
    #cv = colMeans(cv_tmp)
    #plot(60:80, cv, lwd = 2, col = gray(0.4), ylab = "Prediction error", ylim = c(0,0.25),
    #     xlab = "Size of the conjectured Model ", main = paste0(n_folds, 
     #                                                           "-fold Cross-Validation", "\n",
    #                                                            "p = ", 70,","," N = " ,700))
    #lines(60:80, cv, lwd = 2, col = "darkred", lty = 2)
    #cv_sd <- apply(cv_tmp, 2, sd)/sqrt(n_folds)
    #errbar(60:80, cv, cv + cv_sd, cv - cv_sd, add = TRUE, col = "steelblue2", pch = 19, 
    #       lwd = 0.5)
    
    
    ##################################################################
    # p = 60, N = 600, NL = 
    #################################################################
    dat60 = matrix(c(75, 55967.97, 57626.89,   4795.03,
                     75, 52809.85, 54421.78,  4766.58,
                     75, 26574.66, 27723.16, 4360.86,
                     75, 20887.67, 21907.86, 4222.73,
                     48, 96.05, 166.14, 962.10,
                     53, 95.72, 169.13, 964.82,
                     51, 95.62, 167.65, 971.13,
                     50, 95.57, 166.90, 977.47,
                     48, 95.48, 165.39,  983.85,
                     44, 95.29, 162.29,   990.25,
                     75, 96.41, 183.47,   995.92,
                     75, 96.41, 183.47,  1002.32), ncol = 4, byrow = TRUE)
    par( mfrow = c(1,1) )
    #dat601 = scale(dat60)
    
    par(mfrow=c(2,1),oma=c(0,0,3,0),pch=42,font.sub=3,adj=0.5)
    
    par(mar = c(5,5,5,5))
    plot(56:67,dat60[,2], pch=16, col = 1, xlab=NA, ylab= "ERM1 and ERM2", cex=1.2, lty = 1)
    points(56:67,dat60[,3], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(56:67,dat60[,1],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'VC Dimension')
    
    legend("topright", c("ERM1", "ERM2", "VCD"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),lwd = c(1,2,3),
           merge = TRUE, bg = "gray90")
    title(main = "Complexity Methods")
    ####################################
    # Plot for AIC, BIC and CV
    ###################################
    nparam = 60
    n_zero = 10
    Size = 600
    set.seed(10)
    Param_beta = rnorm(nparam, mean = 5, sd = 3)
    Beta = matrix((c(Param_beta, rep(0,n_zero))), ncol = 1)
    data1 = datasim(Beta, mu=0, sigma = 0.4, size = Size,seed=345)
    data1 = data.frame(scale(data1, scale = TRUE, center = TRUE))
    mat_aic_bic = matrix(NA, nrow = ncol(data1)-1, ncol = 2)
    colnames(mat_aic_bic) = c("AIC", "BIC")
    
    for(i in 2:ncol(data1)){
      Model = lm(y ~., data = data1[,1:i])
      mat_aic_bic[i-1,] = c(round(AIC(Model)), round(BIC(Model)))
    }
    
    cv_tmp = Model_cv(nparam = 60, Size = 600, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    out = data.frame(mat_aic_bic, cv = colMeans(cv_tmp))[56:67,]
    
    par(mar = c(5,5,5,5))
    plot(56:67,out[,1], pch=16, col = 1, xlab=NA, ylab= "AIC and BIC", cex=1.2, lty = 1)
    points(56:67,out[,2], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(56:67,out[,3],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'Cross Validation')
    
    legend("topright", c("AIC", "BIC", "CV"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),
           merge = TRUE, bg = "gray90")
    title(main = "Traditional Methods")
    
    mtext(side=3,outer=T,line=0,paste("p = ", 60,","," n = " ,600 , ",", "n_l = ","75 to 600 by 75"),cex=1.3)
    
   # plot(56:67,dat601[,2],xlab = "Size of the nested model", ylab = "ERM2", col = 1, lty = 1, pch = 1, lwd = 1, ylim = c(-1,2))
    #points(56:67,dat601[,3],xlab = "VCD", ylab = "ERM2", col = 2, lty = 2, pch = 2, lwd = 2)
    #points(56:67,dat601[,4], col = 3, lty = 3, pch = 3, lwd = 3)
    #points(30:50,1000*dat_40[,1], col = 4, lty = 4, pch = 4, lwd = 4)
    #legend("topright", c("ERM2", "ERM1", "BIC"), col = c(1,2,3),
    #       text.col = "green4", lty = c(1,2,3), pch = c(1,2,3),lwd = c(1,2,3),
    #       merge = TRUE, bg = "gray90")
    #title(main = paste("ERM1, ERM2, BIC~ VCD,", "p = ", 60,",N = " ,600 , "\n",",NL = ","75 to 600 by 75"  ))
    #plot(56:67,dat60[,1], xlab = "Size of the nested model", ylab = "Estimated Hhat", col = 4, lty = 4, pch = 4, lwd = 4)
    #title(main = paste("p = ", 60,",","N = " ,600 , ",", "\n","NL = ","75 to 600 by 75"  ))
    
    #cv_tmp = Model_cv(nparam = 60, Size = 600, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    #cv_tmp = cv_tmp[,56:67]
    #cv = colMeans(cv_tmp)
    #plot(56:67, cv, lwd = 2, col = gray(0.4), ylab = "Prediction error", ylim = c(0,0.08),
    #     xlab = "Size of the conjectured Model ", main = paste0(n_folds, 
    #                                                            "-fold Cross-Validation", "\n",
    #                                                            "p = ", 60,","," N = " ,600))
    #lines(56:67, cv, lwd = 2, col = "darkred", lty = 2)
    #cv_sd <- apply(cv_tmp, 2, sd)/sqrt(n_folds)
    #errbar(56:67, cv, cv + cv_sd, cv - cv_sd, add = TRUE, col = "steelblue2", pch = 19, 
    #       lwd = 0.5)
    ###################################################################
    # p = 60, N = 2000, NL = 500, 700, 1000, 1500, 2000
    #################################################################
    dat60_2000 = matrix(c(500, 50000, 50000, 15177.13,
                          500, 50000, 50000, 15068.43,
                          500, 50000, 50000, 13705.88,
                          500, 50000, 50000, 13234.45,
                          59, 307.79, 381.07, 2377.64,
                          69, 307.65, 386.20, 2382.83,
                          69, 307.30, 385.80, 2388.11,
                          69, 307.17, 385.66, 2394.87,
                          70, 307.12, 386.14,  2402.16,
                          70, 307.12, 386.13,  2409.73), ncol = 4, byrow = TRUE)
    par( mfrow = c(1,1) )
    #dat60_20001 = scale(dat60_2000)
    
    par(mfrow=c(2,1),oma=c(0,0,3,0),pch=42,font.sub=3,adj=0.5)
    
    par(mar = c(5,5,5,5))
    plot(56:65,dat60_2000[,2], pch=16, col = 1, xlab=NA, ylab= "ERM1 and ERM2", cex=1.2, lty = 1)
    points(56:65,dat60_2000[,3], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(56:65,dat60_2000[,1],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'VC Dimension')
    
    legend("topright", c("ERM1", "ERM2", "VCD"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),lwd = c(1,2,3),
           merge = TRUE, bg = "gray90")
    title(main = "Complexity Methods")
    ####################################
    # Plot for AIC, BIC and CV
    ###################################
    nparam = 60
    n_zero = 10
    Size = 2000
    set.seed(10)
    Param_beta = rnorm(nparam, mean = 5, sd = 3)
    Beta = matrix((c(Param_beta, rep(0,n_zero))), ncol = 1)
    data1 = datasim(Beta, mu=0, sigma = 0.4, size = Size,seed=345)
    data1 = data.frame(scale(data1, scale = TRUE, center = TRUE))
    mat_aic_bic = matrix(NA, nrow = ncol(data1)-1, ncol = 2)
    colnames(mat_aic_bic) = c("AIC", "BIC")
    
    for(i in 2:ncol(data1)){
      Model = lm(y ~., data = data1[,1:i])
      mat_aic_bic[i-1,] = c(round(AIC(Model)), round(BIC(Model)))
    }
    
    cv_tmp = Model_cv(nparam = 60, Size = 2000, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    out = data.frame(mat_aic_bic, cv = colMeans(cv_tmp))[56:65,]
    
    par(mar = c(5,5,5,5))
    plot(56:65,out[,1], pch=16, col = 1, xlab=NA, ylab= "AIC and BIC", cex=1.2, lty = 1)
    points(56:65,out[,2], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(56:65,out[,3],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'Cross Validation')
    
    legend("topright", c("AIC", "BIC", "CV"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),
           merge = TRUE, bg = "gray90")
    title(main = "Traditional Methods")
    
    mtext(side=3,outer=T,line=0,paste("p = ", 60,","," n = " ,2000 , ",", "n_l = ","500, 700, 1000, 1500, 2000"),cex=1.3)
    
    #plot(56:65,dat60_20001[,2],xlab = "Size of the nested model", ylab = "ERM2", col = 1, lty = 1, pch = 1, lwd = 1,ylim = c(-1,2))
    #points(56:65,dat60_20001[,3],xlab = "VCD", ylab = "ERM2", col = 2, lty = 2, pch = 2, lwd = 2)
    #points(56:65,dat60_20001[,4], col = 3, lty = 3, pch = 3, lwd = 3)
    #points(30:50,1000*dat_40[,1], col = 4, lty = 4, pch = 4, lwd = 4)
    #legend("topright", c("ERM2", "ERM1", "BIC"), col = c(1,2,3),
    #       text.col = "green4", lty = c(1,2,3), pch = c(1,2,3),lwd = c(1,2,3),
    #       merge = TRUE, bg = "gray90")
    #title(main = paste("ERM1, ERM2, BIC~ VCD,", "p = ", 60,",N = " ,2000 , "\n",",NL = ","500,700,1000,1500,2000"  ))
    #plot(56:65,dat60_2000[,1], xlab = "Size of the nested model", ylab = "Estimated Hhat", col = 4, lty = 4, pch = 4, lwd = 4, ylim = c(55,500))
    #title(main = paste("p = ", 60,",","N = " ,2000 , ",", "\n","NL = ","500, 700, 1000, 1500, 2000"  ))
    
    #cv_tmp = Model_cv(nparam = 60, Size = 2000, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    #cv_tmp = cv_tmp[,56:65]
    #cv = colMeans(cv_tmp)
    #par(mfrow = c(1,2))
    #plot(56:65, cv, lwd = 2, col = gray(0.4), ylab = "Prediction error", ylim = c(0,0.07),
       #  xlab = "Size of the conjectured Model ", main = paste0(n_folds, 
      #                                                          "-fold Cross-Validation", "\n",
     #                                                           "p = ", 60,","," N = " ,2000 ))
    #lines(56:65, cv, lwd = 2, col = "darkred", lty = 2)
    #cv_sd <- apply(cv_tmp, 2, sd)/sqrt(n_folds)
    #errbar(56:65, cv, cv + cv_sd, cv - cv_sd, add = TRUE, col = "steelblue2", pch = 19, 
    #       lwd = 0.5)
    ##################################################################
    # p = 60, N = 700, NL = 100 to 700 by 100
    #################################################################
    dat60_700 = matrix(c(100, 297003, 300983.80, 6562.45,
                         100, 287528.97, 291446.08, 6546.30,
                         100, 287492.75, 291409.61, 6552.77,
                         100, 265973.71, 269741.87, 6504.86,
                         100, 186316.18, 189473.13, 6262.23,
                         100, 184015.42, 187152.93, 6260.08,
                         100, 139781.55, 142518.55, 6074.17,
                         100, 139329.81, 142062.42, 6078.45,
                         100,  62718.39, 64558.12, 5526.23,
                         100,   4337.61, 4835.85, 3661.76,
                         57,    109.96, 184.46, 1058.24,
                         66, 109.51, 189.27, 1059.55,
                         68, 109.56, 190.48, 1065.97,
                         63, 108.96, 186.81, 1069.64,
                         66, 109.06, 188.69, 1076.16,
                         67, 108.83, 188.96, 1080.88,
                         65, 108.75, 187.70, 1087.34,
                         62, 108.64, 185.80, 1093.87,
                         66, 108.76, 188.30, 1100.30,
                         63, 108.62, 186.36, 1106.58,
                         61, 108.26, 184.72, 1111.20), ncol = 4, byrow = TRUE)
    par( mfrow = c(1,1) )
    
    
    
    par(mfrow=c(2,1),oma=c(0,0,3,0),pch=42,font.sub=3,adj=0.5)
    
    par(mar = c(5,5,5,5))
    plot(50:70,dat60_700[,2], pch=16, col = 1, xlab=NA, ylab= "ERM1 and ERM2", cex=1.2, lty = 1)
    points(50:70,dat60_700[,3], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(50:70,dat60_700[,1],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'VC Dimension')
    
    legend("topright", c("ERM1", "ERM2", "VCD"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),lwd = c(1,2,3),
           merge = TRUE, bg = "gray90")
    title(main = "Complexity Methods")
    ####################################
    # Plot for AIC, BIC and CV
    ###################################
    nparam = 60
    n_zero = 10
    Size = 700
    set.seed(10)
    Param_beta = rnorm(nparam, mean = 5, sd = 3)
    Beta = matrix((c(Param_beta, rep(0,n_zero))), ncol = 1)
    data1 = datasim(Beta, mu=0, sigma = 0.4, size = Size,seed=345)
    data1 = data.frame(scale(data1, scale = TRUE, center = TRUE))
    mat_aic_bic = matrix(NA, nrow = ncol(data1)-1, ncol = 2)
    colnames(mat_aic_bic) = c("AIC", "BIC")
    
    for(i in 2:ncol(data1)){
      Model = lm(y ~., data = data1[,1:i])
      mat_aic_bic[i-1,] = c(round(AIC(Model)), round(BIC(Model)))
    }
    
    cv_tmp = Model_cv(nparam = 60, Size = 700, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    out = data.frame(mat_aic_bic, cv = colMeans(cv_tmp))[50:70,]
    
    par(mar = c(5,5,5,5))
    plot(50:70,out[,1], pch=16, col = 1, xlab=NA, ylab= "AIC and BIC", cex=1.2, lty = 1)
    points(50:70,out[,2], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(50:70,out[,3],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'Cross Validation')
    
    legend("topright", c("AIC", "BIC", "CV"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),
           merge = TRUE, bg = "gray90")
    title(main = "Traditional Methods")
    
    mtext(side=3,outer=T,line=0,paste("p = ", 60,","," n = " ,700 , ",", "n_l = ","100 to 700 by 100"),cex=1.3)
    #dat60_7001 = scale(dat60_700)
    #plot(50:70,dat60_7001[,2],xlab = "Size of the nested model", ylab = "ERM2", col = 1, lty = 1, pch = 1, lwd = 1, ylim = c(-1,2))
    #points(50:70,dat60_7001[,3],xlab = "VCD", ylab = "ERM2", col = 2, lty = 2, pch = 2, lwd = 2)
    #points(50:70,dat60_7001[,4], col = 3, lty = 3, pch = 3, lwd = 3)
    #points(30:50,1000*dat_40[,1], col = 4, lty = 4, pch = 4, lwd = 4)
    #legend("topright", c("ERM2", "ERM1", "BIC"), col = c(1,2,3),
    #       text.col = "green4", lty = c(1,2,3), pch = c(1,2,3),lwd = c(1,2,3),
    #       merge = TRUE, bg = "gray90")
    #title(main = paste("ERM1, ERM2, BIC~ VCD,", "p = ", 60,",N = " ,700 , "\n",",NL = ","100,200,300,400,500,600,700"  ))
    #plot(50:70,dat60_700[,1], xlab = "Size of the nested model", ylab = "Estimated Hhat", col = 4, lty = 4, pch = 4, lwd = 4, ylim = c(50,100))
    #title(main = paste("p = ", 60,",","N = " ,700 , ",", "\n","NL = ","100,200,300,400,500,600,700"  ))
    
    #cv_tmp = Model_cv(nparam = 60, Size = 700, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    #cv_tmp = cv_tmp[,50:70]
    #cv = colMeans(cv_tmp)
    #plot(50:70, cv, lwd = 2, col = gray(0.4), ylab = "Prediction error", ylim = c(0,0.35),
    #     xlab = "Size of the conjectured Model ", main = paste0(n_folds, 
    #                                                            "-fold Cross-Validation", "\n",
    #                                                            "p = ", 60,","," N = " ,700))
    #lines(50:70, cv, lwd = 2, col = "darkred", lty = 2)
    #cv_sd <- apply(cv_tmp, 2, sd)/sqrt(n_folds)
    #errbar(50:70, cv, cv + cv_sd, cv - cv_sd, add = TRUE, col = "steelblue2", pch = 19, 
    #       lwd = 0.5)
    ##################################################################
    # p = 60, N = 2000, NL = 100, 200, 300, 400, 500, 600, 700
    ################################################################
    dat60_2000NL = matrix(c(100, 850014.44, 854513.96, 18175.20,
                             100, 823249.21, 827677.43, 18118.81,
                             100, 823161.98, 827589.97, 18126.20,
                             100, 762207.66, 766468.82, 17979.94,
                             100, 532856.53, 536420.52, 17271.60,
                             100, 526203.77, 529745.48, 17254.07,
                             100, 400283.75, 403373.66, 16714.64,
                             100, 399003.71, 402088.69, 16715.83,
                             100, 180012.14, 182086.56, 15131.51,
                             97,  12722.63, 13272.58, 9839.12,
                             56,    307.63, 379.31, 2377.64,
                             70,    307.63, 386.72, 2382.83,
                             65,    307.15, 383.65, 2388.11,
                             68,    307.10, 385.13, 2394.87,
                             63,    306.93, 382.35, 2402.16,
                             66,    307.00, 383.99, 2409.73,
                             68,    306.94, 384.95, 2416.62,
                             65,    306.84, 383.30, 2424.04,
                             67,    306.81, 384.30, 2431.13,
                             64,    306.21, 382.07, 2435.24,
                             66,    306.26, 383.17, 2442.84),ncol = 4, byrow = TRUE)
    par( mfrow = c(1,1) )
    
    
    par(mfrow=c(2,1),oma=c(0,0,3,0),pch=42,font.sub=3,adj=0.5)
    
    par(mar = c(5,5,5,5))
    plot(50:70,dat60_2000NL[,2], pch=16, col = 1, xlab=NA, ylab= "ERM1 and ERM2", cex=1.2, lty = 1)
    points(50:70,dat60_2000NL[,3], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(50:70,dat60_2000NL[,1],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'VC Dimension')
    
    legend("topright", c("ERM1", "ERM2", "VCD"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),lwd = c(1,2,3),
           merge = TRUE, bg = "gray90")
    title(main = "Complexity Methods")
    ####################################
    # Plot for AIC, BIC and CV
    ###################################
    nparam = 60
    n_zero = 10
    Size = 2000
    set.seed(10)
    Param_beta = rnorm(nparam, mean = 5, sd = 3)
    Beta = matrix((c(Param_beta, rep(0,n_zero))), ncol = 1)
    data1 = datasim(Beta, mu=0, sigma = 0.4, size = Size,seed=345)
    data1 = data.frame(scale(data1, scale = TRUE, center = TRUE))
    mat_aic_bic = matrix(NA, nrow = ncol(data1)-1, ncol = 2)
    colnames(mat_aic_bic) = c("AIC", "BIC")
    
    for(i in 2:ncol(data1)){
      Model = lm(y ~., data = data1[,1:i])
      mat_aic_bic[i-1,] = c(round(AIC(Model)), round(BIC(Model)))
    }
    
    cv_tmp = Model_cv(nparam = 60, Size = 2000, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    out = data.frame(mat_aic_bic, cv = colMeans(cv_tmp))[50:70,]
    
    par(mar = c(5,5,5,5))
    plot(50:70,out[,1], pch=16, col = 1, xlab=NA, ylab= "AIC and BIC", cex=1.2, lty = 1)
    points(50:70,out[,2], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(50:70,out[,3],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'Cross Validation')
    
    legend("topright", c("AIC", "BIC", "CV"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),
           merge = TRUE, bg = "gray90")
    title(main = "Traditional Methods")
    
    mtext(side=3,outer=T,line=0,paste("p = ", 60,","," n = " ,2000 , ",", "n_l = ","100 to 700 by 100"),cex=1.3)
    #dat60_2000NL1 = scale(dat60_2000NL)
    #plot(50:70,dat60_2000NL1[,2],xlab = "Size of the nested model", ylab = "ERM2", col = 1, lty = 1, pch = 1, lwd = 1, ylim = c(-1,2))
    #points(50:70,dat60_2000NL1[,3],xlab = "VCD", ylab = "ERM2", col = 2, lty = 2, pch = 2, lwd = 2)
    #points(50:70,dat60_2000NL1[,4], col = 3, lty = 3, pch = 3, lwd = 3)
    #points(30:50,1000*dat_40[,1], col = 4, lty = 4, pch = 4, lwd = 4)
    #legend("topright", c("ERM2", "ERM1", "BIC"), col = c(1,2,3),
      #     text.col = "green4", lty = c(1,2,3), pch = c(1,2,3),lwd = c(1,2,3),
     #      merge = TRUE, bg = "gray90")
    #title(main = paste("ERM1, ERM2, BIC~ VCD,", "p = ", 60,",N = " ,2000 , "\n",",NL = ","100,200,300,400,500,600,700"  ))
    #plot(50:70,dat60_2000NL[,1], xlab = "Size of the nested model", ylab = "Estimated Hhat", col = 4, lty = 4, pch = 4, lwd = 4, ylim = c(50,100))
    #title(main = paste("p = ", 60,",","N = " ,2000 , ",", "\n","NL = ","100,200,300,400,500,600,700"  ))
    
    #cv_tmp = Model_cv(nparam = 60, Size = 2000, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    #cv_tmp = cv_tmp[,50:70]
    #cv = colMeans(cv_tmp)
    #plot(50:70, cv, lwd = 2, col = gray(0.4), ylab = "Prediction error", ylim = c(0,0.35),
    #     xlab = "Size of the conjectured Model ", main = paste0(n_folds, 
     #                                                           "-fold Cross-Validation", "\n",
      #                                                          "p = ", 60,","," N = " ,2000 ))
    #lines(50:70, cv, lwd = 2, col = "darkred", lty = 2)
    #cv_sd <- apply(cv_tmp, 2, sd)/sqrt(n_folds)
    #errbar(50:70, cv, cv + cv_sd, cv - cv_sd, add = TRUE, col = "steelblue2", pch = 19, 
    #       lwd = 0.5)
    ##################################################################
    # p = 60, N = 2000, NL = 75 to 600 by 75
    ##################################################################
    dat60_2000NL2 = matrix(c( 75, 865665.45, 869722.05, 18211.69,
                              75,836823.97, 840812.51, 18151.53,
                              75,804788.23, 808699.78, 18081.06,
                              75,639641.54, 643129.29, 17629.31,
                              75,378379.55, 381063.24, 16586.88,
                              75,293689.27, 296054.24, 16087.73,
                              75,185570.32, 187451.28, 15177.13,
                              75,175086.88, 176914.08, 15068.43,
                              75,88255.29, 89554.04, 13705.88,
                              75,69458.37, 70611.13, 13234.45,
                             50,307.46, 375.63, 2377.64,
                             57,307.2945, 379.4994, 2382.83,
                             55,306.8884, 377.9203, 2388.11,
                             54,306.7329, 377.1755, 2394.87,
                             52,306.6282, 375.8976, 2402.16,
                             52,306.6232, 375.8919, 2409.73,
                             50,306.4577, 374.5224, 2416.62,
                             48,306.3705, 373.2136, 2424.04,
                             46,306.2312, 371.8198, 2431.13,
                             44,305.6412, 369.9006, 2435.24,
                             77,306.4710, 387.8766, 2442.84), ncol = 4, byrow = TRUE)
    par( mfrow = c(1,1) )
    par(mfrow=c(2,1),oma=c(0,0,3,0),pch=42,font.sub=3,adj=0.5)
    
    par(mar = c(5,5,5,5))
    plot(50:70,dat60_2000NL2[,2], pch=16, col = 1, xlab=NA, ylab= "ERM1 and ERM2", cex=1.2, lty = 1)
    points(50:70,dat60_2000NL2[,3], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(50:70,dat60_2000NL2[,1],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'VC Dimension')
    
    legend("topright", c("ERM1", "ERM2", "VCD"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),lwd = c(1,2,3),
           merge = TRUE, bg = "gray90")
    title(main = "Complexity Methods")
    ####################################
    # Plot for AIC, BIC and CV
    ###################################
    nparam = 60
    n_zero = 10
    Size = 2000
    set.seed(10)
    Param_beta = rnorm(nparam, mean = 5, sd = 3)
    Beta = matrix((c(Param_beta, rep(0,n_zero))), ncol = 1)
    data1 = datasim(Beta, mu=0, sigma = 0.4, size = Size,seed=345)
    data1 = data.frame(scale(data1, scale = TRUE, center = TRUE))
    mat_aic_bic = matrix(NA, nrow = ncol(data1)-1, ncol = 2)
    colnames(mat_aic_bic) = c("AIC", "BIC")
    
    for(i in 2:ncol(data1)){
      Model = lm(y ~., data = data1[,1:i])
      mat_aic_bic[i-1,] = c(round(AIC(Model)), round(BIC(Model)))
    }
    
    cv_tmp = Model_cv(nparam = 60, Size = 2000, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    out = data.frame(mat_aic_bic, cv = colMeans(cv_tmp))[50:70,]
    
    par(mar = c(5,5,5,5))
    plot(50:70,out[,1], pch=16, col = 1, xlab=NA, ylab= "AIC and BIC", cex=1.2, lty = 1)
    points(50:70,out[,2], pch=10, col = 2, xlab=NA, ylab= NA, cex=1, lty = 2)
    par(new = T)
    plot(50:70,out[,3],xlab = "Conjectured Model Size", axes = F, ylab = NA, col = 3, lty = 3, pch = 2)
    axis(side = 4)
    mtext(side = 4, line = 3, 'Cross Validation')
    
    legend("topright", c("AIC", "BIC", "CV"), col = c(1,2,3),
           text.col = "green4", lty = c(1,2,3), pch = c(16,10,2),
           merge = TRUE, bg = "gray90")
    title(main = "Traditional Methods")
    
    mtext(side=3,outer=T,line=0,paste("p = ", 60,","," n = " ,2000 , ",", "n_l = ","75 to 600 by 75"),cex=1.3)
    #dat60_2000NL21 = scale(dat60_2000NL2)
    #plot(50:70,dat60_2000NL21[,3],xlab = "Size of the nested model", ylab = "ERM2", col = 1, lty = 1, pch = 1, lwd = 1, ylim = c(-1,2.5))
    #points(50:70,dat60_2000NL21[,2],xlab = "VCD", ylab = "ERM2", col = 2, lty = 2, pch = 2, lwd = 2)
    #points(50:70,dat60_2000NL21[,4], col = 3, lty = 3, pch = 3, lwd = 3)
    #points(30:50,1000*dat_40[,1], col = 4, lty = 4, pch = 4, lwd = 4)
    #legend("topright", c("ERM2", "ERM1", "BIC"), col = c(1,2,3),
    #       text.col = "green4", lty = c(1,2,3), pch = c(1,2,3),lwd = c(1,2,3),
    #       merge = TRUE, bg = "gray90")
    #title(main = paste("ERM1, ERM2, BIC~ VCD,", "p = ", 60,",N = " ,2000 , "\n",",NL = ","75 to 600 by 75"  ))
    #plot(50:70,dat60_2000NL2[,1], xlab = "Size of the nested model", ylab = "Estimated Hhat", col = 4, lty = 4, pch = 4, lwd = 4)
    #title(main = paste("p = ", 60,",","N = " ,2000 , ",", "\n","NL = ","75 to 600 by 75"  ))
    
    #cv_tmp = Model_cv(nparam = 60, Size = 2000, n_folds <- 10, sigma = 0.4, n_zero = 10)
    
    #cv_tmp = cv_tmp[,50:70]
    #cv = colMeans(cv_tmp)
    #plot(50:70, cv, lwd = 2, col = gray(0.4), ylab = "Prediction error", ylim = c(0,0.35),
      #   xlab = "Size of the conjectured Model ", main = paste0(n_folds, 
     #                                                           "-fold Cross-Validation", "\n",
       #                                                         "p = ", 60, ", ", " N = ", 2000))
    #lines(50:70, cv, lwd = 2, col = "darkred", lty = 2)
    #cv_sd <- apply(cv_tmp, 2, sd)/sqrt(n_folds)
    #errbar(50:70, cv, cv + cv_sd, cv - cv_sd, add = TRUE, col = "steelblue2", pch = 19, 
     #      lwd = 0.5)
    ####################################################################################################################################
    
    # Simulation using Skew distribution for the error and uniform (-1,1) to X
    ##################################################################################################################################
    library(skewt)
    
    datasim_skew = function(nparam, gamma, dof, size, n_zero, seed){
      set.seed(seed)
      Param_beta = rnorm(nparam, mean = 5, sd = 3)
      Beta = matrix((c(Param_beta, rep(0,n_zero))), ncol = 1)
      e.sim<-matrix(data = rskt(n = size, df = dof, gamma = gamma), nrow = size, ncol = 1) 
      x.sim = matrix(runif(n = size*dim(Beta)[1], min = -1, max = 1), nrow = size, ncol = dim(Beta)[1])
      y.sim = x.sim %*% Beta + e.sim
      
      output = data.frame(y=y.sim, x.sim)
    }
    
    dataset = datasim_skew(nparam = 60, gamma = 50, dof = 10, size = 2000, n_zero = 10, seed = 345)
    
    Chxi_skew = function(data){
      NL = seq(from = 75, to = 600, by = 75)
      B = 50
      m = 10
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
            #G1 = as.data.frame(scale(SampleData[Index,], center = TRUE, scale = TRUE))
            #G2 = as.data.frame(scale(SampleData[-Index,], center = TRUE, scale = TRUE))
            
            # Now lets fit a model using the modify dataset.
            Model1 = lm(y ~  ., x = TRUE, y = TRUE, data = G1)
            Model2 = lm(y~  ., x = TRUE, y = TRUE, data = G2)
            #PredSqrt = (Model$residuals)^2
            
            FirstHalfN = matrix(NA, nrow =length(Index), ncol =m)
            FirstHalfQ = matrix(NA, nrow =length(Index), ncol =m)
            colnames(FirstHalfN) =paste("N",1:m,sep="")
            colnames(FirstHalfQ) = paste("Q",1:m,sep="")
            
            X1 = data.frame(Model1$x)
            Y1 = Model1$y
            X2 = data.frame(Model2$x)
            Y2 = Model2$y
            #Pred1 = PredSqrt[Index]
            qstarj = matrix(NA, nrow = 1, ncol = m)
            Pred1 = (predict(Model1, X2)- Y2)^2
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
          MatChixi[bprim,i] = sum (apply(Matxi, 2, max, na.rm = T))
          #MatChixi[bprim,i] = sum (apply(Matxi, 2, mean, na.rm = T))
          
          
          bprim = bprim + 1
        }
        
      }
      
      MeanChi = apply(MatChixi, MARGIN =2, FUN = mean)
    }
    df = 2:ncol(dataset)
    est  <- apply(t(df), 2, function(degf) Chxi_skew(dataset[,1:degf]))
    