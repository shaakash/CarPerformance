#---------------------------------------------------------------------------------------------------------
#R Code for Chapter 17 of:
#
#Field, A. P., Miles, J. N. V., & Field, Z. C. (2012). Discovering Statistics Using R: and Sex and Drugs and Rock 'N' Roll. #London Sage
#
#(c) 2011 Andy P. Field, Jeremy N. V. Miles & Zoe C. Field
#-----------------------------------------------------------------------------------------------------------




#----Set the working directory------
setwd("F:/Sridhar/PGPBA/Courses/Advanced Statistics Course/Pune/Advanced Statistics Pune/R Programs")


#------And then load these packages, along with the boot package.-----

library(corpcor)
library(GPArotation)
library(psych)
library(ggplot2)
library(ggfortify)
library(nFactors)
library(plyr)
library(reshape)

#********************* RAQ Example ********************

#load data and made labels for each car
carData1<-read.table("mbacar.txt", header = TRUE)
ad1<-factor(carData1$Car, levels = c(1:10), labels=c("BMW328i","Ford Explorer", "Infiniti J30", "Jeep Grand Cherokee", "Lexus ES300", "Chrysler Town & Country", "Mercedes C280", "Saab 9000", "Porche Boxster", "Volvo V90"))
carData1<-cbind(carData1,ad1)
ad2<-c(1:10)

ad3<-factor(levels=c(1:10), labels=c("BMW328i","Ford Explorer", "Infiniti J30", "Jeep Grand Cherokee", "Lexus ES300", "Chrysler Town & Country", "Mercedes C280", "Saab 9000", "Porche Boxster", "Volvo V90"))

ad4=cbind(ad2,ad3)


carData<-carData1[3:18]
carData2<-carData1[1:2]
#create a correlation matrix
carMatrix<-cor(carData)
round(carMatrix, 2)
#break down the matrix to make it easier to put in the book
round(carMatrix[,1:8], 2)
round(carMatrix[,9:16], 2)

#Bartlett's test

cortest.bartlett(carMatrix, n = 303)

#KMO test


# KMO Kaiser-Meyer-Olkin Measure of Sampling Adequacy
# Function by G. Jay Kerns, Ph.D., Youngstown State University (http://tolstoy.newcastle.edu.au/R/e2/help/07/08/22816.html)

kmo = function( data ){
  library(MASS) 
  X <- cor(as.matrix(data)) 
  iX <- ginv(X) 
  S2 <- diag(diag((iX^-1)))
  AIS <- S2%*%iX%*%S2                      # anti-image covariance matrix
  IS <- X+AIS-2*S2                         # image covariance matrix
  Dai <- sqrt(diag(diag(AIS)))
  IR <- ginv(Dai)%*%IS%*%ginv(Dai)         # image correlation matrix
  AIR <- ginv(Dai)%*%AIS%*%ginv(Dai)       # anti-image correlation matrix
  a <- apply((AIR - diag(diag(AIR)))^2, 2, sum)
  AA <- sum(a) 
  b <- apply((X - diag(nrow(X)))^2, 2, sum)
  BB <- sum(b)
  MSA <- b/(b+a)                        # indiv. measures of sampling adequacy
  AIR <- AIR-diag(nrow(AIR))+diag(MSA)  # Examine the anti-image of the correlation matrix. That is the  negative of the partial correlations, partialling out all other variables.
  kmo <- BB/(AA+BB)                     # overall KMO statistic
  # Reporting the conclusion 
  if (kmo >= 0.00 && kmo < 0.50){test <- 'The KMO test yields a degree of common variance unacceptable for FA.'} 
  else if (kmo >= 0.50 && kmo < 0.60){test <- 'The KMO test yields a degree of common variance miserable.'} 
  else if (kmo >= 0.60 && kmo < 0.70){test <- 'The KMO test yields a degree of common variance mediocre.'} 
  else if (kmo >= 0.70 && kmo < 0.80){test <- 'The KMO test yields a degree of common variance middling.' } 
  else if (kmo >= 0.80 && kmo < 0.90){test <- 'The KMO test yields a degree of common variance meritorious.' }
  else { test <- 'The KMO test yields a degree of common variance marvelous.' }
  
  ans <- list( overall = kmo,
               report = test,
               individual = MSA,
               AIS = AIS,
               AIR = AIR )
  return(ans)
} 

#To use this function:
kmo(carData)

#Determinent (execute one of these):
det(carMatrix)
det(cor(carData))

#PCA

#pcModel<-principal(dataframe/R-matrix, nfactors = number of factors, rotate = "method of rotation", scores = TRUE)

#On raw data

pc1 <-  principal(carData, nfactors = 16, rotate = "none")
pc1 <-  principal(carData, nfactors = length(carData), rotate = "none")
plot(pc1$values, type = "b") 

pc2 <-  principal(carData, nfactors = 4, rotate = "none")
pc2$loadings

#Explore residuals
factor.model(pc2$loadings)
reproduced<-round(factor.model(pc2$loadings), 3) #format for book
reproduced[,1:9]  #format for book

factor.residuals(carMatrix, pc2$loadings) 
resids<-round(factor.residuals(carMatrix, pc2$loadings), 3) #format for book
resids[,1:9] #format for book

#####A measure of the fit of the model is therefore the sum of the squared residuals divided by the sum of the squared correlations. 
#####As this is considered a measure of fit and sometimes people like measures of fit to go from 0 to 1, we subtract the value from 1.
pc2$fit.off

residuals<-factor.residuals(carMatrix, pc2$loadings)
residuals<-as.matrix(residuals[upper.tri(residuals)])
large.resid<-abs(residuals) > 0.05
sum(large.resid)
sum(large.resid)/nrow(residuals)
sqrt(mean(residuals^2))
hist(residuals)


residual.stats<-function(matrix){
  residuals<-as.matrix(matrix[upper.tri(matrix)])
  large.resid<-abs(residuals) > 0.05
  numberLargeResids<-sum(large.resid)
  propLargeResid<-numberLargeResids/nrow(residuals)
  rmsr<-sqrt(mean(residuals^2))
  
  cat("Root means squared residual = ", rmsr, "\n")
  cat("Number of absolute residuals > 0.05 = ", numberLargeResids, "\n")
  cat("Proportion of absolute residuals > 0.05 = ", propLargeResid, "\n")
  hist(residuals)
}

resids <- factor.residuals(carMatrix, pc2$loadings )
residual.stats(resids)
residual.stats(factor.residuals(carMatrix, pc2$loadings))


#Factor rotation

pc3 <-  principal(carData, nfactors = 4, rotate = "varimax")
print.psych(pc3, cut = 0.3, sort = TRUE)

pc4 <- principal(carData, nfactors = 5, rotate = "oblimin")
print.psych(pc4, cut = 0.3, sort = TRUE)
pc4$loadings%*%pc4$Phi


factor.structure <- function(fa, cut = 0.2, decimals = 2){
  structure.matrix <- fa.sort(fa$loadings %*% fa$Phi)
  structure.matrix <- data.frame(ifelse(abs(structure.matrix) < cut, "", round(structure.matrix, decimals)))
  return(structure.matrix)
}

factor.structure(pc4, cut = 0.3)

###Another method to view factor analysis
pcal<- fa(carData,5,n.obs=303)
pcal
fa.diagram(pcal)
########################
##Scree Plot
########################

ev <- eigen(cor(carData)) # get eigenvalues
ap <- parallel(subject=nrow(carData),var=ncol(carData),
               rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)


#Factor scores

pc5 <- principal(carData, nfactors = 4, rotate = "oblimin", scores = TRUE)

#pc5$scores
head(pc5$scores, 10)
summary(pc5)
loadings(pc5)

factor.structure(pc5, cut = 0.3)
plot(pc5$values,type="lines") # scree plot 
#pc5$scores # the principal components
biplot(pc5)
loadss <- data.frame(pc5$loadings[,1:2])
plot(loadss[,1],loadss[,2],type="b")
text(loadss,labels=names(ad1),cex=.7) # add variable names

### Labeling factors
funandambition1<-cbind(carData1$stylish,carData1$fun,carData1$exciting,carData1$status,carData1$performance,carData1$powerful,carData1$sports,carData1$luxurious)
reliability1<-cbind(carData1$safe,carData1$dependable,carData1$comfortable)
adventurous1<-cbind(carData1$rugged, carData1$outdoorsy)
utility1<-cbind(carData1$practical,carData1$versatile,carData1$family)



##Diagram
#names(pc5$scores)<-c("funandambition","reliability","adventurous","utility")
carData1 <- cbind(carData1, pc5$scores)
rename(carData1, c(TC1="funandambition",TC2="reliability",TC3="adventurous", TC4="utility"))



autoplot(prcomp(carData), data = carData1, colour = 'blue', loadings=TRUE, loadings.colour = 'green',)

ef<-tapply(carData1$TC1,list(carData1$ad1), mean)
fd<-tapply(carData1$TC2,list(carData1$ad1), mean)

meanvector<-cbind(ad2,ad1,ef,fd)
scatter<-ggplot(meanvector, aes(ef, fd))
scatter + geom_point() + labs(x = "Fun and Ambition", y = "Reliability")
scatter+geom_text(aes(labels(ad1)))

interaction.plot(ef, fd, ad2,data=meanvector)
 
###############Discriminant Analysis using R

set.seed(124)
c<-ifelse(runif(303)<.5,1,2)
q<-t(c)
z<-t(q)
carData1<-cbind(carData1,z)
