library(lattice)
library(ggplot2)
library(caret)
library(e1071)
library(nnet)
library(pROC)
library(psych)
library(car)
library(e1071)
library(MVN)
library(tidyverse)
library(crossval)

# remove first variable
glass <- glass[,-1]

count(glass, vars = V11)






#################### Exploratory Analysis #####################






glassex <- glass
glassex$V11 <- as.factor(glassex$V11)


## scatter plot
library(car)
scatterplotMatrix(glass[,1:9],use = "complete.obs",  diagonal="histogram", smooth=FALSE)

cor(glass[,1:9])


## boxplot
par(mfrow=c(1,3))
boxplot(glass$V2 ~ glass$V11, ylab="V2")
boxplot(glass$V3~ glass$V11, ylab="V3")
boxplot(glass$V4 ~ glass$V11, ylab="V4")

par(mfrow=c(1,3))
boxplot(glass$V5 ~ glass$V11, ylab="V5")
boxplot(glass$V6 ~ glass$V11, ylab="V6")
boxplot(glass$V7~ glass$V11, ylab="V7")

par(mfrow=c(1,3))
boxplot(glass$V8~ glass$V11, ylab="V8")
boxplot(glass$V9 ~ glass$V11, ylab="V9")
boxplot(glass$V10 ~ glass$V11, ylab="V10")


## normality test
library(MVN)
result <- mvn(data = glass[ ,1:9], mvnTest = "mardia") 
result$multivariateNormality

hist(glass[ ,1:9])

## normality test for some classes

result <- mvn(data = glass[which(glass$V11 == 1) ,1:9], mvnTest = "mardia") 
result$multivariateNormality

result <- mvn(data = glass[which(glass$V11 == 2) ,1:9], mvnTest = "mardia") 
result$multivariateNormality

result <- mvn(data = glass[which(glass$V11 == 3) ,1:9], mvnTest = "mardia",univariateTest = "AD", univariatePlot = "histogram",
              multivariatePlot = "qq", multivariateOutlierMethod = "adj") 
result$multivariateNormality

result <- mvn(data = glass[which(glass$V11 == 5) ,1:9], mvnTest="energy")
result$multivariateNormality

## power

# powera<- boxCox( 
#          glass$V11~cbind(glass$V2,glass$V3,glass$V4,glass$V5,glass$V6,glass$V7,
#              glass$V8,glass$V9,glass$V10)
#        , family="yjPower", plotit = TRUE)
# pl<- cbind(powera$x,powera$y)
# 
# detans <- yjPower(cbind(glass$V2,glass$V3,glass$V4,glass$V5,glass$V6,glass$V7,
#                         glass$V8,glass$V9,glass$V10),-0.1,jacobian.adjusted=FALSE)
# 
# 
# glasstrans <-detans
# 
# colnames(glasstrans) <- c("V2","V3","V4","V5","V6","V7","V8","V9","V10")
# 
# result <- mvn(data = glasstrans[ ,1:9], mvnTest = "hz") 
# result$multivariateNormality


## Manva

glmn <-manova(cbind(glass$V2,glass$V3,glass$V4,glass$V5,glass$V6,glass$V7,
                    glass$V8,glass$V9,glass$V10)
              ~ 
                glass$V11)

summary(glmn, test="Wilks")



########### SVM ##################



glassexst <- glassex

glassexst[,1:9] <-  scale(glassex[,1:9],center=T,scale=T)

## split dataset into train and test

ind <- sample(2, nrow(glassexst), replace=TRUE, prob=c(0.8, 0.2))
traindata <-glassexst[ind==1,]#train data
testdata <- glassexst[ind==2,]#test data

traindata<-transform(traindata,V11=as.factor(V11))
testdata<-transform(testdata,V11=as.factor(V11))

svm.model<-svm(traindata$V11~., traindata[,-10],cross =5)
summary(svm.model)

confusion.sv=table(traindata$V11,predict(svm.model,traindata,type="V11"))
accuracy.sv=sum(diag(confusion.sv))/sum(confusion.sv)
confusion.sv
accuracy.sv
1- accuracy.sv




############ logistic #############

Y1 = as.factor(glass[,10])
X1 = as.matrix(glass[,c(-10)])




# train_sub = sample(nrow(glass),8/10*nrow(glass))
# train_data = glass[train_sub,]
# test_data =glass[-train_sub,]


## predicting function 
predfun.lg = function(train.x, train.y, test.x, test.y)

{

train_data <- as.data.frame(cbind(train.x,train.y))
colnames(train_data)[10] <-"V11"

test_data <- cbind(test.x,test.y)
colnames(test_data)[10] <-"V11"

# train_data$class2<-relevel(as.factor(train.y),ref = "1")\
mult.model<-multinom(V11~V2+V3+V4+V5+V6+V7+V8+V9,maxit=1000,data = train_data)

# #train_data[,1]+train_data[,2]+
#   train_data[,3]+train_data[,4]+
#   train_data[,5]+train_data[,6]+train_data[,7]+
#   train_data[,8]+train_data[,9]

# summary(mult.model)


# z <- summary(mult.model)$coefficients/summary(mult.model)$standard.errors
# p <- (1 - pnorm(abs(z), 0, 1))*2
# p

pre_logistic<-predict(mult.model,newdata = test.x) 

# table(test_data$class,pre_logistic)

out <-  table(test.y,pre_logistic)
# confusionMatrix(factor(pre_logistic),factor(test.y),negative=negative)

return(out)

}


# cross validation
cv.out.lg = crossval(predfun.lg, X1, Y1, K=5, B=1)

# calculate accruacy 

x <- rep(0,180)
dim(x) <- c(6,6,5)

mat <- as.matrix(cv.out.lg$stat.cv[,1:6])

rownames(mat)=NULL

for (i in 1:5 ) {
  x[,,i] <- mat[(6*i-5):(6*i),]
}
sum(x,3)

atable <- apply(x,c(1,2),sum)


# accuracy
acc.lg <- tr(atable)/sum(atable)
error.lg <-1-acc.lg
error.lg 


############  random forest ##############

library(randomForest)

rf.glass <- glass
rf.glass$V11 <- as.factor(rf.glass$V11)

# select parameter mtry
training <- rf.glass

rf.e <- c(0,0,0,0,0,0,0,0)

for (i in 2:9) {

mt = i

rf_classifier = randomForest(V11 ~ ., data=training, ntree=1000, mtry = mt, importance=TRUE)

rf.e[i-1] <- 1-tr(rf_classifier$confusion[,-7])/sum(rf_classifier$confusion[,-7])

}

## draw the mtry plot
plot(rf.e)

m_rf_im<-train(V11~.,data=glass,method='rf',metric='Kappa',
               tuneGrid=data.frame(.mtry=c(2,4,6,8)))

rf_classifier = randomForest(V11 ~ ., data=training, ntree=1000, mtry = 3, importance=TRUE)


# importance plot
importance(rf_classifier,type=1) 

importance(rf_classifier,type=2)

varImpPlot(rf_classifier)   



################## LDA ###########################
library(crossval)
c1 = glass[,10]
X1 = as.matrix(glass[,c(-10)])
Y1 = as.factor(c1)
levels(Y1)


predfun.lda = function(train.x, train.y, test.x, test.y, negative)
{
  require("MASS") # for lda function
  lda.fit = lda(train.x, grouping=train.y, prior=c(1,1,1,1,1,1)/6)
  ynew = predict(lda.fit, test.x)$class
  # count TP, FP etc.
  out =  table(test.y,ynew)
  return( out )
}


## 5_fold Cross validation
cv.out.l = crossval(predfun.lda, X1, Y1, K=5, B=1, negative="1")



## calculate error rate
x <- rep(0,180)
dim(x) <- c(6,6,5)

mat <- as.matrix(cv.out.l$stat.cv[,1:6])

rownames(mat)=NULL

for (i in 1:5 ) {
  x[,,i] <- mat[(6*i-5):(6*i),]
}

atable <- apply(x,c(1,2),sum)


# accuracy
acc.l <- tr(atable)/sum(atable)
error.l <- 1-acc.l
error.l
