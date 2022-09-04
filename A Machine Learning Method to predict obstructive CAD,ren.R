

#############
rm(list=ls())   
options(stringsAsFactors = F)
library("pROC")
library(caret)
library(randomForest)
library(rminer) 
library(glmnet)
library(neuralnet)
library(e1071)
library(dplyr)

#######################################
## L# Data preparation for models ##
#######################################
data<-read.csv(...)
set.seed(516) 
smp.size = floor(0.7*nrow(data)) 
train.ind = sample(seq_len(nrow(data)), smp.size)
traindata<-data[train.ind, ] # 70%
testdata<-data[-train.ind, ]


# Set Training Control
myTrainingControl <- trainControl(method = "cv", 
                                  number = 5,
                                  savePredictions = TRUE, 
                                  classProbs = TRUE, 
                                  verboseIter = FALSE)# Train RF
###################
## Random forest ##
###################

fit_RF <- caret::train(CHD ~ .,   
                       data = traindata, 
                       method = "ranger", 
                       tuneLength = 3,     
                       importance = "permutation",
                       trControl = myTrainingControl)
print(fit_RF)
# Prediction in test set
pred_RF_prob <- predict(fit_RF, testdata) 



###################
## CACS ##
###################

fit_glm <- caret::train(CHD ~CACS,   
                        data = traindata, 
                        method = "glm",
                        family = "binomial",
                        trControl = myTrainingControl)
print(fit_glm)
# Prediction in test set
pred_glm_prob1 <- predict(fit_glm, testdata, type="prob")



###################
## clinical ##
###################

fit_glm2 <- caret::train(CHD~sex+Age+Diabetes+Hypertention+Current_smoker+Total_cholesterol
                         +LDL_C+HDL_C+CACS,   
                         data = traindata, 
                         method = "glm",
                         family = "binomial", 
                         trControl = myTrainingControl)
print(fit_glm2)
# Prediction in test set
pred_glm_prob2 <- predict(fit_glm2, testdata, type="prob")


###################
## ROC ##
###################
df<-cbind(testdata$CHD,pred_RF_prob$yes,pred_glm_prob1$yes,pred_glm_prob2$yes)

mycol <- c('slateblue','seagreen3','dodgerblue','firebrick1','lightgoldenrod','magenta','orange2')
x1<-plot.roc(df[,1],df[,2],
             smooth=F,
             lwd=2,
             ylim=c(0,1),
             xlim=c(1,0),
             legacy.axes=T,
             main='',
             col=mycol[2],print.thres=TRUE)
auc(x1) 
x2<-plot.roc(df[,1],df[,3],
             smooth=F,
             add=T,
             lwd=2,
             ylim=c(0,1),
             xlim=c(1,0),
             legacy.axes=T,
             main='',
             col=mycol[3],print.thres=TRUE)
auc(x2)
x3<-plot.roc(df[,1],df[,4],
             smooth=F,
             add=T,
             lwd=2,
             ylim=c(0,1),
             xlim=c(1,0),
             legacy.axes=T,
             main='',
             col=mycol[4],print.thres=TRUE)
auc(x3)
legend.name <- c(paste('RF model','AUC',0.841,sep=' '),paste('LR:CACS model','AUC',0.746,sep=' '),paste('LR:Clinical model','AUC',0.810,sep=' ')
)
legend('bottomright', 
       legend=legend.name,
       col = mycol[2:4],
       lwd = 2,
       bty='n')
###################
## SE, SP, PPV,NPV ##
###################


library(OptimalCutpoints)
df<-data.frame(df)
optimal.cutpoint.Youden <- optimal.cutpoints(X = "pred_RF_prob", status = "V1", tag.healthy = 0,  
                                             methods = "Youden", data = df, pop.prev = NULL,  
                                             control = control.cutpoints(), ci.fit = TRUE, conf.level = 0.95, trace = FALSE)
optimal.cutpoint.Youden <- optimal.cutpoints(X = "pred_glm_prob1", status = "V1", tag.healthy = 0,  
                                             methods = "Youden", data = df, pop.prev = NULL,  
                                             control = control.cutpoints(), ci.fit = TRUE, conf.level = 0.95, trace = FALSE)
optimal.cutpoint.Youden <- optimal.cutpoints(X = "pred_glm_prob2", status = "V1", tag.healthy = 0,  
                                             methods = "Youden", data = df, pop.prev = NULL,  
                                             control = control.cutpoints(), ci.fit = TRUE, conf.level = 0.95, trace = FALSE)
summary(optimal.cutpoint.Youden)
plot(optimal.cutpoint.Youden)


###################
## NRI,IDI ##
###################

library('PredictABEL')
library('nricens')

pred1<-df$pred_glm_prob2
pred2<-df$pred_RF_prob
pred3<-df$pred_glm_prob1
lables<-df$V1
####
NRIb = nribin(event = df$V1,  p.std = pred1, p.new = pred2, updown ='category', cut=c(0.2,0.4),niter = 1000) 
NRIb = nribin(event = df$V1,  p.std = pred1, p.new = pred2, updown ='diff', cut=0,niter = 1000) 
NRIb = nribin(event = df$V1,  p.std = pred3, p.new = pred2, updown ='diff', cut=0,niter = 1000) 
NRIb = nribin(event = df$V1,  p.std = pred3, p.new = pred1, updown ='diff', cut=0,niter = 1000) 
reclassification(data=df, cOutcome = 1, predrisk1 = pred3, predrisk2 = pred1, cutoff = c(0,0.2,0.4,1))


###################
## DCA ##
###################

library(rmda)
model_1<- decision_curve(V1~pred_RF_prob,data=df,
                         family = binomial(link = "logit"),
                         thresholds = seq(0,1,by=0.01),
                         confidence.intervals = 0.95,
                         study.design = "case-control",
                         population.prevalence = 0.3) 
model_2<- decision_curve(V1~pred_glm_prob1,data=df,
                         family = binomial(link = "logit"),
                         thresholds = seq(0,1,by=0.01),
                         confidence.intervals = 0.95,
                         study.design = "case-control",
                         population.prevalence = 0.3) 

model_3<- decision_curve(V1~pred_glm_prob2,data=df,
                         family = binomial(link = "logit"),
                         thresholds = seq(0,1,by=0.01),
                         confidence.intervals = 0.95,
                         study.design = "case-control",
                         population.prevalence = 0.3) 

plot_decision_curve(model_1,curve.names = c("RF model"),xlim=c(0,0.8),
                    cost.benefit.axis = FALSE,col =c("red"),
                    confidence.intervals = FALSE,
                    standardize = FALSE)    
plot_decision_curve(model_2,curve.names = c("CACS model"),xlim=c(0,0.8),
                    cost.benefit.axis = FALSE,col=c("blue"),
                    confidence.intervals = FALSE,
                    standardize = FALSE) 
plot_decision_curve(model_3,curve.names = c("Clinical model"),xlim=c(0,0.8),
                    cost.benefit.axis = FALSE,col=c("green"),
                    confidence.intervals = FALSE,
                    standardize = FALSE) 

model_all<-list(model_1,model_2,model_3)

plot_decision_curve(model_all,curve.names = c("RF model","CACS model","Clinical model"),xlim=c(0,0.8),
                    cost.benefit.axis = FALSE,col=c("red","blue","green"),
                    confidence.intervals = FALSE,
                    standardize = FALSE)
