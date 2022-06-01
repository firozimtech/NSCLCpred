# This code is use for developing a LASSO model for predicting NSCLC
# Dr. Firoz Ahmed, University of Jeddah

library(glmnet)
library(tidyverse) # for data manupulation
library(plotmo) # for plot_glmnet
library(ROCR) # for ROC curve
library(caret) # for confusion Matrix

set.seed(123) # seed for reproducibility
# https://www.statology.org/lasso-regression-in-r/
# https://bookdown.org/tpinto_home/Regularisation/lasso-regression.html
# http://www.sthda.com/english/articles/36-classification-methods-essentials/149-penalized-logistic-regression-essentials-in-r-ridge-lasso-and-elastic-net/#r-functions

# Set the working directory path
# setwd("/Users/firozahmed/Documents/LASSO")

# ----------------Step 1: Upload the the expression Data----------------
# Expression Data: 1st column is the class name (eg Cancer, Normal).
# Rest columns are gene name and their expression values. 
# Each row contains expression values of different genes either in each sample (cancer or normal).
#library(readr)
#library(tidyverse)
#library(vroom)


# read the gene expression data having columns as sample names, rows as gene names
lungExpData<-read.table(file = "GeneExp.txt", header	 = TRUE, sep = "\t", row.names = 1, check.names=FALSE ) 
# check.names=FALSE (avoid converting "-" in column name to "." eg "TCGA-12-34" TCGA.12.34")
View(head(lungExpData))
str(lungExpData)

# selected gene names with log2FC and highest node degree
Gene_DEG<-read.table("DEGs_degree.txt", header=T, sep="\t", check.names=FALSE)

# take only important genes (40 DEGs) expression values from the lungExpData
lungExpData_2<-lungExpData[Gene_DEG$Gene,]
View(head(lungExpData_2))

dim(lungExpData_2)
# read the sample name and its associated sample type (Cancer vs Normal)
SampleName<-read.table("Sample_Name_3.csv", header=T, sep=",", row.names = 1, check.names=FALSE)

# Transpose the gene expression table
lungExpData_3 <-t(lungExpData_2)

# Merge the two data column by column based upon common row names
lungExpData_4<-merge(as.data.frame(SampleName), as.data.frame(lungExpData_3), by='row.names')

#row names is in the First column having sample names in lungExpData_4
dim(SampleName)
dim(lungExpData_4)

#typeof(lungExpData1)
View(head(lungExpData_4))

rownames((SampleName))
# Adding new column (Class) after column "sample_type"  
# If "sample_type" is "Cancer" then "1"
# If "sample_type" is "Normal" then "0"

library(tibble)
library(dplyr)

lungExpData_5<-lungExpData_4 %>%  mutate(Class = case_when(grepl("Cancer", sample_type) ~ 1, 
                                              grepl("Normal", sample_type) ~ 0), .after="sample_type" ) 
   
View(head(lungExpData_5))
rownames(lungExpData_4)
NSCSdata<-lungExpData_5
dim(NSCSdata)

# number of rows(observation) in gene matrix data
n<-nrow(NSCSdata)
n
set.seed(198)  # Set seed for reproducibility
# Split data into train (80%: 4/5) and test (20%: 1/5) sets
train_rows <- sample(1:n, .80*n)
train <- NSCSdata[train_rows, ] # matrix of training data
test <- NSCSdata[-train_rows, ] # matrix of testing data

dim(train)
dim(test)
library(stringr)
# number of different data in Train data
sum(str_count(train$sample3,"Normal Tissue")) #GTEx
sum(str_count(train$sample3,"Primary Tumor")) #Cancer
sum(str_count(train$sample3,"Recurrent Tumor")) #Cancer
sum(str_count(train$sample3,"Solid Tissue Normal")) #Normal

# number of different data in Test data
sum(str_count(test$sample3,"Normal Tissue")) #GTEx
sum(str_count(test$sample3,"Primary Tumor")) #Cancer
sum(str_count(test$sample3,"Recurrent Tumor")) #Cancer
sum(str_count(test$sample3,"Solid Tissue Normal")) #Normal

# number of different sample type in training and testing 
train %>% group_by(sample3) %>% summarize(count=n())
test %>% group_by(sample3) %>% summarize(percent_rank=n())

View(head(test))

# remove the row name of training and testing
#rownames(train)<-NULL
#rownames(test)<-NULL

x.train <- data.matrix(train[,-c(1:5)]) # Take expression values from all genes.
y.train <- data.matrix(train[,5]) # the response variable

# Save the training data
saveRDS(x.train, file="x.train.RData")
saveRDS(y.train, file="y.train.RData")


x.test <- data.matrix(test[,-c(1:5)]) # Take expression values from all genes.
y.test <- data.matrix(test[,5]) # the response variable


dim(x.train)
dim(y.train)

View(x.test)
View(y.test)
# number of Cancer Sample (1) in Training data
length(y.train[y.train == 1])

# number of Normal Sample (0) in Training data
length(y.train[y.train == 0])

# number of Cancer Sample (0) in Test data
length(y.test[y.test == 1])

# number of Normal Sample (1) in Test data
length(y.test[y.test == 0])



# ----------------Step 2: fit the model using all genes glmnet for LASSO ----------------
# fit is an object of class glmnet that contains all the relevant information of the fitted model for further use.
set.seed(199) 
fit <- glmnet(x.train, y.train, family= "binomial", alpha = 1, standardize = TRUE) #  pmax=20,
fit

# pmax= if we want to have maximum 20 genes
# alpha= “1”: for lasso regression
# alpha= “0”: for ridge regression
# alpha= a value between 0 and 1 (say 0.3) for elastic net regression.

# ----------------Step 3: visualize the coefficients of all genes by executing the plot method ----------------
par(mfrow=c(3,2)) # If you want figures in three rows, 2 column 
par(mar=c(1,1,1,1)) # If need big figure margin
#plot(fit, xvar = "lambda", label = TRUE)

pdf("LASSO_coefficient_profiles_40_gene.pdf", width = 6, height=4)
plot_glmnet(fit)
dev.off()


#--------------------Step 4(A) For computing penalized logistic regression-----------
# This results in shrinking the coefficients of the less contributing variables toward zero.
# In penalized regression, you need to specify a constant lambda to adjust the amount of the coefficient shrinkage.
# The best lambda value can be automatically get using the function cv.glmnet().

#-------------------- (Step A1): Find the best penalty lambda using cv.glmnet()  ---------------
set.seed(212) 
cv.lasso <- cv.glmnet(x.train, y.train, alpha = 1, family = "binomial", type.measure = "auc", standardize = TRUE) #type.measure = "auc"; type.measure = "deviance",
# (type.measure = "auc"), The default is type.measure="deviance",
# by default (intercept = T),
#AUC 
summary(cv.lasso$cvm)

# minimum MSE
min(cv.lasso$cvm)



# lambda.min returns the value of lambda that gives minimum mean cross-validated error. # lambda for this min MSE
cv.lasso$lambda.min

# lambda.1se is the value of lambda that gives the most regularized model such that the cross-validated error is within one standard error of the minimum.
cv.lasso$lambda.1se

# No. of coef at Min MSE
cv.lasso$nzero[cv.lasso$lambda == cv.lasso$lambda.min] 

# The best lambda value is stored inside 'cv.lasso$lambda.min'.
# identify the lambda value that produces the lowest test mean squared error (MSE).
best_lambda <- cv.lasso$lambda.min
best_lambda


#-------------------- (Step A2): Plot give Mean Square Error(MSE) across all the  λ values------------
# The numbers across the top of the plot refer to the number of features in the model.

par(mar=c(1,1,1,1))
pdf("LASSO_Plot_cross-validation_model_40_gene.pdf", width = 6, height=4)
plot(cv.lasso)
dev.off()


#-------- (Step A3A): Using lambda.min to get the regression coefficients-------
df_coef_min<- round(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min)), 5)

# See all contributing variables
coef_output_min<-as.matrix(df_coef_min[df_coef_min[, 1] != 0, ])
write.table(coef_output_min, file = "CV_Lasso_lambda-min_Model_coefficient.txt",  sep = "\t", quot=FALSE)

GeneName_min<-rownames(coef_output_min)
GeneName_min<-GeneName_min[-1]


#-------- (Step A3B): Using lambda.1se to get the regression coefficients-------
df_coef_se<- round(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.1se)), 5)

# See all contributing variables
coef_output_se<-as.matrix(df_coef_se[df_coef_se[, 1] != 0, ])
write.table(coef_output_se, file = "CV_Lasso_lambda-se_Model_coefficient.txt",  sep = "\t", quot=FALSE)

GeneName_se<-rownames(coef_output_se)
GeneName_se<-GeneName_se[-1]

#-------(Step A3C): Comparing the gene between min and 1se
#Find Elements that Exist Only in min, But Not in 1se
setdiff(GeneName_min, GeneName_se)

#----------Coefficient table and bar plot --------------------

# create a function to transform coefficient of glmnet and cvglmnet to data.frame
coeff2dt <- function(fitobject, s) {
  coeffs <- coef(fitobject, s)
  coeffs.dt <- data.frame(name = coeffs@Dimnames[[1]][coeffs@i + 1], coefficient = coeffs@x)

  # reorder the variables in term of coefficients
  return(coeffs.dt[order(coeffs.dt$coefficient, decreasing = T),])
}

coeff2dt(fitobject = cv.lasso, s = best_lambda) %>% head(20)


coeffs.table <- coeff2dt(fitobject = cv.lasso, s = best_lambda)

pdf("LASSO_coeffient_bargraph.pdf", width = 6, height=4)

ggplot(data = coeffs.table) +
  geom_col(aes(x = name, y = coefficient, fill = {coefficient > 0})) +
  xlab(label = "") +
  ggtitle(expression(paste("Lasso Coefficients with ", lambda.min, " 0.0005101641"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

dev.off()


# ------------- heatmap  on training dataset------------------
heatmapdataTrain<-(select(train, c(Class, GeneName_min))) # select the expression data according to Gene name, and Class
colnames(heatmapdataTrain)
m<-(as.matrix(heatmapdataTrain[,-1]))
rownames(m) <- heatmapdataTrain$Class


library("gplots")
# Show effect of z-score scaling within columns(single gene across different samples), blue-red color scale
# heatmap.2(m)

par(mar=c(1,2,2,2))
pdf("Heatmap_training_data_selected_17_gene.pdf", width = 12, height=6)
heatmap.2(m, col=bluered, scale="column", tracecol="#303030")
dev.off()



#--------- final model on selected 17 genes based on lambda.min--------
x2.train<-as.data.frame(x.train)
x3.train<-(select(x2.train, GeneName_min)) # select the expression data according to important Gene name
x4.train <-as.matrix(x3.train) # for cv.glment, data should be in matrix (not in data frame)
#cv.lasso_selected <- cv.glmnet(x4.train, y.train, alpha = 1, family = "binomial", type.measure = "auc", standardize = TRUE)
#plot(cv.lasso_selected)
#class(x4.train)
colnames(x3.train)
lasso.final.model.selectedgene <- glmnet(x4.train, y.train, alpha = 1, family = "binomial",
                lambda = cv.lasso$lambda.min, standardize = TRUE)

#----------- save and load the lasso model -----------------------------------
saveRDS(lasso.final.model.selectedgene, file="LASSO_model.rds")
LASSO_model<-readRDS("LASSO_model.rds")
              
#----------- predict of independt datasets ----------------------             
x2.test<-as.data.frame(x.test)
x3.test<-(select(x2.test, GeneName_min)) # select the expression data according to important Gene name
x4.test <-as.matrix(x3.test) # for cv.glment, data should be in matrix (not in data frame)                

#predict_class_selectedgene <- predict(lasso.final.model.selectedgene, type = "class", s=best_lambda, newx=x4.test)
predict_probability_selectedgene <- predict(lasso.final.model.selectedgene,  type = "response", s=best_lambda, newx=x4.test)           
predict2 <- predict(LASSO_model,  type = "response", s=best_lambda, newx=x4.test)           

###################

saveRDS(x4.train, file="x4.train.rds")
saveRDS(y.train, file="y.train.rds")
saveRDS(x4.test, file="x4.test.rds")
saveRDS(y.test, file="y.test.rds")

write.table(x4.train, file = "x_training_data.txt", sep="\t")
write.table(y.train, file = "y_training_data.txt", sep="\t")
write.table(x4.test, file = "x_test_data.txt", sep="\t")
write.table(y.test, file = "y_test_data.txt", sep="\t")

View(y.train)
#Package pROC includes function coords for calculating best threshold:
#load necessary packages
library(ggplot2)
library(pROC)

rocobj<-roc(y.test,predict_probability_selectedgene)
saveRDS(x4.test, file="x4.test.RData")
saveRDS(y.test, file="y.test.RData")

auc <- round(auc(y.test, predict_probability_selectedgene),4)

#create ROC plot for independent dataset1 (TD1)
pdf("Plot_ROC_indep_TD11.pdf", width = 6, height=6)

ggroc(rocobj, colour = 'steelblue', size = 2) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')'))
dev.off()

# creat roc plot
#ggroc(rocobj)



# spec and Sensitivity at different threshold
#coords(rocobj, seq(0, 1, 0.1))
table_TD1<-coords(rocobj, seq(0, 1, 0.1), ret=c("threshold","specificity", "sensitivity", "accuracy",
                                      "tn", "tp", "fn", "fp", "npv", "ppv", "1-specificity",
                                      "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
 
 
 write.table(table_TD1, file = "AccuracyTable_independentdataset_TD1.txt",  sep = "\t", quote=FALSE)                                     
# to get best threshold                                      
#best.coords <- coords(rocobj, "best", best.method="youden", 
#                      transpose = TRUE)
#
#coords(rocobj, "best", ret = "threshold")

#------- following is very good for AUC graphics -------
library(precrec)
precrec_obj <- evalmod(scores = predict_probability_selectedgene, labels = y.test)
par(mar=c(1,1,1,1))
pdf("ROC_curve_precrec.pdf", width = 12, height=6)
autoplot(precrec_obj)
dev.off()

# Get a data frame with AUC scores
aucs <- auc(precrec_obj)
# Use knitr::kable to display the result in a table format
knitr::kable(aucs)
write.table(aucs, file = "ROC_indepdataset_1.txt",  sep = "\t")

### Computing simple logistic regression on whole gene expression data

