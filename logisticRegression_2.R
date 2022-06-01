# This program is to develop logistic regression model on training data
# and check the performance of the model on test data
# Training and Test data contains the class as first column and 
# gene name as rest of the columns. Each row contains gene expression value from patients/normal sample.

#------[1] data using 17 gene expression ----------------
# upload the data
train_17<-readRDS("train_17.rds")
test_17<-readRDS("test_17.rds")

str(train_17)
str(test_17)

#------[2] data using 40 gene expression -----------------
# upload the data
train_40<-readRDS("train_40.rds")
test_40<-readRDS("test_40.rds")

str(train_40)
str(test_40)

#------[3] logistic regression model on 40 genes data------

# (a) Model construction using train data   
model_glm_40 <- glm( class ~ ., data = train_40, family = binomial)

# (b) Summarize the model
summary(model_glm_40)
summary(model_glm_40)$coef

# print the summary
sink("model_glm_40_summary.txt")
print(summary(model_glm_40))
sink()  # returns output to the console


# (c) Make predictions using test data
probabilities_glm_40 <- model_glm_40 %>% predict(test_40, type = "response")
#predicted.classes <- ifelse(probabilities > 0.5, "pos", "neg")

library(ggplot2)
library(pROC)
library(precrec)

precrec_obj_glm_40 <- evalmod(scores = probabilities_glm_40, labels = test_40$class)
par(mar=c(1,1,1,1))
pdf("ROC_test_glm_40.pdf", width = 12, height=6)
autoplot(precrec_obj_glm_40)
dev.off()

# Get a data frame with AUC scores
aucs_glm_40 <- auc(precrec_obj_glm_40)
# Use knitr::kable to display the result in a table format
knitr::kable(aucs_glm_40)
write.table(aucs_glm_40, file = "ROC_glm_test_40.txt",  sep = "\t")

#------END of [3] logistic regression model on 40 genes data------


#------[4] logistic regression model on 17 genes data------

# (a) Model construction using train data   
model_glm_17 <- glm( class ~ ., data = train_17, family = binomial)

# (b) Summarize the model
summary(model_glm_17)
summary(model_glm_17)$coef

# print the summary
sink("model_glm_17_summary.txt")
print(summary(model_glm_17))
sink()  # returns output to the console

# (c) Make predictions using test data
probabilities_glm_17 <- model_glm_17 %>% predict(test_17, type = "response")
#predicted.classes <- ifelse(probabilities > 0.5, "pos", "neg")

library(ggplot2)
library(pROC)
library(precrec)

precrec_obj_glm_17 <- evalmod(scores = probabilities_glm_17, labels = test_17$class)
par(mar=c(1,1,1,1))
pdf("ROC_test_glm_17.pdf", width = 12, height=6)
autoplot(precrec_obj_glm_17)
dev.off()

# Get a data frame with AUC scores
aucs_glm_17 <- auc(precrec_obj_glm_17)
# Use knitr::kable to display the result in a table format
knitr::kable(aucs_glm_17)
write.table(aucs_glm_17, file = "ROC_glm_test_17.txt",  sep = "\t")

#------END of [4] logistic regression model on 17 genes data------
