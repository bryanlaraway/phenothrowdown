install.packages("ROCR")
install.packages("rjson")
install.packages("ggplot2")
library(ROCR)
library(rjson)
library(ggplot2)





mydata_mix <- read.csv("/Volumes/Time Machine/PycharmProjects/phenothrowdown/inter/omim/morbid_disease_predictions_true_false.csv")  # read csv file 
mydata_mix

mydata_ranks <- read.csv("/Volumes/Time Machine/PycharmProjects/phenothrowdown/inter/omim/morbid_disease_predictions_with_rankings.csv")  # read csv file 
mydata_ranks

newdata <- mydata_ranks[order(mydata_ranks$top_phenolog_max_rank),]
newdata$top_phenolog_max_rank_percentage <- mydata_ranks$top_phenolog_max_rank.rank/nrow(newdata$top_phenolog_max_rank)

mydata <- read.csv("/Volumes/Time Machine/PycharmProjects/phenothrowdown/inter/omim/morbid_disease_predictions.csv")  # read csv file 
mydata
pred <- prediction(mydata$top_owlsim_iccs_score, mydata$category)
perf1 <- performance(pred, "sens", "spec")
plot(perf1)



"mydata$category"
# mydata[,category] <- [, mydata$category] != 0
pred <- prediction(mydata_mix$top_owlsim_iccs_score, mydata_mix$category)
perf1 <- performance(pred, "sens", "spec")
plot(perf1)
perf <- performance(pred,"tpr","fpr")
plot(perf)
newdata <- mydata_mix[order(mydata_mix$category),]
newdata


json_file <- "/Volumes/Time Machine/PycharmProjects/phenothrowdown/out/roc/top_phenolog_max_score_list.json"
predictions <- fromJSON(file=json_file)
json_file2 <- "/Volumes/Time Machine/PycharmProjects/phenothrowdown/out/roc/list_of_ones.json"
labels <- fromJSON(file=json_file2)
labels <- labels != 1
pred <- prediction(predictions, labels)
perf1 <- performance(data2, "prec", "rec")
plot(perf1)
print(predictions)
print(labels)
print(ROCR.simple)
data(ROCR.simple)
pred <- prediction( ROCR.simple$predictions, ROCR.simple$labels)
perf1 <- performance(pred, "sens", "spec")
plot(perf1)
perf <- performance(pred,"tpr","fpr")
plot(perf)

library(ggplot2)

diamonds$is_expensive <- diamonds$price > 2400
is_test <- runif(nrow(diamonds)) > 0.75
train <- diamonds[is_test==FALSE,]
test <- diamonds[is_test==TRUE,]

summary(fit <- glm(is_expensive ~ carat + cut + clarity, data=train))
prob <- predict(fit, newdata=test, type="response")
pred <- prediction(prob, test$is_expensive)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
# I know, the following code is bizarre. Just go with it.
auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]

roc.data <- data.frame(fpr=unlist(perf@x.values),
                       tpr=unlist(perf@y.values),
                       model="GLM")
ggplot(roc.data, aes(x=fpr, ymin=0, ymax=tpr)) +
  geom_ribbon(alpha=0.2) +
  geom_line(aes(y=tpr)) +
  ggtitle(paste0("ROC Curve w/ AUC=", auc))
