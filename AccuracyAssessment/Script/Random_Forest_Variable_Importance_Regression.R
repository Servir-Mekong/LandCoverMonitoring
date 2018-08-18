library(randomForest)
library(caret)
library(ggplot2)
setwd('C:\\Users\\joshuagoldstein\\Downloads')

ModelDataset = read.delim('Training_PercentImpervious_TOA_20170714.csv',header = TRUE,sep = ',')

ModelDataset$drycool_count <- NULL
ModelDataset$dryhot_count <- NULL
ModelDataset$rainy_count <- NULL

correlationMatrix <- cor(ModelDataset[,c(c(1:103),c(105:160))]) # others
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.9, names = TRUE, exact = TRUE)
print(highlyCorrelated)
ModelDataset <- ModelDataset[ , -which(names(ModelDataset) %in% highlyCorrelated)]

ModelDataset$percent <- NULL
RF_ModelTOA = randomForest(logratio~., data=ModelDataset,importance=TRUE,ntree=500,na.action=na.omit)
varImpPlot(RF_ModelTOA)
print(RF_ModelTOA)

RF_import <- importance(RF_ModelTOA)
RF_import <- as.data.frame(RF_import)
RF_import <- RF_import[order(-RF_import[["%IncMSE"]]),]
RF_import <- cbind(rownames(RF_import), RF_import)
rownames(RF_import) <- NULL
colnames(RF_import)[1] <- 'covariates'


ggplot(data = RF_import, aes(x = reorder(covariates,-RF_import[["%IncMSE"]]), y = RF_import[["%IncMSE"]])) + 
  xlab('covariates') + ylab('%IncMSE') + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data = RF_import, aes(x = reorder(covariates,-IncNodePurity), y = IncNodePurity)) + 
  xlab('covariates') + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data = RF_import, aes(x = RF_import[["%IncMSE"]], y = IncNodePurity)) +
  geom_point()

ggplot(data = RF_import, aes(x = RF_import[["%IncMSE"]], y = IncNodePurity, label = covariates)) +
  geom_point() + geom_text(aes(label = covariates),hjust=0, vjust=0)

pc_results <- princomp(RF_import[,c('%IncMSE','IncNodePurity')], cor=TRUE)
pc_scores <- as.data.frame(pc_results$scores)
pc_scores$covariates <- RF_import$covariates
pc_scores <- pc_scores[order(-pc_scores$Comp.1),]


ggplot(data = pc_scores, aes(x = reorder(covariates,-Comp.1), y = Comp.1)) + 
  xlab('covariates') + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data = pc_scores, aes(x = Comp.1, y = Comp.2, label = covariates)) +
  geom_point() + geom_text(aes(label = covariates),hjust=0, vjust=0)

out <- as.character(pc_scores$covariates[pc_scores$Comp.1 > 0])

cat(paste(shQuote(out, type="cmd"), collapse=", "))

ModelFinal <- ModelDataset[ , which(names(ModelDataset) %in% c(out,'logratio'))]

RF_ModelFinal = randomForest(logratio~., data=ModelFinal,importance=TRUE,ntree=500,na.action=na.omit)
varImpPlot(RF_ModelFinal)
print(RF_ModelFinal)
