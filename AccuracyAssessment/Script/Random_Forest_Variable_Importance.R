library(randomForest)
library(caret)
library(ggplot2)
setwd('/home/ate/rrd/')

# ModelDataset = read.delim('Training_Imperv-Crop-Other_SR_20170713.csv',header = TRUE,sep = ',',na.strings=c(""))
# ModelDataset = read.delim('Training_Imperv-Crop-Other_TOA_20170714.csv',header = TRUE,sep = ',')
# ModelDataset = read.delim('Training_Phenology_TOA_20170714.csv',header = TRUE,sep = ',')
# ModelDataset = read.delim('Training_Grass-Shrub-Tree_TOA_20170714.csv',header = TRUE,sep = ',')
# ModelDataset = read.delim('Training_Leaf_Type_TOA_20170718.csv',header = TRUE,sep = ',')
ModelDataset = read.delim('sampledata4.csv',header = TRUE,sep = ',')



ModelDataset$rand <- sample(100, size = nrow(ModelDataset), replace = TRUE)


# ModelDataset = read.delim('Training_Autocorr.csv',header = TRUE,sep = ',')
# ModelDataset = read.delim('Training_agg.csv',header = TRUE,sep = ',')
#drops <- c("ndwi","sarVH","sarVV","blue","evi","ndvi","nir","swir","nd2")
#ModelDataset = ModelDataset[ , !(names(ModelDataset) %in% drops)]
ModelDataset = ModelDataset[duplicated(ModelDataset), ]
ModelDataset = ModelDataset[complete.cases(ModelDataset), ]
ModelDataset = ModelDataset[ModelDataset$rand < 40,]

#for (i in seq(1,length(ModelDataset$land_class),1)){
#  if (ModelDataset$land_class[i] < 0){
#    ModelDataset$land_class[i] = 0
#  }
#  if (ModelDataset$land_class[i] > 120){
#    ModelDataset$land_class[i] = 120
#  }
#  }

#ModelDataset = ModelDataSet[ModelDataSet$land_class > 0,]
#ModelDataset = ModelDataSet[ModelDataSet$land_class <150,]
ModelDataset$land_class = round(ModelDataset$land_class / 10)
# ModelDataset <- ModelDataset[ModelDataset$land_class != 2,] #phenology
# ModelDataset <- ModelDataset[ModelDataset$land_class != 3,] #leaf type
ModelDataset$land_class <- factor(ModelDataset$land_class)


# correlationMatrix <- cor(ModelDataset[,c(c(1:103),c(105:160))]) # others
#correlationMatrix <- cor(ModelDataset[,c(2:12)]) # mangrove, autocorr
# correlationMatrix <- cor(ModelDataset[,c(c(1:106),c(108:163))]) # agg
#highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.9, names = TRUE, exact = TRUE)
#print(highlyCorrelated)
#ModelDataset <- ModelDataset[ , -which(names(ModelDataset) %in% highlyCorrelated)]


RF_ModelTOA = randomForest(land_class~., data=ModelDataset,importance=TRUE,ntree=100,na.action=na.omit)
varImpPlot(RF_ModelTOA)
print(RF_ModelTOA)


RF_import <- importance(RF_ModelTOA)
RF_import <- as.data.frame(RF_import)
RF_import <- RF_import[order(-RF_import$MeanDecreaseAccuracy),]
RF_import <- cbind(rownames(RF_import), RF_import)
rownames(RF_import) <- NULL
colnames(RF_import)[1] <- 'covariates'

ggplot(data = RF_import, aes(x = reorder(covariates,-MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) + 
  xlab('covariates') + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data = RF_import, aes(x = reorder(covariates,-MeanDecreaseGini), y = MeanDecreaseGini)) + 
  xlab('covariates') + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data = RF_import, aes(x = MeanDecreaseAccuracy, y = MeanDecreaseGini)) +
  geom_point()

ggplot(data = RF_import, aes(x = MeanDecreaseAccuracy, y = MeanDecreaseGini, label = covariates)) +
  geom_point() + geom_text(aes(label = covariates),hjust=0, vjust=0)

pc_results <- princomp(RF_import[,c('MeanDecreaseAccuracy','MeanDecreaseGini')])
pc_scores <- as.data.frame(pc_results$scores)
pc_scores$covariates <- RF_import$covariates
# pc_scores$Comp.1 <- pc_scores$Comp.1*(-1) # autocorr, agg
pc_scores <- pc_scores[order(-pc_scores$Comp.1),]

ggplot(data = pc_scores, aes(x = reorder(covariates,-Comp.1), y = Comp.1)) + 
  xlab('covariates') + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data = pc_scores, aes(x = Comp.1, y = Comp.2, label = covariates)) +
  geom_point() + geom_text(aes(label = covariates),hjust=0, vjust=0)

out <- as.character(pc_scores$covariates[pc_scores$Comp.1 > 0])

cat(paste(shQuote(out, type="cmd"), collapse=", "))

ModelFinal <- ModelDataset[ , which(names(ModelDataset) %in% c(out,'land_class'))]

RF_ModelFinal = randomForest(land_class~., data=ModelFinal,importance=TRUE,ntree=500,na.action=na.omit)
varImpPlot(RF_ModelFinal)
print(RF_ModelFinal)

write.csv(ModelDataset,"/home/ate/rrd/riceData.csv")
