WorkingDirectory <-'\\\\166.2.126.227\\International_Programs\\Mekong 2016\\githubAccount\\LandCoverMonitoring\\AccuracyAssessment'

setwd(WorkingDirectory)

AllCountryfiles <- list.files(path = '.\\Data\\CountryFiles', pattern = ".csv", full.names = FALSE)
AllData <- NULL

for (i in 1:length(AllCountryfiles)){
  dataTemp <- read.csv(paste0('.\\Data\\CountryFiles\\',AllCountryfiles[i]))
  dataTemp$Country<-substr(AllCountryfiles[i], 11, nchar(AllCountryfiles[i])-15)
  
  if (i == 1){
    AllData<- dataTemp
  }
  
  AllData <-rbind(AllData, dataTemp)
}

head(AllData)
unique(AllData$Country)

## Change column names before importing into EE.
colnames(myfiles)<-c("PLOTID", "longitude", "latitude", "SIZEM", "SHAPE", 
                     "FLAGGED", "ANALYSES", "SAMPLEPOINTS", "USERID", 
                     "BUILTSURFACE", "BUILTVEGNONTREE", "BUILTTREE", "MINING", 
                     "MUDFLATBEACH", "BARRENOTHER", "TREEPLANTATIONORCHARD", 
                     "TREEMANGROVE", "TREEOTHER", "SHRUB", "GRASS", "CROP", 
                     "AQUACULTUREPOND", "AQUATICVEGOTHER", "WATER", 
                     "SNOWICE", "UNKNOWN", "OTHER", "Country")


## Export as csv.
write.csv(AllData, '.\\Data\\ceo-rlcms-Mekong.csv', row.names = F)