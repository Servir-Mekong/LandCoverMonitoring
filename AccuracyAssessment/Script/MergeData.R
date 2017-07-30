## Pull in the CEO mapathon files.
files = list.files(path = '.\\AccuracyAssessment\\Data\\CountryFiles', pattern="*.csv")

## Merge files into one csv.
myfiles = do.call(rbind, lapply(files, function(x) read.csv(x, stringsAsFactors = FALSE)))

## Change column names before importing into EE.
colnames(myfiles)<-c("PLOTID", "longitude", "latitude", "SIZEM", "SHAPE", 
                     "FLAGGED", "ANALYSES", "SAMPLEPOINTS", "USERID", 
                     "BUILTSURFACE", "BUILTVEGNONTREE", "BUILTTREE", "MINING", 
                     "MUDFLATBEACH", "BARRENOTHER", "TREEPLANTATIONORCHARD", 
                     "TREEMANGROVE", "TREEOTHER", "SHRUB", "GRASS", "CROP", 
                     "AQUACULTUREPOND", "AQUATICVEGOTHER", "WATER", 
                     "SNOWICE", "UNKNOWN", "OTHER")

## Export one csv.
write.csv(myfiles, 'allPoints.csv', row.names = F)