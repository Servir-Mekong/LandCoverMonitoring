#################################################
## Accuracy Assessment
#################################################
library(MASS)
library(RColorBrewer)
library(ggplot2)
library(segmented)

## Load in data.
CEOfull <- read.csv('.\\AccuracyAssessment\\Data\\Validation_20170729.csv')
colnamesfull(CEO)[49]<-'ephemeral'

colnames(CEOfull)
CEOfull$ID <- seq(1:length(CEOfull$ANALYSES))
  
xdata <- CEO$tcc
ydata <- CEO$TREEPLANTATIONORCHARD + CEO$TREEOTHER + CEO$TREEMANGROVE + CEO$BUILTTREE
xlabel <- "Predicted: Fractional Tree Cover"
ylabel <- "Validation Data: Fractional Tree Cover"

xdata <- CEO$barren
ydata <- CEO$BARRENOTHER+CEO$MINING+CEO$MUDFLATBEACH
ID <- CEOfull$ID

xlabel <- "Predicted: Prob of Barren (including mudflat, beaches, mining)"
ylabel <- "Validation Data: Fractional Barren Cover"

#############################################
## sample the data by strata
CEO_sampleLower <- cbind(ydata[xdata<=50], 
                         xdata[xdata<=50],
                         ID[xdata<=50])
colnames(CEO_sampleLower) <- c('ydata','xdata','ID')

CEO_sampleUpper <- cbind(ydata[xdata>50], 
                         xdata[xdata>50],
                         ID[xdata>50])

colnames(CEO_sampleUpper) <- c('ydata','xdata','ID')

sampleSize <- min(length(CEO_sampleLower[ , 1]),
                  length(CEO_sampleUpper[ , 1]))

LowerID <- sample(CEO_sampleLower[ ,'ID'], sampleSize)
UpperID <- sample(CEO_sampleUpper[ ,'ID'], sampleSize)
sampleData<-c(LowerID, UpperID)
length(sampleData)
sampleData<-cbind(sampleData,sampleData)
dim(sampleData)
colnames(sampleData)<-c('sampleData','sampleData2')
head(sampleData)

dataUpper <-  merge(CEO_sampleUpper, sampleData, 
                    by.x = 'ID', by.y = 'sampleData')
dataLower <-  merge(CEO_sampleLower, sampleData, 
                    by.x = 'ID', by.y = 'sampleData')

CEO_sample <- rbind(dataUpper[,-1], dataLower[,-1])
head(CEO_sample)

colnames(CEO_sample)[3]<-'ID'

#
##############################
# Segmented regression
# y ~ x
lin.mod <- lm(ydata ~ xdata, data = CEO_sample)

segmented.mod <- segmented(lin.mod, seg.Z = ~xdata, psi=20)
lines(segmented.mod, col = "blue")
segmented.mod

segmented.mod <- segmented(lin.mod, seg.Z = ~xdata, psi=c(20,40))
lines(segmented.mod, col = 'red')
segmented.mod

segmented.mod$psi[ , 2]

#fit models to each segment.
seg1Subset <- CEO_sample[CEO_sample$xdata <= segmented.mod$psi[1,2], 1:2]
head(seg1Subset)
lmMod1 <- lm(ydata ~ xdata, data = seg1Subset)
segment1_y1 <- lmMod1$coefficients[1] + lmMod1$coefficients[2]*1
segment1_y2 <- lmMod1$coefficients[1] + lmMod1$coefficients[2]*segmented.mod$psi[1, 2]

if(length(segmented.mod$psi[ , 2]) == 1){
  seg2Subset <- CEO_sample[CEO_sample$xdata > segmented.mod$psi[1, 2], 
                           c('ydata' , 'xdata')]
  head(seg2Subset)
  
  lmMod2 <- lm(ydata ~ xdata, data = seg2Subset)
  segment2_y1 <- lmMod2$coefficients[1] + lmMod2$coefficients[2]*segmented.mod$psi[1, 2]
  segment2_y2 <- lmMod2$coefficients[1] + lmMod2$coefficients[2]*100
  
  } else {
    seg2Subset <- CEO_sample[CEO_sample$xdata > segmented.mod$psi[1, 2]
                             & CEO_sample$xdata <= segmented.mod$psi[2, 2], 
                             c('ydata' , 'xdata')]

    lmMod2 <- lm(ydata ~ xdata, data = seg2Subset)
    segment2_y1 <- lmMod2$coefficients[1] + lmMod2$coefficients[2]*segmented.mod$psi[1, 2]
    segment2_y2 <- lmMod2$coefficients[1] + lmMod2$coefficients[2]*segmented.mod$psi[2, 2]
  
    seg3Subset <- CEO_sample[CEO_sample$xdata > segmented.mod$psi[2, 2], 
                             c('ydata' , 'xdata')]
    
    lmMod3 <- lm(ydata ~ xdata, data = seg3Subset)
    segment3_y1 <- lmMod3$coefficients[1] + lmMod3$coefficients[2]*segmented.mod$psi[2, 2]
    segment3_y2 <- lmMod3$coefficients[1] + lmMod3$coefficients[2] * 100
  }

###############################
# plot data
#commonTheme = list(labs(color="Density",fill="Density", x=xlabel, y=ylabel),
#                   theme_bw(), theme(legend.position=c(0,1),
#                         legend.justification=c(0,1)))

commonTheme = list(labs(color="Density",fill="Density",
                        x=xlabel,
                        y=ylabel),
                   theme_bw(),
                   theme(legend.position='none'))

plot(CEO_sample$ydata~ CEO_sample$xdata)
ggplot(data = CEO_sample ,aes(CEO_sample$xdata, CEO_sample$ydata)) + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  #stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE) + 
#  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="lightblue",high="darkred") +
  guides(alpha="none") +
  geom_vline(xintercept = segmented.mod$psi[1, 2], colour = 'black', lwd = 1.5, lty = 2) +
  #geom_vline(xintercept = segmented.mod$psi[2, 2], colour = 'black', lwd = 1.5, lty = 2) +
  commonTheme +coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) 
#+
#  geom_segment(aes(x = 1, y = segment1_y1, 
#                   xend = segmented.mod$psi[1, 2], yend = segment1_y2)) +
#  geom_segment(aes(x = segmented.mod$psi[1, 2], y = segment2_y1, 
#                   xend = segmented.mod$psi[2, 2], yend = segment2_y2))+
#  geom_segment(aes(x = segmented.mod$psi[2, 2], y = segment3_y1, 
#                 xend = 100, yend = segment3_y2))+
#  stat_smooth(method=lm) + geom_point()

summary(lm(ydata ~ xdata, data = CEO_sample))
summary(lmMod1)
summary(lmMod2)
summary(lmMod3)












########################################################
########################################################
CEO_sample$partition<- 1
CEO_sample$partition[CEO_sample$xdata > segmented.mod$psi[1, 2] & 
                       CEO_sample$xdata <= segmented.mod$psi[2, 2]]<- 2
CEO_sample$partition[CEO_sample$xdata > segmented.mod$psi[2, 2]]<- 3

CEO_sample <- merge(CEO_sample, 
                    CEOfull[,c('ID',"WATER", "AQUACULTUREPOND", 
                               "AQUATICVEGOTHER", "BARRENOTHER", 
                               "MINING", "MUDFLATBEACH", 
                              "BUILTSURFACE","BUILTTREE", "BUILTVEGNONTREE", 
                              "CROP", "GRASS", "SHRUB", "OTHER", 
                              "SNOWICE", 
                              "TREEOTHER", "TREEPLANTATIONORCHARD", "TREEMANGROVE", 
                              "UNKNOWN")], by.x = 'ID', by.y = 'ID')
CEOfull$tree<-CEOfull$TREEOTHER + CEOfull$TREEPLANTATIONORCHARD + CEOfull$TREEMANGROVE
CEOfull$wet <-  CEOfull$WATER + CEOfull$AQUACULTUREPOND
CEOfull$built<-CEOfull$BUILTSURFACE + CEOfull$BUILTTREE + CEOfull$BUILTVEGNONTREE 

CEOfull$barren<- CEOfull$BARRENOTHER + CEOfull$MINING + CEOfull$MUDFLATBEACH
CEOfull$grassShrub <- CEOfull$GRASS + CEOfull$SHRUB

# start at 5
i = 5
colnames(CEO_sample)[i]
CEO_sample[CEO_sample$partition == 1 & CEO_sample$CROP>50, i]
i = i + 1

# start at 5
i = 5
colnames(CEO_sample)[i]
CEO_sample[CEO_sample$partition == 3 & CEO_sample$CROP<50, i]
i = i + 1

dim(CEO_sample)
"SNOWICE", 
"AQUATICVEGOTHER", 

y2label <- 'Validation Data: Fractional Cover of All Built, including gardens'
commonTheme = list(labs(color="Density",fill="Density",
                        x=xlabel, y=y2label), theme_bw(),
                   theme(legend.position='none'))

ggplot(data = CEOfull[CEOfull$built>0,],
       aes(CEOfull$builtup[CEOfull$built>0], CEOfull$built[CEOfull$built>0])) + 
  #stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE) + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  scale_fill_continuous(low="red",high="black") +
  guides(alpha="none") +
  geom_vline(xintercept = segmented.mod$psi[1, 2], colour = 'black', lwd = 1.5, lty = 2) +
  geom_vline(xintercept = segmented.mod$psi[2, 2], colour = 'black', lwd = 1.5, lty = 2) +
  commonTheme +coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) 


y2label <- 'Validation Data: Fractional Cover of Cropland Surfaces'
commonTheme = list(labs(color="Density",fill="Density",
                        x=xlabel, y=y2label), theme_bw(),
                   theme(legend.position='none'))

ggplot(data = CEOfull[CEOfull$CROP>0,],
       aes(CEOfull$builtup[CEOfull$CROP>0], 
           CEOfull$CROP[CEOfull$CROP>0])) + 
  #stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE) + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  scale_fill_continuous(low="yellow",high="darkred") +
  guides(alpha="none") +
  geom_vline(xintercept = segmented.mod$psi[1, 2], colour = 'black', lwd = 1.5, lty = 2) +
  geom_vline(xintercept = segmented.mod$psi[2, 2], colour = 'black', lwd = 1.5, lty = 2) +
  commonTheme +coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) 
