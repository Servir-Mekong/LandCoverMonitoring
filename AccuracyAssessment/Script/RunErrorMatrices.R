#################################################
## Accuracy Assessment
#################################################
CEO <- read.csv('.\\AccuracyAssessment\\Data\\Validation_20170729.csv')
colnames(CEO)[49]<-'ephemeral'

colnames(CEO)

## Cropland, Built, Rice, Other
hist(CEO$cropland)
hist(CEO$CROP)

plot(CROP~cropland, data = CEO, 
     ylab = 'Fractional Crop Cover', xlab = 'Prob of Crop Cover')
lw1 <- loess(CROP ~ cropland, data = CEO, span=0.2)

my.count <- seq(from=1, to=101, by=1)
predd <- predict(lw1, my.count, se=TRUE) 
lines(predd$fit, lty = 'solid', col = 'darkred', lwd=3)

summary(lm(CROP ~ cropland, data = CEO))
