#################################################
## Accuracy Assessment
#################################################
library(MASS)
library(RColorBrewer)
library(ggplot2)
library(segmented)

## Load in data.
CEO <- read.csv('.\\AccuracyAssessment\\Data\\Validation_20170729.csv')
colnames(CEO)[49]<-'ephemeral'

colnames(CEO)

## Cropland, Built, Rice, Other
hist(CEO$surface_water)
hist(CEO$WATER)

xdata <- CEO$cropland
#ydata <- CEO$TREEPLANTATIONORCHARD + CEO$TREEOTHER + CEO$TREEMANGROVE + CEO$BUILTTREE
ydata <- CEO$CROP
xlabel <- "predicted tree canopy cover"
ylabel <- "reference data Tree Canopy Cover"

plot(ydata ~ xdata, data = CEO, 
     ylab = ylabel, xlab = xlabel)
lw1 <- loess(ydata ~ xdata)

my.count <- seq(from=1, to=101, by=1)
predd <- predict(lw1, my.count, se=TRUE) 
lines(predd$fit, lty = 'solid', col = 'darkred', lwd=3)

summary(lm(ydata ~ xdata, data = CEO))

###############################
# Segmented regression
# y ~ x
lin.mod <- lm(ydata ~ xdata, data = CEO)
segmented.mod <- segmented(lin.mod, seg.Z = ~xdata, psi=c(20,40))
lines(segmented.mod, col = 'red')
segmented.mod
segmented.mod <- segmented(lin.mod, seg.Z = ~xdata, psi=20)
lines(segmented.mod, col = "blue")
segmented.mod
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

ggplot(data = CEO ,aes(xdata, ydata)) + 
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="green",high="red") +
  guides(alpha="none") +
  geom_point() + commonTheme #+ theme(legend.position="none")


