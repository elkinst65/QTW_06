# MSDS7333 - Quantifying the World
# Matt Baldree
# Tom Elkins
# Ben Brock
# Austin Kelly
#
# Code from Data Science in R, Case Study 1

#################################################
# 1.2 Raw Data
#################################################
options(digits = 2)
# readlines for examination
txt = readLines("offline.final.trace.txt")
# how many lines begin with #? 5312
sum(substr(txt, 1, 1) == "#")
# how many lines in the file? 151,392
length(txt)
# split element on ';' then ='
unlist(lapply(strsplit(txt[4],";")[[1]],function(x) sapply(strsplit(x,"=")[[1]],strsplit,",")))
# prereview tokens from 4th line in file split on ';'
strsplit(txt[4], ";")[[1]]
# create tokens of line 4
tokens = strsplit(txt[4], "[;=,]")[[1]]
# view the first 10 tokens of line 4
tokens[1:10]
# grab the values
tokens[c(2, 4, 6:8, 10)]
# recorded signals within an observation
tokens[-(1:10)]
# unravel signals and turn into a 11x10 matrix
tmp = matrix(tokens[-(1:10)], ncol = 4, byrow = TRUE)
mat = cbind(matrix(tokens[c(2, 4, 6:8, 10)], nrow = nrow(tmp), ncol = 6, byrow = TRUE), tmp)
dim(mat)

# create function to process a line in the file
processLine = function(x)
{
  tokens = strsplit(x, "[;=,]")[[1]]
  tmp = matrix(tokens[-(1:10)], ncol = 4, byrow = TRUE)
  cbind(matrix(tokens[c(2, 4, 6:8, 10)], nrow = nrow(tmp), ncol = 6, byrow = TRUE), tmp)
}
# apply function to several lines
tmp = lapply(txt[4:20], processLine)
# simplify
sapply(tmp, nrow)
# stack matrices together
offline = as.data.frame(do.call("rbind", tmp))
dim(offline)

# grab all lines that do not begin with '#'
lines = txt[ substr(txt, 1, 1) != "#" ]
# process them; warning messages will appear in output that data exceed size of matrix
tmp = lapply(lines, processLine)

# modify function
processLine = function(x)
{
  tokens = strsplit(x, "[;=,]")[[1]]

  # return NULL if only 10 tokens
  if (length(tokens) == 10)
    return(NULL)

  tmp = matrix(tokens[-(1:10)], ncol = 4, byrow = TRUE)
  mat = cbind(matrix(tokens[c(2, 4, 6:8, 10)], nrow(tmp), 6, byrow = TRUE), tmp)
}

tmp = lapply(lines, processLine)
offline = as.data.frame(do.call("rbind", tmp), stringsAsFactors = FALSE)
dim(offline)

#################################################
# 1.3 Cleaning the Data and Building a Representation for Analysis
#################################################

# variable names
names(offline) = c("time", "scanMac", "posX", "posY", "posZ", "orientation", "mac",
                   "signal", "channel", "type")
# numeric variables
numVars = c("time", "posX", "posY", "posZ", "orientation", "signal")
offline[numVars] =  lapply(offline[numVars], as.numeric)

# drop all adhoc measurements and remvoe type variable
# over a 100,000 records will be removed
offline = offline[ offline$type == "3", ]
offline = offline[ , "type" != names(offline) ]
dim(offline)

# save off original time (ms)
offline$rawTime = offline$time
# convert to seconds (s)
offline$time = offline$time/1000
# calendar time
class(offline$time) = c("POSIXt", "POSIXct")
# check variable types
unlist(lapply(offline, class))
# data sanity check
summary(offline[, numVars])
# convert character variables to factors
summary(sapply(offline[, c("mac", "channel", "scanMac")],
               as.factor))
# anomalies found: only one value for scanMac, all values for posZ are 0
# discard scanMac and posZ
offline = offline[ , !(names(offline) %in% c("scanMac", "posZ"))]

###############################################
# 1.3.1 Exploring Orientation

# check number of orientations
length(unique(offline$orientation))
# look at distribution of orientation
plot(ecdf(offline$orientation))

# empirical cumulative distribution function plot
oldPar = par(mar = c(4, 4, 1, 1))
plot(ecdf(offline$orientation), pch = 19, cex = 0.3,
     xlim = c(-5, 365), axes = FALSE,
     xlab = "orientation", ylab = "Empirical CDF", main = "")
box()
axis(2)
axis(side = 1, at = seq(0, 360, by = 45))
par(oldPar)
dev.off()

# density plot
oldPar = par(mar = c(4, 4, 1, 1))
plot(density(offline$orientation, bw = 2),
     xlab = "orientation",
     main = "")
par(oldPar)
dev.off()

# round orientation to 8 orientations
roundOrientation = function(angles) {
  refs = seq(0, by = 45, length  = 9)
  q = sapply(angles, function(o) which.min(abs(o - refs)))
  c(refs[1:8], 0)[q]
}
offline$angle = roundOrientation(offline$orientation)
# check rounded values to ensure they are correct
with(offline, boxplot(orientation ~ angle, xlab = "nearest 45 degree angle", ylab = "orientation"))

###############################################
# 1.3.2 Exploring MAC Addresses

# how many unique addresses and channels?
c(length(unique(offline$mac)), length(unique(offline$channel)))
# check counts of observations for the various mac addresses
table(offline$mac)

# keep records from top seven devices
subMacs = names(sort(table(offline$mac), decreasing = TRUE))[1:7]
offline = offline[ offline$mac %in% subMacs, ]

# create a table of counts for remaining mac x channel combinations and confirm
# there is one non-zero entry in each row
macChannel = with(offline, table(mac, channel))
apply(macChannel, 1, function(x) sum(x > 0))

# since there is a one-to-one correspondence between max and channel, we can
# eliminate channel variable
offline = offline[ , "channel" != names(offline)]

###############################################
# 1.3.3 Exploring the Position of the Hand-Held Device

# how many different locations do we have data?
locDF = with(offline,
             by(offline, list(posX, posY), function(x) x))
length(locDF)

# how many locations are empty?
sum(sapply(locDF, is.null))

# drop locations that were not observed, null
locDF = locDF[ !sapply(locDF, is.null) ]
length(locDF)

# determine the number of observations recorded at each location
locCounts = sapply(locDF, nrow)

# keep position information with location
locCounts = sapply(locDF,
                   function(df)
                     c(df[1, c("posX", "posY")], count = nrow(df)))

# confirm it is a matrix
class(locCounts)
dim(locCounts)

# examine a few items: 5,500 recordings at each position
# 8 orientations x 110 replications x 7 access points = 6,160 signal strengths
locCounts[ , 1:8]

# plot total signals recorded at access point
oldPar = par(mar = c(3.1, 3.1, 1, 1))
# visualize all 166 locations
# transpose matrix
locCounts = t(locCounts)
plot(locCounts, type = "n", xlab = "", ylab = "")
text(locCounts, labels = locCounts[,3], cex = .8, srt = 45)
par(oldPar)
dev.off()

###############################################
# 1.3.4 Creating a Function to Prepare the Data

# function to process data based on previous work
readData = function(filename = 'offline.final.trace.txt',
                    subMacs = c("00:0f:a3:39:e1:c0", "00:0f:a3:39:dd:cd", "00:14:bf:b1:97:8a",
                       "00:14:bf:3b:c7:c6", "00:14:bf:b1:97:90", "00:14:bf:b1:97:8d",
                       "00:14:bf:b1:97:81")) {
  txt = readLines(filename)
  lines = txt[substr(txt, 1, 1) != "#"]
  tmp = lapply(lines, processLine)
  offline = as.data.frame(do.call("rbind", tmp),
                          stringsAsFactors = FALSE)

  names(offline) = c("time", "scanMac", "posX", "posY", "posZ", "orientation",
                     "mac", "signal", "channel", "type")

  # convert numeric values
  numVars = c("time", "posX", "posY", "orientation", "signal")
  offline[numVars] = lapply(offline[numVars], as.numeric)

  # keep only signals from access points
  offline = offline[offline$type == "3",]
  offline = offline[, "type" != names(offline)]

  # convert time to POSIX
  offline$rawTime = offline$time
  offline$time = offline$time / 1000
  class(offline$time) = c("POSIXt", "POSIXct")

  # drop scanMac, posZ, channel, and type - no info in them
  offline = offline[, !(names(offline) %in% c("scanMac", "posZ"))]

  # round orientations to nearest 45
  offline$angle = roundOrientation(offline$orientation)

  # drop more unwanted access points
  #subMacs = names(sort(table(offline$mac),decreasing=TRUE))[1:7]
  offline = offline[offline$mac %in% subMacs,]

  # drop channel
  offline = offline[, "channel" != names(offline)]

  return(offline)
}

# create data frame
offlineRedo = readData()
# sanity check
identical(offline, offlineRedo)

# what variables are global?
library(codetools)
findGlobals(readData, merge=FALSE)$variables

#################################################
# 1.4 Signal Strength Analysis
#################################################

# we have visualized and looked at stat summaries to clean and format data
# now, investigate response variable and signal strength

# compare signal strength at different orientations. higher values are stronger signals
# mac address dropped because it was identified as an extra address
oldPar = par(mar = c(3.1, 3, 1, 1))
library(lattice)
bwplot(signal ~ factor(angle) | mac, data = offline,
       subset = posX == 2 & posY == 12 & mac != "00:0f:a3:39:dd:cd",
       layout = c(2,3))

par(oldPar)
dev.off()
# summary
summary(offline$signal)

# density plot for x=24, y=4 show conditioning on angle and mac is required
oldPar = par(mar = c(3.1, 3, 1, 1))
densityplot( ~ signal | mac + factor(angle), data = offline,
             subset = posX == 24 & posY == 4 & mac != "00:0f:a3:39:dd:cd",
             bw = 0.5, plot.points = FALSE)
par(oldPar)
dev.off()

# create factor for all unique xy combinations
offline$posXY = paste(offline$posX, offline$posY, sep = "-")
# create list of data frames for each combination (xy, angle, access point)
byLocAngleAP = with(offline, by(offline, list(posXY, angle, mac), function(x) x))
# create summary stat for each data frame
signalSummary =
  lapply(byLocAngleAP,
         function(oneLoc) {
           ans = oneLoc[1, ]
           ans$medSignal = median(oneLoc$signal)
           ans$avgSignal = mean(oneLoc$signal)
           ans$num = length(oneLoc$signal)
           ans$sdSignal = sd(oneLoc$signal)
           ans$iqrSignal = IQR(oneLoc$signal)
           ans
           })
offlineSummary = do.call("rbind", signalSummary)

# box plots of stddev for signal strength
# weakest signals have smaller stddev
oldPar = par(mar = c(3.1, 3, 1, 1))
breaks = seq(-90, -30, by = 5)
bwplot(sdSignal ~ cut(avgSignal, breaks = breaks),
       data = offlineSummary,
       subset = mac != "00:0f:a3:39:dd:cd",
       xlab = "Mean Signal", ylab = "SD Signal")
par(oldPar)
dev.off()

# examine skewness of signal strength
# plot difference against number of observations
oldPar = par(mar = c(4.1, 4.1, 1, 1))
with(offlineSummary,
     smoothScatter((avgSignal - medSignal) ~ num,
                   xlab = "Number of Observations",
                   ylab = "mean - median"))
abline(h = 0, col = "#984ea3", lwd = 2)
# locally smooth differences between mean and median
lo.obj =
  with(offlineSummary,
       loess(diff ~ num,
             data = data.frame(diff = (avgSignal - medSignal),
                               num = num)))
# predict difference for each value of num and add these to the scatter plot
lo.obj.pr = predict(lo.obj, newdata = data.frame(num = (70:120)))
lines(x = 70:120, y = lo.obj.pr, col = "#4daf4a", lwd = 2)
par(oldPar)
dev.off()

###############################################
# 1.4.2 The Relationship between Signal and Distance

# one way to examine the relationship between distance and signal strength is to smooth
# the signal strength over the region where it is measured and create a contour plot.
# need to control for access point and orientation
# select one mac and one orientation to examine
oneAPAngle = subset(offlineSummary, mac == subMacs[5] & angle == 0)

# make topo map with heatmap
# fields uses splines to fit a surface to teh signal strength values at the observed locations.
# it also provides plotting routines for visualizing the surface with a heatmap.
library(fields)
# fit a smooth surface to mean signal strength
smoothSS = Tps(oneAPAngle[, c("posX","posY")], oneAPAngle$avgSignal)
# predict the value for the fitted surface at a grid of xy values
vizSmooth = predictSurface(smoothSS)
# plot predicted signal strength
plot.surface(vizSmooth, type = "C")
# add locations where the measurements were taken
points(oneAPAngle$posX, oneAPAngle$posY, pch=19, cex = 0.5)

# create contour plotting function
surfaceSS = function(data, mac, angle = 45) {
  require(fields)
  oneAPAngle = data[ data$mac == mac & data$angle == angle, ]
  smoothSS = Tps(oneAPAngle[, c("posX","posY")], oneAPAngle$avgSignal)
  vizSmooth = predictSurface(smoothSS)
  plot.surface(vizSmooth, type = "C", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  points(oneAPAngle$posX, oneAPAngle$posY, pch=19, cex = 0.5)
}
# save current settings for plotting parameters
parCur = par(mfrow = c(2,2), mar = rep(1, 4))
# make four calls to sufraceSS for macs and angle
# these plots show the following:
#   location of access point at the top of contour
#   effect of the orientation on signal strength
#   corridor effect; i.e., the signal strength is stronger relative to distance along corridors
#     where teh signals are not blocked by walls.
mapply(surfaceSS, mac = subMacs[ rep(c(5, 1), each = 2) ],
       angle = rep(c(0, 135), 2),
       data = list(data = offlineSummary))
# reset plotting parameters
par(parCur)

offlineSummary = subset(offlineSummary, mac != subMacs[2])
# create a matrix of relevant positions for the six access points on the floor plan by mac
AP = matrix( c( 7.5, 6.3, 2.5, -.8, 12.8, -2.8,
                1, 14, 33.5, 9.3,  33.5, 2.8),
            ncol = 2, byrow = TRUE,
            dimnames = list(subMacs[ -2 ], c("x", "y") ))
AP

# examine relationship between signal strength and distance from access point
# distances from locations of device emitting signal to AP receving the signal
diffs = offlineSummary[ , c("posX", "posY")] - AP[ offlineSummary$mac, ]
# find Euclidean distance between the position of the hand-held device and AP
offlineSummary$dist = sqrt(diffs[ , 1]^2 + diffs[ , 2]^2)
# scatter plots for each AP and device orientation
# curvature shows in plot; log transform might help but need to be careful with log transform of negative values
library(lattice)
oldPar = par(mar = c(3.1, 3.1, 1, 1))
xyplot(signal ~ dist | factor(mac) + factor(angle), data = offlineSummary, pch = 19, cex = 0.3,
       xlab ="distance")
par(oldPar)
dev.off()

#################################################
# 1.5 Nearest Neighbor Methods to Predict Location
#################################################

macs = unique(offlineSummary$mac)
online = readData("online.final.trace.txt", subMacs = macs)
# locations where test measurements were taken to assess accuracy
# create a unique location identifier
online$posXY = paste(online$posX, online$posY, sep = "-")
length(unique(online$posXY))
# tally the signal strength recorded at each location
tabonlineXYA = table(online$posXY, online$angle)
# output of six rows shows signal strengths were recorded at one orientation for each location
tabonlineXYA[1:6, ]

# organize data with six columns of mean signal strengths; i.e., one for each AP
keepVars = c("posXY", "posX","posY", "orientation", "angle")
byLoc = with(online,
             by(online, list(posXY),
                function(x) {
                  ans = x[1, keepVars]
                  avgSS = tapply(x$signal, x$mac, mean)
                  y = matrix(avgSS, nrow = 1, ncol = 6,
                             dimnames = list(ans$posXY, names(avgSS)))
                  cbind(ans, y)
                }))
onlineSummary = do.call("rbind", byLoc)

# data structure details
dim(onlineSummary)
names(onlineSummary)

###############################################
# 1.5.2 Choice of Orientation

# in kNN, we want to find records that have similar orientations to our new observation because orientation
# impacts strength of signal.
# m is number of neighboring same degree angles, angleNewObs is angle of new observation
m = 3; angleNewObs = 230
# reference angles
refs = seq(0, by = 45, length  = 8)
# round to one of our angles
nearestAngle = roundOrientation(angleNewObs)
# handle odd and even number of angles
if (m %% 2 == 1) {
  angles = seq(-45 * (m - 1) /2, 45 * (m - 1) /2, length = m)
} else {
  m = m + 1
  angles = seq(-45 * (m - 1) /2, 45 * (m - 1) /2, length = m)
  if (sign(angleNewObs - nearestAngle) > -1)
    angles = angles[ -1 ]
  else
    angles = angles[ -m ]
}
# map angles to values in refs
angles = angles + nearestAngle
angles[angles < 0] = angles[ angles < 0 ] + 360
angles[angles > 360] = angles[ angles > 360 ] - 360
# select observations to analyze
offlineSubset = offlineSummary[ offlineSummary$angle %in% angles, ]

# reshape signal strength help function
# aggregate signal strengths from these angles and create a data structure similar to onlineSummary
reshapeSS = function(data,
                     varSignal = "signal",
                     keepVars = c("posXY", "posX", "posY")) {
  byLocation = with(data, by(data, list(posXY),
                             function(x) {
                               ans = x[1, keepVars]
                               avgSS = tapply(x[, varSignal], x$mac, mean)
                               y = matrix(avgSS, nrow = 1, ncol = 6, dimnames = list(ans$posXY, names(avgSS)))
                               cbind(ans, y)
                             }))
  newDataSS = do.call("rbind", byLocation)
  return(newDataSS)
}
# summarize and reshape offlineSubset
trainSS = reshapeSS(offlineSubset, varSignal = "avgSignal")

# function from code above that select angles and reshape
# averages signal strengths for the different angles to produce
# one set of signal strengths for each locations in training data
selectTrain = function(angleNewObs, signals = NULL, m = 1){
  # m is the number of angles to keep between 1 and 5
  refs = seq(0, by = 45, length  = 8)
  nearestAngle = roundOrientation(angleNewObs)

  # handle odd and even m number
  if (m %% 2 == 1)
    angles = seq(-45 * (m - 1) /2, 45 * (m - 1) /2, length = m)
  else {
    m = m + 1
    angles = seq(-45 * (m - 1) /2, 45 * (m - 1) /2, length = m)
    if (sign(angleNewObs - nearestAngle) > -1)
      angles = angles[ -1 ]
    else
      angles = angles[ -m ]
  }
  # map angles to values in refs
  angles = angles + nearestAngle
  angles[angles < 0] = angles[ angles < 0 ] + 360
  angles[angles > 360] = angles[ angles > 360 ] - 360
  angles = sort(angles)
  # select observations to analyze
  offlineSubset = signals[ signals$angle %in% angles, ]
  # reshape signal strength
  reshapeSS(offlineSubset, varSignal = "avgSignal")
}

# test function for angle of 130 deg and three neighbors
train130 = selectTrain(130, offlineSummary, m = 3)

# examine output
head(train130)
# 166 locations
# what if we didn't want to collapse the signal strengths across the m angles, and
# instead return a set of mx166 signals for each access point?
length(train130[[1]])

###############################################
# 1.5.3 Finding the Nearest Neighbors

# want to look the distance in terms of signal strengths from these training data to the new
# data point. we need to calculate teh distrance from the new point to all observations in the
# training set with findNN().
# returns locations of the training observations in order of closeness to the new observation's signal strength.
findNN = function(newSignal, trainSubset) {
  diffs = apply(trainSubset[ , 4:9], 1, function(x) x - newSignal)
  dists = apply(diffs, 2, function(x) sqrt(sum(x^2)) )
  closest = order(dists)
  return(trainSubset[closest, 1:3 ])
}

# make location prediction
predXY = function(newSignals, newAngles, trainData, numAngles = 1, k = 3){
  closeXY = list(length = nrow(newSignals))

  for (i in 1:nrow(newSignals)) {
    trainSS = selectTrain(newAngles[i], trainData, m = numAngles)
    closeXY[[i]] = findNN(newSignal = as.numeric(newSignals[i, ]), trainSS)
  }
  # for some k of nearest neighbors, simply average the first k locations
  # could have used weights in the average that are inversely proportional to the distance in signal strength
  # from the test observation. will need to call findNN() to get distance to point. the weights
  # are 1/d / sum(1/d).
  # could also use a different metric besides Euclidean like Manhattan.
  # could use medians instead of averages when combining neighbors if the distribution of values are quite skewed.
  estXY = lapply(closeXY, function(x) sapply(x[ , 2:3], function(x) mean(x[1:k])))
  estXY = do.call("rbind", estXY)
  return(estXY)
}

# test prediction function with three nearest neighbors and three orientations
estXYk3 = predXY(newSignals = onlineSummary[ , 6:11],
                 newAngles = onlineSummary[ , 4],
                 offlineSummary, numAngles = 3, k = 3)
# test prediction function with one nearest neighbors and three orientations
estXYk1 = predXY(newSignals = onlineSummary[ , 6:11],
                 newAngles = onlineSummary[ , 4],
                 offlineSummary, numAngles = 3, k = 1)

# map of actual and predicted locations
# redline segments connect test locations (black dots) to their predicted locations (asterisks).
floorErrorMap = function(estXY, actualXY, trainPoints = NULL, AP = NULL){
    plot(0, 0, xlim = c(0, 35), ylim = c(-3, 15), type = "n", xlab = "", ylab = "", axes = FALSE)
    box()
    # plot AP
    if ( !is.null(AP) ) points(AP, pch = 15)
    # plot training points
    if ( !is.null(trainPoints) )
      points(trainPoints, pch = 19, col="grey", cex = 0.6)
    # plot actual point (black dot)
    points(x = actualXY[, 1], y = actualXY[, 2], pch = 19, cex = 0.8 )
    # plot predicted point (asterisk)
    points(x = estXY[, 1], y = estXY[, 2], pch = 8, cex = 0.8 )
    # redline from actual to predicted
    segments(x0 = estXY[, 1], y0 = estXY[, 2], x1 = actualXY[, 1], y1 = actualXY[ , 2], lwd = 2, col = "red")
}
# training points of average signal strengths from each of teh 166 offline locations to the six access points
trainPoints = offlineSummary[ offlineSummary$angle == 0 & offlineSummary$mac == "00:0f:a3:39:e1:c0",
                              c("posX", "posY")]

# three neighbor plot
oldPar = par(mar = c(1, 1, 1, 1))
floorErrorMap(estXYk3, onlineSummary[ , c("posX","posY")], trainPoints = trainPoints, AP = AP)
par(oldPar)
dev.off()
# one neighbor plot
oldPar = par(mar = c(1, 1, 1, 1))
floorErrorMap(estXYk1, onlineSummary[ , c("posX","posY")], trainPoints = trainPoints, AP = AP)
par(oldPar)
dev.off()

# compare fit numerically with sum of squared errors
calcError = function(estXY, actualXY) {
  sum(rowSums((estXY - actualXY) ^ 2))
}

# actual xy
actualXY = onlineSummary[ , c("posX", "posY")]
# errors: k=1 is 659 and k=3 is 307
sapply(list(estXYk1, estXYk3), calcError, actualXY)

###############################################
# 1.5.4 Cross-Validation and Choice of k

# selecting the correct number of k is a model selection problem.
# need to ensure we don't overfit the model.
# we can use v-fold cross-validation to help with the overfitting concern.
# divide data into v non-overlapping subsets of equal size.
# for each subset, build models without subset and use subset to assess model.
# repeat for all folds and aggregate errors across folds.

# we have eight orientations and six mac addresses with each location.
# we need to cross-validate on 166 locations.

# choose a v-fold of 11
v = 11
# number of locations
permuteLocs = sample(unique(offlineSummary$posXY))
# matrix locations into folds
permuteLocs = matrix(permuteLocs, ncol = v, nrow = floor(length(permuteLocs)/v))
# offline posXY in first fold
onlineFold = subset(offlineSummary, posXY %in% permuteLocs[ , 1])
# reshape data
reshapeSS = function(data, varSignal = "signal",
                     keepVars = c("posXY", "posX","posY"),
                     sampleAngle = FALSE,
                     refs = seq(0, 315, by = 45)) {
  byLocation = with(data, by(data, list(posXY), function(x) {
                    # select one angle at random for each location
                    if (sampleAngle) {
                      x = x[x$angle == sample(refs, size = 1), ]}
                    ans = x[1, keepVars]
                    avgSS = tapply(x[ , varSignal ], x$mac, mean)
                    y = matrix(avgSS, nrow = 1, ncol = 6, dimnames = list(ans$posXY, names(avgSS)))
                    cbind(ans, y)}))

  newDataSS = do.call("rbind", byLocation)
  return(newDataSS)
}

offline = offline[ offline$mac != "00:0f:a3:39:dd:cd", ]
keepVars = c("posXY", "posX","posY", "orientation", "angle")
# summarize and format offline
onlineCVSummary = reshapeSS(offline, keepVars = keepVars, sampleAngle = TRUE)
# first fold
onlineFold = subset(onlineCVSummary, posXY %in% permuteLocs[ , 1])
# training data
offlineFold = subset(offlineSummary, posXY %in% permuteLocs[ , -1])
# prediction data
estFold = predXY(newSignals = onlineFold[ , 6:11],
                 newAngles = onlineFold[ , 4],
                 offlineFold, numAngles = 3, k = 3)
# actual data
actualFold = onlineFold[ , c("posX", "posY")]
# error
calcError(estFold, actualFold)

# for each fold, we want to find k-NN estimates from 1..K and aggregate errros over the v folds.
# wrap code above in loops over folds and number of neighbors
K = 20
# errors array
err = rep(0, K)

# loop through folds
for (j in 1:v) {
  onlineFold = subset(onlineCVSummary, posXY %in% permuteLocs[ , j])
  offlineFold = subset(offlineSummary, posXY %in% permuteLocs[ , -j])
  actualFold = onlineFold[ , c("posX", "posY")]

  # loop through neighbors
  for (k in 1:K) {
    estFold = predXY(newSignals = onlineFold[ , 6:11],
                     newAngles = onlineFold[ , 4],
                     offlineFold, numAngles = 3, k = k)
    err[k] = err[k] + calcError(estFold, actualFold)
  }
}

# plot results as sum of squared errors as function of k
oldPar = par(mar = c(4, 3, 1, 1))
plot(y = err, x = (1:K),  type = "l", lwd= 2,
     ylim = c(1200, 2100),
     xlab = "Number of Neighbors",
     ylab = "Sum of Square Errors")

rmseMin = min(err)
kMin = which(err == rmseMin)[1]
segments(x0 = 0, x1 = kMin, y0 = rmseMin, col = gray(0.4), lty = 2, lwd = 2)
segments(x0 = kMin, x1 = kMin, y0 = 1100,  y1 = rmseMin, col = grey(0.4), lty = 2, lwd = 2)
text(x = kMin - 2, y = rmseMin + 40, label = as.character(round(rmseMin)), col = grey(0.4))
par(oldPar)
dev.off()

# use k=5 and apply to original training and test data
estXYk5 = predXY(newSignals = onlineSummary[ , 6:11],
                 newAngles = onlineSummary[ , 4],
                 offlineSummary, numAngles = 3, k = 5)
calcError(estXYk5, actualXY)

# revision to prediction function that might speed it up
# return all k estimates and use sumsum() to provide k means; i.e., cumsum(x[1-:K])/(1:K)
predXY = function(newSignals, newAngles, trainData, numAngles = 1, k = 3){
  closeXY = list(length = nrow(newSignals))

  for (i in 1:nrow(newSignals)) {
    trainSS = selectTrain(newAngles[i], trainData, m = numAngles)
    closeXY[[i]] = findNN(newSignal = as.numeric(newSignals[i, ]), trainSS)
  }

  estXY = lapply(closeXY, function(x){ sapply(x[ , 2:3], function(x) mean(x[1:k]))})
  estXY = do.call("rbind", estXY)
  return(estXY)
}

# did not use cross-validation to determine optimal angles to use.
