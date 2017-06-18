# MSDS7333 - Quantifying the World
# Matt Baldree
# Tom Elkins
# Ben Brock
# Austin Kelly
#
# Assignment:
# 1. Conduct a more thorough data analysis into these two MAC addresses including determining locations
# by using data corresponding to both MAC addresses.
#   Which of these two MAC addresses should be used and which should not be used for RTLS?
#   Which MAC address yields the best prediction of location?
#   Does using data for both MAC addresses simultaneously yield more, or less, accurate prediction of location?
#
# 2. Implement alternative k-nearest prediction method using weights on received signal strength.
#   For what range of values of weights are you able to obtain better prediction values than for the unweighted approach?
#     Use calcError() to compare approaches.
#
# Write up assignment in an Ipython Notebook with an introduction, explanation of approaches and outputs.
#

# two digit precision for printing and formatting
options(digits = 2)
options(error = recover, warn = 1)
# include necessary functions
source("MSDS7333-baldree-case6-fx.r", print.eval=TRUE)

#### Part 1: Access point MAC address analysis of 00:0f:a3:39:e1:c0 and 00:0f:a3:39:dd:cd ####

# read offline data
offline = readData()

# create factor for all unique xy combinations
offline$posXY = paste(offline$posX, offline$posY, sep = "-")

# access points (AP) addresses in question
# note that both devices are Alpha Networks devices because they begin with '00:0f'
chosenAP = "00:0f:a3:39:e1:c0"
rejectedAP = "00:0f:a3:39:dd:cd"
remainingAP = unique(offline[!offline$mac %in% c(chosenAP, rejectedAP),]$mac)

#### Signal Collection Map ####

# plot the locations and count of signals recorded for both access points
# we would plot the APs on the map but we do not know their X,Y locations
# result: map looks equivalent
plotSignalMaps(c(chosenAP, rejectedAP))

#### Signal Strength Distribution ####

# plot signal strength distribution for each angle per AP
# we are using the same stationary point of 2, 12 as the book
# result: rejected AP has a lower signal strength for the X,Y than chosen AP.
df = subset(offline, posX == 2 & posY == 12 & !(mac %in% remainingAP), c(signal, angle, mac))
plotBoxplotSignalStrength(df)

# examine AP strength from a fixed location at opposite end of building say 33, 3
# result: rejected AP has a lower signal strength for X,Y than chosen AP.
df = subset(offline, posX == 33 & posY == 3 & !(mac %in% remainingAP), c(signal, angle, mac))
plotBoxplotSignalStrength(df)

#### Signal Density ####

# plot signal strength density for each angle per AP
# we are using the same stationary point of 24, 4 as the book
# i'm not sure why they are using this point.
# result: rejected AP has more density at lower strength.
#         chosen AP has more density coverage at higher strength.
#         many distributions look normal. several are skewed, secondary modes, and departure from norm.
#         median of distribution vary by angle. angle makes a difference in strength.
#         it may be best to keep use both points to have coverage for all range of strengths.
df = subset(offline, posX == 24 & posY == 4 & !(mac %in% remainingAP), c(signal, angle, mac))
plotDensitySignalStrength(df)


#### Average Signal Strength Distribution ####

# create offline summary data frame
offlineSummary = createOfflineSummary(offline)

# plot stddev of average signal strength per targeted APs
# result: reinforces that both addresses are important for full range coverage
df = subset(offlineSummary, mac %in% c(chosenAP, rejectedAP), c(sdSignal, avgSignal, mac))
plotStdDevSignalStrength(df)

#### k-NN ####

# The statistical technique used to determine signal detection location is the k-nearest neighbors
# or k-nn. The model will be trained with existing signal plus angle to coordinate position and used to predict
# a new coordinate given a signal strength and angle.
#
# For part 1, we will train the model for the chosen, rejected, and combined APs. The online data
# points will be used to predict their coorindations and the error prediction will be determined
# for each model.

# look at online data
macs = unique(offlineSummary$mac)
macs
online = readData("online.final.trace.txt", subMacs = macs)
online$posXY = paste(online$posX, online$posY, sep="-")
length(unique(online$posXY))
tabonlineXYA = table(online$posXY, online$angle)
# output shows that measurements were taken at all angles throughout the floor
tabonlineXYA

# organize the data where each AP is in a column
onlineSummary = castOnline(online)
# confirm rejected AP is now included
dim(onlineSummary)
names(onlineSummary)

# for k-NN, we will include training data with angles close to point in question since angle matters.
# if we want one angle, then include angles that match the rounded orientation of new observation.
# if we want two angles, then pick two multiples of 45 degrees that flank the new observation's orientation.
# if we want three angles, then pick the closest 45 degree increment and on either side of it.
# what if we didn't want to collapse the signal strengths across the m angles, and
# instead return a set of mx166 signals for each access point?

# iterate through 1 to 3 neighbors with 1 to 3 angles for each scenario to determine best combination
# for cross validation
# result: it is evident that we should include both access points for location predictions.
for (ap in c("None", rejectedAP, chosenAP)){
  if (ap == "None"){
    test = onlineSummary
    train = offlineSummary
    newSignals = test[ , 6:12]
  } else {
    test = onlineSummary[, !(names(onlineSummary) %in% ap)]
    train = subset(offlineSummary, mac != ap)
    newSignals = test[ , 6:11]
  }

  estXY = predXY(newSignals = newSignals,
                 newAngles = test[ , 4],
                 train,
                 numAngles = 3, k = 3)
  actXY = test[ , c("posX", "posY")]
  err = calcError(estXY, actXY)

  print(paste(ap, err))
}

# perform cross validation and plot for desired AP for a range of k's
# for each fold, we want to find k-NN estimates from 1..K and aggregate errros over the v folds.
# choose a v-fold of 11
v = 11
# number of locations
permuteLocs = sample(unique(offlineSummary$posXY))
# matrix locations into folds
permuteLocs = matrix(permuteLocs, ncol = v, nrow = floor(length(permuteLocs)/v))

# summarize and format offline
keepVars = c("posXY", "posX","posY", "orientation", "angle")
onlineCVSummary = reshapeSS(offline, keepVars = keepVars, sampleAngle = TRUE)

# number of k's
K = 10
# errors array
err = rep(0, K)

# loop through folds
for (j in 1:v) {
  testFold = subset(onlineCVSummary, posXY %in% permuteLocs[ , j])
  trainFold = subset(offlineSummary, posXY %in% permuteLocs[ , -j])
  actFold = testFold[ , c("posX", "posY")]

  # loop through neighbors
  for (k in 1:K) {
    estFold = predXY(newSignals = testFold[ , 6:12],
                     newAngles = testFold[ , 4],
                     trainFold, numAngles = 3, k = k)
    err[k] = err[k] + calcError(estFold, actFold)
  }
}
# plot k to sum of squared errors
plotSSErrors(err, K)

# experiment with a few k's based on above analysis to determine k
# that results in lowest sum squared error on test data.
estXY = predXY(newSignals = onlineSummary[ , 6:12],
               newAngles = onlineSummary[ , 4],
               offlineSummary, numAngles = 3, k = 5)
actXY = onlineSummary[ , c("posX", "posY")]
calcError(estXY, actXY)

#### Part 2: Alternatvie k-nearest ####


