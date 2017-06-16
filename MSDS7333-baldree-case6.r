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
# 2. Implement alternative k-nearest precition method using weights on received signal strength.
#   For what range of values of weights are you able to obtain better prediction values than for the unweighted approach?
#     Use calcError() to compare approaches.
#
# Write up assignment in an Ipython Notebook with an introduction, explanation of approaches and outputs.
#

# two digit precision for printing and formatting
options(digits = 2)
options(error = recover, warn = 2)
# include necessary functions
#if (!exists("processLine", mode="function")) source("MSDS7333-baldree-case6-fx.r")
source("MSDS7333-baldree-case6-fx.r", print.eval=TRUE)

#### Part 1: Access point MAC address analysis of 00:0f:a3:39:e1:c0 and 00:0f:a3:39:dd:cd ####

# read offline data
offline = readData()

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

# In order to examine distribution for all locations, angles, and are two interested APs, we
# will create a summary statistics for all location-orientation-AP combinations with a new factor.
# For each combination there are around 100 observations.
offline$posXY = paste(offline$posX, offline$posY, sep="-")

# create data frames for each combination
byLocAngleAP = with(offline, by(offline, list(posXY, angle, mac), function (x) x))

# summary statistic
signalSummary = lapply(byLocAngleAP, function(oneLoc) {
  ans = oneLoc[1, ]
  ans$medSignal = median(oneLoc$signal)
  ans$avgSignal = mean(oneLoc$signal)
  ans$num = length(oneLoc$signal)
  ans$sdSignal = sd(oneLoc$signal)
  ans$iqrSignal = IQR(oneLoc$signal)
  ans
})
offlineSummary = do.call("rbind", signalSummary)

# plot stddev of average signal strength per targeted APs
# result: reinforces that both addresses are important for full range coverage
df = subset(offlineSummary, mac %in% c(chosenAP, rejectedAP), c(sdSignal, avgSignal, mac))
plotStdDevSignalStrength(df)

#### k-NN ####


