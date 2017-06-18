processLine = function(x)
{
  # Process a line of text returning a 11x10 matrix of tokens
  #
  # Args:
  #   x: The line of text
  #
  tokens = strsplit(x, "[;=,]")[[1]]

  # return NULL if only 10 tokens
  if (length(tokens) == 10)
    return(NULL)

  tmp = matrix(tokens[-(1:10)], ncol = 4, byrow = TRUE)
  mat = cbind(matrix(tokens[c(2, 4, 6:8, 10)], nrow(tmp), 6, byrow = TRUE), tmp)
}

roundOrientation = function(angles) {
  # Round an array of degree angles to the nearest 45 degree angle.
  #
  # Args:
  #   angles: The line of text
  #
  refs = seq(0, by = 45, length  = 9)
  q = sapply(angles, function(o) which.min(abs(o - refs)))
  c(refs[1:8], 0)[q]
}


readData = function(filename = 'offline.final.trace.txt',
                    subMacs = c("00:0f:a3:39:e1:c0", "00:0f:a3:39:dd:cd", "00:14:bf:b1:97:8a",
                                "00:14:bf:3b:c7:c6", "00:14:bf:b1:97:90", "00:14:bf:b1:97:8d",
                                "00:14:bf:b1:97:81")) {
  # Process the csv file peforming the following steps:
  #   Read all non-comment lines.
  #   Call processLine() for each line of text.
  #   Convert list of text into a data frame
  #   Name columns of data frame
  #   Define column types
  #   Convert millisecond time to seconds
  #   Round angle orientations
  #   Drop data
  #
  # Args:
  #   filename: The name of the csv file to process. The default is 'offline.final.trace.txt.'
  #   subMacs: The list of access point MAC addresses to filter the dataset down to. A default list is provided.
  #
  # Returns:
  #   The data frame
  #

  # process csv into matrix of tokens
  txt = readLines(filename)
  lines = txt[substr(txt, 1, 1) != "#"]
  tmp = lapply(lines, processLine)
  # convert into tokens into a data frame
  df = as.data.frame(do.call("rbind", tmp), stringsAsFactors = FALSE)

  # variable names
  names(df) = c("time", "scanMac", "posX", "posY", "posZ", "orientation", "mac", "signal", "channel", "type")

  # convert numeric values
  numVars = c("time", "posX", "posY", "orientation", "signal")
  df[numVars] = lapply(df[numVars], as.numeric)

  # convert time to POSIX
  df$rawTime = df$time
  df$time = df$time / 1000
  class(df$time) = c("POSIXt", "POSIXct")

  # round orientations to nearest 45
  df$angle = roundOrientation(df$orientation)

  # keep only signals from access points
  df = df[df$type == "3",]
  df = df[, "type" != names(df)]

  # drop scanMac, posZ, channel, and type - no info in them
  df = df[, !(names(df) %in% c("scanMac", "posZ"))]

  # drop more unwanted access points
  df = df[df$mac %in% subMacs,]

  # drop channel
  df = df[, "channel" != names(df)]

  return(df)
}

plotBoxplotSignalStrength = function(df) {
  # Box plot signal strength
  #
  # Args:
  #   df: data frame of data to analyze
  #
  oldPar = par(mar = c(3.1, 3.1, 1, 1), mfrow = c(1,1))
  par(mai=c(.8, .8, .5, .25))
  library (lattice)

  print(bwplot(signal ~ factor(angle) | mac, data = df, layout = c(2,1),
               main="Signal Strength Distribution for Access Points",
               xlab="Angle Measured (deg)", ylab="Signal Strengh (dBM)"))
  par(oldPar)
}

plotDensitySignalStrength = function(df) {
  # Density plot signal strength
  #
  # Args:
  #   df: data frame of data to analyze
  #
  oldPar = par(mar = c(3.1, 3.1, 1, 1), mfrow = c(1,1))
  par(mai=c(.8, .8, .5, .25))
  library (lattice)

  print(densityplot( ~ signal | mac + factor(angle), data = df, bw = 0.5, plot.points = FALSE,
        main="Signal Strength Distribution for Access Points",
        xlab="Angle Measured (deg)", ylab="Signal Strengh (dBM)"))
  par(oldPar)
}

plotStdDevSignalStrength = function(df) {
  # Standard Deviation of signal strength average
  #
  # Args:
  #   df: data frame of data to analyze
  #
  oldPar = par(mar = c(3.1, 3.1, 1, 1), mfrow = c(1,1))
  par(mai=c(.8, .8, .5, .25))
  library (lattice)
  breaks = seq(-90, -30, by=5)

  print(bwplot(sdSignal ~ cut(avgSignal, breaks=breaks) | mac, data=df,
               xlab = "Mean Signal (dBM)", ylab = "Signal Standard Deviation (dBM)",
               main = "Standard Deviation of Average Signal Strength from Detector to Access Point",
               cex.main=.8, cex.axis=.8, cex.lab=.8, layout = c(1, 2)))
  par(oldPar)
}


plotSignalMaps = function(macAddresses) {
  # Plot the number of signals per XY location at XY location.
  #
  # Args:
  #   macAddresses: a list of macs to plot signals
  #
  oldPar = par(mar = c(3.1, 3.1, 1, 1), mfrow = c(1,1))
  par(mfrow = c(2, 1), mai=c(.8, .8, .4, .25))

  for (i in macAddresses) {
    # filter data set
    df = offline[offline$mac %in% i, ]
    # location of datapoints filtered
    locDF = with(df, by(df, list(posX, posY), function(x) x))
    # drop locations that were not observed, null
    locDF = locDF[!sapply(locDF, is.null)]

    # determine the number of observations recorded at each location
    locCounts = sapply(locDF, nrow)

    # keep position information with location
    locCounts = sapply(locDF, function(df) c(df[1, c("posX", "posY")], count = nrow(df)))

    # transpose matrix
    locCounts = t(locCounts)
    plot(locCounts, type = "n", xlab = "Position X", ylab = "Position Y",
         main = paste("Total Signals from Detector to Access Point: ", i),
         cex.main=.8, cex.axis=.8, cex.lab=.8)
    text(locCounts, labels = locCounts[, 3], cex = .6, srt = 45)
  }
  par(oldPar)
}

castOnline = function(df) {
  # Cast dataframe of online data so each access point (AP) is in its own column.
  # Matrix is now a 1x7 since we are including rejected AP.
  #
  # Args:
  #   df: online data frame
  #
  keepVars = c("posXY", "posX","posY", "orientation", "angle")
  numOfAPs = length(unique(df$mac))
  byLoc = with(df,
               by(df, list(posXY),
                  function(x) {
                    ans = x[1, keepVars]
                    avgSS = tapply(x$signal, x$mac, mean)
                    y = matrix(avgSS, nrow = 1, ncol = numOfAPs,
                               dimnames = list(ans$posXY, names(avgSS)))
                    cbind(ans, y)
                  }))
  return(do.call("rbind", byLoc))
}

createOfflineSummary = function(df) {
  # In order to examine distribution for all locations, angles, and are two interested APs, we
  # will create a summary statistics for all location-orientation-AP combinations with a new factor.
  # For each combination there are around 100 observations.
  #
  # Args:
  #   df: offline data frame
  #
  df$posXY = paste(df$posX, df$posY, sep="-")

  # create data frames for each combination
  byLocAngleAP = with(df, by(df, list(posXY, angle, mac), function (x) x))

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
  return(do.call("rbind", signalSummary))
}

reshapeSS = function(df,
                     varSignal = "signal",
                     keepVars = c("posXY", "posX", "posY"),
                     sampleAngle = FALSE,
                     refs = seq(0, 315, by=45)) {
  # reshape signal strength help function
  # aggregate signal strengths from these angles and create a data structure similar to onlineSummary
  # Args:
  #   df: dataframe to reshape
  #   varSignal: variable to aggregate
  #   keepVars: variables to retain
  #   sampleAngle: true if we want to select one angle at random from df
  #   refs: angle references
  #
  numOfAPs = length(unique(df$mac))
  byLocation = with(df, by(df, list(posXY),
                             function(x) {
                               # select one angle at random for each location
                               if (sampleAngle) {
                                 x = x[x$angle == sample(refs, size = 1), ]}
                               ans = x[1, keepVars]
                               avgSS = tapply(x[, varSignal], x$mac, mean)
                               y = matrix(avgSS, nrow = 1, ncol = numOfAPs,
                                          dimnames = list(ans$posXY, names(avgSS)))
                               cbind(ans, y)
                             }))
  return(do.call("rbind", byLocation))
}

selectTraininingData = function(newObservationAngle, df = NULL, m = 1){
  # select observations from df to analyze aggregate signal strengths
  # from these angles and create a data structure similar to onlineSummary.
  #
  # Args:
  #   newObservationAngle: the angle of the new observation
  #   df: offline summary data frame
  #   m: number, between 1 and 5, of angles to keep

  # m is the number of angles to keep between 1 and 5
  refs = seq(0, by = 45, length  = 8)
  nearestAngle = roundOrientation(newObservationAngle)

  # handle odd and even m number
  if (m %% 2 == 1)
    angles = seq(-45 * (m - 1) /2, 45 * (m - 1) /2, length = m)
  else {
    m = m + 1
    angles = seq(-45 * (m - 1) /2, 45 * (m - 1) /2, length = m)
    if (sign(newObservationAngle - nearestAngle) > -1)
      angles = angles[ -1 ]
    else
      angles = angles[ -m ]
  }
  # map angles to values in refs
  # negative angles and angles greater than 360 are mapped to appropriate angles; e.g., -45 maps to 335 and 405 maps to 45
  angles = angles + nearestAngle
  angles[angles < 0] = angles[ angles < 0 ] + 360
  angles[angles > 360] = angles[ angles > 360 ] - 360
  angles = sort(angles)
  # select observations to analyze
  subset = df[df$angle %in% angles, ]
  # reshape signal strength
  return(reshapeSS(subset, varSignal = "avgSignal"))
}

findNN = function(newSignal, trainingSubset) {
  # want to look the distance in terms of signal strengths from these training data to the new
  # data point. we need to calculate teh distrance from the new point to all observations in the
  # training set with findNN().
  # returns locations of the training observations in order of closeness to the new observation's signal strength.
  #
  # Args:
  #   newSignal: signal of new observation
  #   trainingSubset: training data to find neighbors
  # Returns:
  #   training neighbors
  #
  cols = length(trainingSubset)
  diffs = apply(trainingSubset[ , 4:cols], 1, function(x) x - newSignal)
  dists = apply(diffs, 2, function(x) sqrt(sum(x^2)) )
  closest = order(dists)
  return(trainingSubset[closest, 1:3 ])
}

predXY = function(newSignals, newAngles, training, numAngles = 1, k = 3){
  # Predict the XY coordinates given a list of signals along with their angles measured and training data.
  # k neighbors will be found to make the prediction of coordinate.
  #
  # Args:
  #   newSignals: list of signals to predict
  #   newAngles: list of angles for signals to predict
  #   trainingData: data needed for k-nn
  #   numAngles: angles we want to use for prediction
  #   k: number of neighbors to use for prediction
  #
  closeXY = list(length = nrow(newSignals))

  for (i in 1:nrow(newSignals)) {
    trainSS = selectTraininingData(newAngles[i], training, m = numAngles)
    closeXY[[i]] = findNN(newSignal = as.numeric(newSignals[i, ]), trainSS)
  }
  # for some k of nearest neighbors, simply average the first k locations
  # could have used weights in the average that are inversely proportional to the distance in signal strength
  # from the test observation. will need to call findNN() to get distance to point. the weights
  # are 1/d / sum(1/d).
  # could also use a different metric besides Euclidean like Manhattan.
  # could use medians instead of averages when combining neighbors if the distribution of values are quite skewed.
  estXY = lapply(closeXY, function(x) sapply(x[ , 2:3], function(x) mean(x[1:k])))

  return(do.call("rbind", estXY))
}

calcError = function(estXY, actualXY) {
  # compare fit numerically with sum of squared errors
  #
  # Args:
  #   estXY:
  #   actualXY:
  # Returns:
  #   calculated residual sum of squares
  return(sum(rowSums((estXY - actualXY) ^ 2)))
}

plotSSErrors = function(err, K){
  # Plot results as sum of squared errors as function of k
  #
  # Args:
  #   err: Error array
  #   K: number of neighbors
  #

  oldPar = par(mar = c(3.1, 3.1, 1, 1), mfrow = c(1,1))
  par(mai=c(.8, .8, .4, .25))

  plot(y = err, x = (1:K),  type = "l", lwd= 2, ylim = c(900, 2100), cex.main=.8, cex.axis=.8, cex.lab=.8,
       xlab = "Number of Neighbors",
       ylab = "Sum of Squared Errors",
       main = "Cross Validation Section of k")

  rmseMin = min(err)
  kMin = which(err == rmseMin)[1]
  segments(x0 = 0, x1 = kMin, y0 = rmseMin, col = gray(0.4), lty = 2, lwd = 2)
  segments(x0 = kMin, x1 = kMin, y0 = 800,  y1 = rmseMin, col = grey(0.4), lty = 2, lwd = 2)
  text(x = kMin - 2, y = rmseMin + 40, label = as.character(round(rmseMin)), col = grey(0.4), cex = .6)

  par(oldPar)
}

