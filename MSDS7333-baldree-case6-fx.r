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
