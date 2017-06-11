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


#######################################################################################
# Part 1: Access point MAC address analysis of 00:0f:a3:39:e1:c0 and 00:0f:a3:39:dd:cd
#######################################################################################

# two digit precision for printing and formatting
options(digits = 2)


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

# read offline data
offline = readData()
