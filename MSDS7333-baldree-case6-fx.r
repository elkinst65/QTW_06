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