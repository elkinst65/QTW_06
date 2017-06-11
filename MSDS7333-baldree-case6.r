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
source("MSDS7333-baldree-case6-fx.r")

#### Part 1: Access point MAC address analysis of 00:0f:a3:39:e1:c0 and 00:0f:a3:39:dd:cd ####

# read offline data
offline = readData()

# addresses in question
chosenMac = "00:0f:a3:39:e1:c0"
rejectedMac = "00:0f:a3:39:dd:cd"

# signal map of both mac addresses
# result: map looks equivalent
plotSignalMap(offline[offline$mac %in% c(chosenMac), ])
plotSignalMap(offline[offline$mac %in% c(rejectedMac), ])

# box plot of two mac addresses
# result: shows a tighter spread with better mean signal for chosenMac but a number of outliers
boxplotSS(offline[offline$mac %in% c(chosenMac, rejectedMac), ])

### SANDBOX ####
# our mac addresses have the most entries
subMacs = names(sort(table(offline$mac), decreasing = TRUE))[1:7]
identical(subMacs[1:2], c(chosenMac, rejectedMac))

# show heatmap for both mac addresses

# box plot two macs similar to 1.4.1

# plot signal distribution

# SD of signal strength by mean

### SANDBOX ^^^^^^




# relationship between signal strength and distance

# floor error map

# training with each mac address

