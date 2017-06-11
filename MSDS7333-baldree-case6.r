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
# include necessary functions
if (!exists("processLine", mode="function")) source("MSDS7333-baldree-case6-fx.r")

#### Part 1: Access point MAC address analysis of 00:0f:a3:39:e1:c0 and 00:0f:a3:39:dd:cd ####

# read offline data
offline = readData()

# signal strength analysis of two mac addresses

