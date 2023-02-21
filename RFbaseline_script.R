######################################################################
#                                                                    #
#                 R Script to calculate the baseline of              #
#                 atmospheric data using RFbaseline                  #
#                                                                    #
#                 Created by P.A.Pickers on 20Jan2014                #
#                                                                    #
######################################################################

# This script explains how to use RFbaseline (also called REBS - robust
# extraction of baseline signal). RFbaseline is a tool for fitting a
# baseline to atmospheric mole fraction data. Unlike some other 
# baseline fitting routines, RFbaseline uses an assymetric distribution
# of the residuals of the LOESS fit, which means that the baseline is
# not biased by repeated, large pollution spikes that occur in the
# same direction. 

# Before running this script, you must download the IDPmisc package 
# from the CRAN repository (if you haven't already).

# For more information about RFbaseline, please see Ruckstuhl et al. 
# (2001) and Ruckstuhl et al. (2012).

# First, check the working directory:

getwd()

# To change you working directory, use the following command:


# Note: R doesn't recognise '\', so they must be replaced with '/'.

# Now load the IDPmisc package (note: skip this step if the package 
# is already loaded):

library("IDPmisc")
library(ncdf4)

# Next, import the data. For this script, the input data file 
# should be in text format with a date column as the first column 
# and the mole fraction data for the second column.

file_dir <- file.path("/Users", "vf20487", "Documents", "Data", "Timeseries") #, fsep="\\")
setwd(file_dir)

apo_ts <- nc_open("WAO_APO_timeseries_2015.nc")
time <- array(ncvar_get(apo_ts, 'time'))

# make a data frame to add the fits into
RFoutput <- data.frame(matrix(ncol = 5, nrow = length(time)))
RFoutput$X1 <- time
# rename the data frame columns
names(RFoutput) <- c("time", 'no_ocean', 'ecco', 'jena', 'nemo')

# loop through the ocean models
oceans = c('no_ocean', 'ecco', 'jena', 'nemo')
for (ocean in oceans) {
  # get the APO model
  apo <- ncvar_get(apo_ts, ocean)
  
  if (ocean=='nemo') {b=3.0} else if (ocean=='ecco') {b=2.0} else {b=1.0}
  
  # run the fit
  RF <- rfbaseline(time,apo*-1,span=0.08,maxit=c(5,0), b = b)
  # add the fit into the data frame
  RFoutput[ocean] <- RF$fit * -1
  
  # plot the model and fit 
  plot((apo) ~ time)
  lines(time,RF$fit*-1,col="purple",lwd=3)
}

# The important variable is the span, which should be between 0 and
# 1, with higher values giving a smoother baseline. Also, the 'b' term affects
# how well the baseline fits along the bottom of the data. 

# For more imformation about how to implement the routine, please
# refer to the IDPmisc package help pdf.

# save to csv
write.csv(RFoutput,file="WAO_APOmodel_REBSbaseline_2015.csv", row.names = FALSE)

#### END OF SCRIPT ####
