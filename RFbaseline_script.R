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

apo_ts <- nc_open("apo_timeseries_2015.nc")

apo_ecco <- ncvar_get(apo_ts, 'ecco')
apo_jena <- ncvar_get(apo_ts, 'jena')
apo_nemo <- ncvar_get(apo_ts, 'nemo')
time <- ncvar_get(apo_ts, 'time')


#inputdata <- read.csv("WAO_May2010-Jan2021_interpolated_APO.txt", header = TRUE, sep = "\t")
#inputdata$date <- as.POSIXct(strptime(inputdata$date, format = "%d/%m/%Y %H:%M", tz = "GMT"))
#inputdata.2 <- inputdata[,c(1,2)]
#inputdata.2 <- na.omit(inputdata.2)

# Now split the input data file into the two variables:

date <- array(time) # inputdata.2[,1]
molefract <- apo_nemo # inputdata.2[,2] *-1

# Now implement the RFbaseline routine:

RF <- rfbaseline(date,molefract,span=0.12,maxit=c(5,0), b = 0.6)

# The important variable is the span, which should be between 0 and
# 1, with higher values giving a smoother baseline. Also, the 'b' term affects
# how well the baseline fits along the bottom of the data. 

# For more imformation about how to implement the routine, please
# refer to the IDPmisc package help pdf.

# To check whether the fit is correct, plot the results:

#plot((inputdata.2$APO) ~ inputdata.2$date)
plot((molefract) ~ date)
lines(date,-RF$fit*-1,col="purple",lwd=3)

# Now export the data as a csv file with two columns - the date 
# and the RF fit.

RFoutput <- data.frame(matrix(ncol = 2, nrow = length(RF$fit)))
RFoutput$X1 <- date
RFoutput$X2 <- RF$fit
names(RFoutput) <- c("date", "APO.baseline")
RFoutput$APO.baseline <- RFoutput$APO.baseline *-1
RFoutput.2 <- as.data.frame(timeAverage(RFoutput, avg.time = "hour", fill = TRUE))
#RFoutput <- cbind(date,RF$fit)
write.csv(RFoutput,file="WAO_model_REBSbaseline_2015.csv", row.names = FALSE)

#### END OF SCRIPT ####
