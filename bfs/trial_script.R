###########################################################################
#CALIBRATION OF BASEFLOW SEPARATION (BFS) MODEL
###########################################################################
#The BFS model can be calibrated in R using the function bfs_calibrate(SiteID,AREA,Qin,dys). Source code files should be in a sub-directory "source" under the current working directory. The suggested directory structure for calibrating multiple sites is:

#     /[working directory]
#         /source
#         /out
#              /ci_tables
#              /hydrographs
#              /timeseries
#              /site_sum
setwd("C:/Users/ewebsteresho/OneDrive - The University of Alabama/University of Alabama/Spring 2025/new/Trial/calibration_template/calibration_template")

#There are four functions that must be sourced:
  
source('source/Rfunction.bfs.R')
source('source/Rfunctions.bfs_utilities.R')
source('source/Rfunctions.bfs_calibration_sub.R')
source('source/Rfunction.bfs_calibrate.R')

# To use the plotting function (not required for calibration),

source('source/Rfunction.bfs_plot.R')

#BFS calls a function "flow_metrics" that requires the quantreg package:

install.packages(c('quantreg'), dependencies=T)
library(quantreg)

#Arguments for the main calibration function, bfs_calibrate(SiteID, AREA, Qin, dys):

#SiteID: a character string that identifies a site
#AREA: an numeric value of the drainage area for the site with units of square meters
#Qin: a numeric vector of observed streamflow with units of cubic meters per day
#dys: a character vector of days (YYYY-MM-DD)


#bfs_calibrate returns a list with four elements.
#bf_params: a vector with parameters
#bff: a summary of the baseflow, surface flow, and direct runoff fractions
#ci_table: a table of credible intervals (5% to 95%0 for streamflow given baseflow
#bfs_out: the time series output from bfs

#The function can be run as
#out=bfs_calibrate(SiteID,AREA,Qin,dys)
#where out[[1]]=bf_params, out[[2]]=bff, out[[3]]=ci_table, and out[[4]]=bfs_out

##########################################################################################################
#SPECIFY ARGUMENTS
#SiteID = 
#AREA = 
#Qin =
dys= seq.Date(from = as.Date("2013-10-01"), to = as.Date("2015-09-30"), by = "day")
##########################################################################################################
#THESE LINES ARE CALIBRATING MULTIPLE SITES

#SITE FILE FOR CALIBRATING MULTIPLE SITES. MUST HAVE COLUMNS "SiteId" (character) AND "Area.M2" (numeric)
sites=read.csv("siteinfo.csv")

#STREAMFLOW TABLE FOR CALIBRATING MULTIPLE SITES. MUST HAVE UNITS OF M3/DAY#
qdv=read.csv("streamflow.csv")

ns = nrow(site_info) #NUMBER OF SITES

for(s in 1:ns){
  SiteID=sites$SiteID[s]
  AREA=sites$Area.M2[s]
  Qin=qdv[,s]
##########################################################################################################
out=bfs_calibrate(SiteID,AREA,Qin,dys)

#WRITE OUTPUT AS FILES USING THE SUGGESTED DIRECTORY STRUCTURE  
write.csv(out[[1]],paste('out/site_sum/bfs_params_',SiteID,'.csv',sep=''),row.names=F)
write.csv(out[[2]],paste('out/site_sum/bff_',SiteID,'.csv',sep=''),row.names=F)
write.csv(out[[3]],paste('out/ci_tables/ci_table_',SiteID,'.csv',sep=''),row.names=F)
write.csv(out[[4]],paste('out/timeseries/bfs_',SiteID,'.csv',sep=''),row.names=F)

#CREATE A HYDROGRAPH
pdf(file=paste('out/hydrographs/bfs_',SiteID,'.pdf',sep=''), height=7.5, width=10)  
bfs_plot(siteid=SiteID,bfs_out=out[[4]],qthresh=out[[1]]$Qthresh)
dev.off()
} #CLOSE MULTI-SITE CALIBRATION FOR LOOP
###########################################################################
