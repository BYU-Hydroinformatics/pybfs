###########################################################################
# CALIBRATION SCRIPT FOR SITE 12167000
###########################################################################

# Set working directory (adjust if needed)
# setwd('/Users/njones/python_projects/pybfs-1/bfs/calibration_R_original')

# Set CRAN mirror for package installation
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install and load required package
if (!require(quantreg, quietly = TRUE)) {
  install.packages('quantreg', dependencies=TRUE)
}
library(quantreg)

# Source required functions
source('source/Rfunction.bfs.R')
source('source/Rfunctions.bfs_utilities.R')
source('source/Rfunctions.bfs_calibration_sub.R')
source('source/Rfunction.bfs_calibrate.R')
source('source/Rfunction.bfs_plot.R')

# Read site information
site_info <- read.csv('siteinfo_12167000.csv')
SiteID <- as.character(site_info$SiteID[1])
AREA <- as.numeric(site_info$Area.M2[1])

# Read dates and streamflow
dys <- readLines('dates_12167000.csv')
Qin <- as.numeric(readLines('streamflow_12167000.csv'))

# Print information
cat('Running BFS calibration for site:', SiteID, '\n')
cat('Drainage area:', AREA, 'square meters\n')
cat('Number of days:', length(dys), '\n')
cat('Date range:', dys[1], 'to', dys[length(dys)], '\n\n')

# Run calibration
cat('Starting calibration...\n')
out <- bfs_calibrate(SiteID, AREA, Qin, dys)

cat('Calibration complete!\n\n')

# Create output directories if they don't exist
dir.create('out/site_sum', showWarnings = FALSE, recursive = TRUE)
dir.create('out/ci_tables', showWarnings = FALSE, recursive = TRUE)
dir.create('out/timeseries', showWarnings = FALSE, recursive = TRUE)
dir.create('out/hydrographs', showWarnings = FALSE, recursive = TRUE)

# Write output files
cat('Writing output files...\n')
write.csv(out[[1]], paste('out/site_sum/bfs_params_', SiteID, '.csv', sep=''), row.names=FALSE)
write.csv(out[[2]], paste('out/site_sum/bff_', SiteID, '.csv', sep=''), row.names=FALSE)
write.csv(out[[3]], paste('out/ci_tables/ci_table_', SiteID, '.csv', sep=''), row.names=FALSE)
write.csv(out[[4]], paste('out/timeseries/bfs_', SiteID, '.csv', sep=''), row.names=FALSE)

# Create hydrograph
cat('Creating hydrograph...\n')
pdf(file=paste('out/hydrographs/bfs_', SiteID, '.pdf', sep=''), height=7.5, width=10)
bfs_plot(siteid=SiteID, bfs_out=out[[4]], qthresh=out[[1]]$Qthresh)
dev.off()

cat('\nCalibration complete! Output files written to out/ directory.\n')

