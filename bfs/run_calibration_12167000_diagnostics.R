# R script to run BFS calibration with detailed diagnostics
# This matches the Python diagnostic output for comparison

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install required packages if needed
if (!require(quantreg, quietly = TRUE)) {
  install.packages('quantreg', dependencies=TRUE)
}
library(quantreg)

# Source all required functions (matching run_calibration_12167000.R)
source('source/Rfunction.bfs.R')
source('source/Rfunctions.bfs_utilities.R')
source('source/Rfunctions.bfs_calibration_sub.R')
# Use diagnostic version instead of regular version
source('source/Rfunction.bfs_calibrate_diagnostics.R')

# Read input data
siteinfo <- read.csv('siteinfo_12167000.csv')
dates <- read.csv('dates_12167000.csv', header=FALSE, col.names='Date')
streamflow <- read.csv('streamflow_12167000.csv', header=FALSE, col.names='Streamflow')

# Extract values
tmp.site <- as.character(siteinfo$SiteID[1])
# Handle both AREA and Area.M2 column names
if ('AREA' %in% names(siteinfo)) {
  tmp.area <- as.numeric(siteinfo$AREA[1])
} else if ('Area.M2' %in% names(siteinfo)) {
  tmp.area <- as.numeric(siteinfo$Area.M2[1])
} else {
  stop("Could not find AREA or Area.M2 column in siteinfo")
}
tmp.q <- as.numeric(streamflow$Streamflow)
dys <- as.Date(dates$Date)

# Create output directories
dir.create('out', showWarnings = FALSE)
dir.create('out/site_sum', showWarnings = FALSE, recursive = TRUE)
dir.create('out/ci_tables', showWarnings = FALSE, recursive = TRUE)
dir.create('out/timeseries', showWarnings = FALSE, recursive = TRUE)
dir.create('out/hydrographs', showWarnings = FALSE, recursive = TRUE)

cat("\n============================================================\n")
cat("Running R calibration with diagnostics...\n")
cat("============================================================\n\n")

# Run calibration with diagnostics
result <- bfs_calibrate_diagnostics(tmp.site, tmp.area, tmp.q, dys)

# The diagnostics are printed within the modified bfs_calibrate function
# For now, we'll create a modified version that prints diagnostics

cat("\n============================================================\n")
cat("R CALIBRATION COMPLETE\n")
cat("============================================================\n\n")

