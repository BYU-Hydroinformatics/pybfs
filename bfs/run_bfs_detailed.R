# Run BFS and save detailed output for comparison with Python

library(quantreg)
source('source/Rfunctions.bfs_utilities.R')
source('source/Rfunction.bfs.R')

# Load parameters
bfs_params <- read.csv('out/site_sum/bfs_params_12167000.csv')

# Load streamflow data
streamflow_data <- read.csv('12167000.csv', fileEncoding='UTF-8-BOM')
streamflow_data$Date <- as.Date(streamflow_data$Date, format='%m/%d/%Y')

tmp.q <- streamflow_data$mean_daily_streamflow
dys <- as.character(streamflow_data$Date)

# Get flow metrics
flow <- flow_metrics(tmp.q, timestep='day', fr4rise=0.05)

# Parameters
tmp.area <- bfs_params$tmp.area[1]
Lb <- bfs_params$Lb[1]
X1 <- bfs_params$X1[1]
Wb <- bfs_params$Wb[1]
POR <- bfs_params$POR[1]
ALPHA <- bfs_params$ALPHA[1]
BETA <- bfs_params$BETA[1]
Ks <- bfs_params$Ks[1]
Kb <- bfs_params$Kb[1]
Kz <- bfs_params$Kz[1]

basin_char <- c(tmp.area, Lb, X1, Wb, POR)
gw_hyd <- c(ALPHA, BETA, Ks, Kb, Kz)

# Run BFS
Error <- bfs(tmp.q, dys, timestep='day', error_basis='total', basin_char, gw_hyd, flow)

# Get the output from global environment
bfs_out <- get('bfs_out', envir=.GlobalEnv)

# Print first 20 time steps
cat("\nFirst 20 time steps - R BFS results:\n")
cat(sprintf("%-5s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s\n", 
    "ts", "Qob", "Qsim", "Baseflow", "SurfaceFlow", "DirectRunoff", "Eta", "StSur", "StBase", "I"))
cat(paste(rep("-", 150), collapse=""), "\n")
for(i in 1:min(20, nrow(bfs_out))) {
    cat(sprintf("%-5d %-15.2f %-15.2f %-15.2f %-15.2f %-15.2f %-15.2f %-15.2f %-15.2f %-15.6f\n",
        i, 
        bfs_out$Qob.L3[i],
        bfs_out$Qsim.L3[i],
        bfs_out$Baseflow.L3[i],
        bfs_out$SurfaceFlow.L3[i],
        bfs_out$DirectRunoff.L3[i],
        bfs_out$Eta.L3[i],
        bfs_out$StSur.L3[i],
        bfs_out$StBase.L3[i],
        bfs_out$Impulse.L[i]))
}

# Save to CSV
write.csv(bfs_out, 'r_bfs_detailed.csv', row.names=FALSE)
cat("\nSaved R BFS results to r_bfs_detailed.csv\n")

