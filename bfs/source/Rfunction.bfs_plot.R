#PLOT HYDROGRAPHS FROM BF_SEP OUTPUT
bfs_plot=function(siteid,bfs_out,qthresh) {

dy=bfs_out$Date
nr=nrow(bfs_out)
wyst=c(1:nr)
wyst=wyst[substr(dy,6,10)=='10-01'] #INDEX FOR START OF WATER YEARS
wyst=wyst[order(wyst,decreasing=TRUE)] #REORDER STARTING WITH LAST PERIOD
if(length(wyst)>1) {plst=wyst[2];plen=plst+729} #START AND END INDICES FOR PLOT
if(length(wyst)<2) {plst=match(TRUE, as.Date(bfs_out$Date)>as.Date(bfs_out$Date[nr])-730);plen=nr}

for(w in wyst) {if(sum(bfs_out$Qob.L3[w:(w+729)],na.rm=TRUE)>0){
    if(sum(bfs_out$Qob.L3[plst:plen],na.rm=TRUE)==0){plst=w; plen=w+729} else {
    if(sum(is.na(bfs_out$Qob.L3[plst:plen]))>60) {if(sum(is.na(bfs_out$Qob.L3[w:(w+729)]))<sum(is.na(bfs_out$Qob.L3[plst:plen]))) {plst=w; plen=w+729}}}}}

ts=as.numeric(substr(bfs_out$Date[plst:plen],6,7))
l=plen-plst+1
tmp=(ts[2:l]==ts[1:(l-1)])
ts[c(FALSE,tmp)]=0
ts_index=c(plst:plen)
ts_index=ts_index[ts>0]
if(length(ts_index)>12) {x=seq(1,length(ts_index),12)
ts_lab_index=ts_index[x]} else {ts_lab_index=ts_index}
ts_lab=substr(dy[ts_lab_index],1,10)

yl=expression(paste('Streamflow [',m^3,day^-1,']'))

layout(matrix(c(1,2), nrow=2, ncol=1))
par(mar=c(2,4,4,4), cex.lab=0.75)

#LINEAR AXIS FOR STREAMFLOW
ylm=c(0,1.25*max(bfs_out$Qob.L3[plst:plen],na.rm=TRUE))

plot(c(plst:plen),bfs_out$Qob.L3[plst:plen], xlim=c(plst,plen), ylim=ylm, type='l', xaxt='n', yaxt='n', xlab=NA, ylab=yl, cex.axis=0.75,lwd=0.75, bty='n')

axis(1,at=ts_index,labels=NA,tck=-0.02,cex.axis=0.75)
axis(1,at=ts_lab_index,labels=ts_lab,tck=-0.03,cex.axis=0.75)

ytck=seq(0,signif(ylm[2],2),signif(ylm[2],2)/4)
axis(2,at=ytck,labels=format(ytck,scientific=TRUE,digits=2),tck=0.03,cex.axis=0.75)
#lines(c(plst:plen),bfs_out$SurfaceFlow.L3[plst:plen]+bfs_out$Baseflow.L3[plst:plen], col=rgb(0.5,0.5,1),lwd=1)
lines(c(plst:plen),bfs_out$Baseflow.L3[plst:plen], col=rgb(0.1,0.6,0.2),lwd=2)

text(0.9*(plen-plst)+plst,ylm[2],'Streamflow (OBSERVED)',pos=2,cex=0.75)
text(0.9*(plen-plst)+plst,0.9*ylm[2],'Baseflow (CALCULATED)',pos=2,cex=0.75)
lines(c(0.9*(plen-plst)+plst,0.95*(plen-plst)+plst),rep(ylm[2],2),col=gray(0))
lines(c(0.9*(plen-plst)+plst,0.95*(plen-plst)+plst),rep(0.9*ylm[2],2),col=rgb(0.1,0.6,0.2),lwd=2)
#text(0.9*(plen-plst)+plst,0.8*ylm[2],'Surface Flow (CALCULATED)',pos=2,cex=0.75)
#lines(c(0.9*(plen-plst)+plst,0.95*(plen-plst)+plst),rep(0.8*ylm[2],2),col=rgb(0.5,0.5,1),lwd=1)


mtext(paste('Station',siteid),3,2)

#LOG AXIS FOR STREAMFLOW
ylm=c(min(c(bfs_out$Qob.L3[bfs_out$Qob.L3>0],0.6*qthresh),na.rm=TRUE),max(bfs_out$Qob.L3[plst:plen],na.rm=TRUE))

plot(c(plst:plen),bfs_out$Qob.L3[plst:plen], xlim=c(plst,plen), ylim=ylm, type='l', xaxt='n',yaxt='n',xlab=NA, ylab=yl,log='y',cex.axis=0.75,lwd=0.75,bty='n')

axis(1,at=ts_index,labels=NA,tck=-0.02,cex.axis=0.75)
axis(1,at=ts_lab_index,labels=ts_lab,tck=-0.03,cex.axis=0.75)
ytck=10^seq(floor(log(ylm[1],10)),ceiling(log(ylm[2],10)))
axis(2,at=ytck,labels=format(ytck,scientific=TRUE,digits=2),tck=0.03,cex.axis=0.75)

ytck_minor=c()
for(y in ytck){ytck_minor=c(ytck_minor,y*c(1:9))}
axis(2,at=ytck_minor,labels=NA,tck=0.02,cex.axis=0.75)

#lines(c(plst:plen),bfs_out$SurfaceFlow.L3[plst:plen]+bfs_out$Baseflow.L3[plst:plen], #col=rgb(0.5,0.5,1),lwd=1)
lines(c(plst:plen),bfs_out$Baseflow.L3[plst:plen], col=rgb(0.1,0.6,0.2),lwd=2)
lines(c(plst,plen),rep(qthresh,2),lty=2,col=gray(0.5))

ycfs=ylm*35.31/24/3600
ytckcfs=10^(seq(floor(log(ycfs[1],10)),ceiling(log(ycfs[2],10))))
axis(4,at=ytckcfs/35.31*24*3600,labels=ytckcfs,tck=0.03,cex.axis=0.75)
ytck_minor=c()
for(y in ytckcfs){ytck_minor=c(ytck_minor,y*c(1:9))}
axis(4,at=ytck_minor/35.31*24*3600,labels=NA,tck=0.02,cex.axis=0.75)

text(0.2*(plen-plst)+plst,qthresh,'Threshold for Calibration',pos=1,cex=0.75)
mtext('Streamflow [cfs]',4,2,cex=0.75)

}
