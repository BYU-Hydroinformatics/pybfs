##############################################################################################################
#THIS SCRIPT DEFINES FIVE FUNCTIONS FOR CALIBRATING THE BASEFLOW SEPARATION MODEL

#THEY ARE INTENDED TO BE USED IN ORDER
#ini_params: SETS INITIAL VALUES OF Lb, Wb, AND Kb TO MATCH RECESSION RATE AND ALLOW BASEFLOW TO EQUAL MEAN FLOW WHILE Kb HAS PHYSICALLY POSSIBLE VALLUE
#cal_initial: INITIAL CALIBRATION OF Lb, X1, Wb, ALPHA, Ks, Kb, Kz
#cal_basetable: ASSIGN X1, Kb, AND B FOR BASE RESERVOIR
#cal_base: CALIBRATE X1,Wb,Kb,Kz
#cal_surface: CALIBRATE Wb,ALPHA,Ks

#CONDITIONS PREVENT BF_SEP FROM BEING CALLED AND, INSTEAD, RETURN NOMINAL VALUE OF 100 FOR OBJECTIVE FUNCTION
###################################################################################
#OBJECTIVE FUNCTION FOR cal_base AND cal_surface
###################################################################################
objective=function(bfs_out,prec) {

APE=bfs_out$AdjPctEr
Weight=bfs_out$Weight
Weight[APE>prec/(bfs_out$Qob.L3+prec)]=0 #DAYS WHEN ERROR IS POSITIVE (OVER PREDICTION) DO NOT COUNT TOWARD MINIMIZING THE OBJECTIVE FUNCTION
Weight[bfs_out$Qob.L3==0 & bfs_out$Baseflow.L3==0]=0 #LIMITS THE INFLUENCE OF HYDROPERIOD FOR NON-PERENNIAL STREAMS
OBJ=sum(Weight*(-1+APE^2),na.rm=T)
OBJ}
##############################################################################################################
ini_params=function(area,lb,x1,wb,por,beta,rb1,tmp.q) {
#ASSIGN INITIAL VALUES OF Lb, Wb, AND Kb
#Kb=Rb1*POR*X1^BETA*(Lb-Xb/2)/BETA/(Xb^(BETA-1)) #HORIZONTAL CONDUCTIVITY OF BASE [L/T] TO MATCH Rb1 WHEN Xb=Lb/2
#Kb=-Rb1*POR*(3/4*Lb)*X1 #HORIZONTAL CONDUCTIVITY OF BASE [L/T] TO MATCH Rb1 WHEN Xb=Lb/2

xb=seq(0,lb,lb/1000)
z=(xb/x1)^beta
dzdx=beta/(x1^beta)*xb^(beta-1)
sb=wb*por*(lb-xb/2)*(xb/x1)^beta

#SET Kb USING Qmean
qmean=mean(tmp.q,na.rm=T)
kb=qmean/(wb*z*dzdx)

#SET Kb USING Rb1
#q/sb=-rb1=(kb*beta*xb^(beta-1))/(por *(lb-xb/2) * x1^beta) 

#LIMIT Kb FROM 1e-7 TO 1e4,
#Kb=-Rb1*(POR*(Lb-Xb/2)*X1^BETA)/(BETA*Xb^(BETA-1))

#BOTH STATEMENTS SHOULD BE TRUE FOR ALL XB 
iterate=TRUE
while(iterate){iterate=FALSE;xb=seq(0,lb,lb/1000)
if(any(-rb1*(por*(lb-xb/2)*x1^beta)/(beta*xb^(beta-1)) < 1e-7)) {lb=lb*1.1;iterate=TRUE} #Kb < 1e-7
if(any(-rb1*(por*(lb-xb/2)*x1^beta)/(beta*xb^(beta-1)) > 1e4)) {lb=lb*0.9;iterate=TRUE}} #Kb > 1e4

#ADJUST wb SO THAT Xb IS BETWEEN 0.1 AND 0.9 X Lb WHEN BASEFLOW IS EQUAL TO MEAN STREAMFLOW AT Xb=Lb/2
xb=seq(0,lb,lb/1000)
z=(xb/x1)^beta
dzdx=beta/(x1^beta)*xb^(beta-1)
sb=wb*por*(lb-xb/2)*(xb/x1)^beta
kb=-rb1*(por*(lb-xb/2)*x1^beta)/(beta*xb^(beta-1))

#BASEFLOW VECTOR CORRESPONDING TO Z VECTOR
iterate=TRUE

while(iterate){iterate=FALSE
  qb=wb*z*kb*dzdx
  qmean.i=match(TRUE,qb>qmean)
  if(is.na(qmean.i)){qmean.i=length(qb)}
  if(xb[qmean.i]<(0.1*lb)){wb=0.9*wb #REDUCE Wb SO Xb IS FURTHER FROM OUTLET
    if(wb>10){iterate=TRUE}}
  if(xb[qmean.i]>(0.9*lb)){wb=1.1*wb #INCREASE Wb SO Xb IS CLOSER TO OUTLET
    if(wb*area>lb){lb=0.9*area/wb; iterate=TRUE}}}

kb=qb[qmean.i]/(wb*z[qmean.i]*dzdx[qmean.i])

data.frame('Lb'=lb,'Wb'=wb,'Kb'=kb)}
##############################################################################################################
cal_initial=function(logx,tmp.q,dys,timestep,error_basis,basin_char,gw_hyd,flow,qmean) {
#INITIAL CALIBRATION
#LOGX=log(c(Lb,Wb,ALPHA,Ks,Kb,Kz),10)

#RETURN 100 FOR BAD PARAMETER VALUES: WIDTH X LENGTH > AREA, SURFACE SLOPE > 10%, HYDRAULIC CONDUCTIVITIES <10^-8 or > 10^5 M/DAY
bad_X=c(10^(logx[1]+logx[2])>basin_char[1],logx[3]>(-1),logx[4]<(-8),(logx[4]>5),logx[5]<(-8),logx[5]>5,logx[6]<(-8),logx[6]>5)

if(any(bad_X)){obj=100} else {
  lb=10^logx[1];wb=10^logx[2]
  basin_char[2]=lb;basin_char[4]=wb

  a=10^logx[3];ks=10^logx[4];kb=10^logx[5];kz=10^logx[6]
  gw_hyd[1]=a;gw_hyd[3]=ks;gw_hyd[4]=kb;gw_hyd[5]=kz

  out=bfs(tmp.q,dys,timestep,error_basis,basin_char,gw_hyd,flow)
  obj=objective(bfs_out,prec=flow[5])}
  
obj}
##############################################################################################################
cal_basetable=function(x,b,params,tmp.q){
#CALIBRATE Lb, X1 AND Kb FOR NON-LINEAR SURFACE OF BASE RESERVOIR
#X=c(Lb,X1,Wb,Kb)
#params=c(AREA,POR,Qmean,Qthresh,RbI,RbS)
area=params[1];por=params[2];qmean=params[3];qthresh=params[4];rbi=params[5];rbs=params[6]

bad_X=c(is.infinite(x),x[1]<=0,x[1]>area/x[3],x[2]<=0,x[3]>area/x[1],x[4]<10^-8,x[4]>10^5)

#RETURN 1000 IF Lb X Wb > AREA
if(any(bad_X)) {obj=1000} else {
  lb=x[1];x1=x[2];wb=x[3];kb=x[4]
  sbt=base_table(lb,x1,wb,b,kb,tmp.q,por)

#CHECK FOR VIABLE SOLUTION, RETURN 1000 IF NOT VIABLE
  if(max(sbt$Q,na.rm=T)<qthresh){obj=1000} else {

#RECESSION COEFFICIENT FOR EACH INTERVAL
    r=sbt$Q/sbt$S
    r[is.na(r)]=0
    r[r==Inf]=0

#INDEX FOR DECILES OF RANGE BETWEEN QMIN AND QMEAN
    tmp.qq=data.frame(10^(seq(log(min(tmp.q[tmp.q>0],na.rm=T)/10,10),log(qmean,10),(log(qmean,10)-log(min(tmp.q[tmp.q>0],na.rm=T)/10,10))/100)))
    m=function(x,q) {sum(q<=x)}
    tmp=apply(tmp.qq,1,m,q=sbt$Q)
    tmp=tmp[sbt$Q[tmp]!=0]
    tmp=tmp[rbi+rbs*log(sbt$Q[tmp],10)<0]
    tmp=unique(tmp)

#OBJECTIVE IS MATCHING RECESSION RATES FOR DECILES FROM QTHRESH TO QMEAN
    if(length(tmp)<10){obj=100} else {#RETURN 1000 IF LESS THAN 10 VALUES IN RANGE
      obj=mean(abs((r[tmp]+(rbi+rbs*log(sbt$Q[tmp],10))))) 
      if(is.infinite(obj)) {obj=100}}}}
      
  obj}
###############################################################################
cal_base=function(x,tmp.q,dys,timestep,error_basis,basin_char,gw_hyd,flow){
#CALIBRATE BASIN LENGTH, BASE WIDTH, BASE CONDUCTIVITY, AND RECHARGE
#X=c(Lb,Wb,Kb,Kz)

#RETURN 1000 IF CONDUCTIVITIES ARE UNREALISTIC
bad_X=c(x[1]*x[2]>basin_char[1],x[1]<0,x[2]<0,x[3]<(10-8),x[3]>10^5,x[4]<10^-8,x[4]>10^5)

if(any(bad_X)){obj=100} else {
  lb=x[1];wb=x[2]
  basin_char[2]=lb;basin_char[4]=wb

  kb=x[3];kz=x[4]
  gw_hyd[4]=kb;gw_hyd[5]=kz

  out=bfs(tmp.q,dys,timestep,error_basis,basin_char,gw_hyd,flow)
  obj=objective(bfs_out,prec=flow[5])}
  
obj}
###############################################################################
cal_surface=function(x,tmp.q,dys,timestep,error_basis,basin_char,gw_hyd,flow){
#CALIBRATE SURFACE
#X=c(Wb,ALPHA,Ks,Kz)
bad_X=c(10^x[1]>basin_char[1]/basin_char[2],x[2]>(-1),x[3]<(-8),x[3]>5,x[4]<(-8),x[4]>5)

if(any(bad_X)){obj=100} else {
  wb=10^x[1]
  basin_char[4]=wb

  a=10^x[2];ks=10^x[3];kz=10^x[4]
  gw_hyd[1]=a;gw_hyd[3]=ks;gw_hyd[5]=kz

  out=bfs(tmp.q,dys,timestep,error_basis,basin_char,gw_hyd,flow)
  obj=objective(bfs_out,prec=flow[5])}
  
obj}
###############################################################################################
