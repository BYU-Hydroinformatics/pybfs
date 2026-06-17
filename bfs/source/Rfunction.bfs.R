bfs=function(qin,dys,timestep,error_basis,basin_char,gw_hyd,flow){

#HYDROGRAPH IS SEPARATED INTO THREE COMPONENTS:
#1. DIRECT RUNOFF OF PPT IMPULSE ON SATURATED SURFACE (SATURATED OVERLAND FLOW) AND PPT-INFILTRATION (INFILTRATION LIMITED RUNOFF)
#2. SUBSURFACE FLOW LATERALLY THROUGH SURFACE ZONE, UNIFORM SATURATED THICKNESS AND HYDRAULIC GRADIENT
#3. BASE FLOW LONGITUDINALLY THROUGH SUBSURFACE ZONE, VARIABLE SATURATED THICKNESS AND HYDRAULIC GRADIENT DEPENDING ON WATER LEVEL

#ARGUMENTS - CONSISTENT LENGTH UNITS FOR LENGTH (E.G., METERS) FOR ALL VARIABLES AND PARAMETERS
#qin is a numeric vector of streamflow volume for each time step
#dys is a vector of dates (YYYY-MM-DD) or POSIX date-time for each time step
#timestep is 'day' or 'hour'
#error_basis is 'base' or 'total'
#basin_char is vector with 5 parameters for basin characteristics
#gw_hyd is a vector with 5 parameters for groundwater hydraulics
#flow is a vector with 6 parameters for streamflow

#REMOVE NEGATIVE FLOW VALUES
qin[qin<0]=NA

#REQUIRED UTILITY FUNCTIONS
#base_table; sur_store; sur_z; sur_q; dir_q; infiltration; recharge; bf_ci

#BASIN CHARACTERISTICS
area=basin_char[1] #DRAINAGE AREA
lb=basin_char[2] #BASIN LENGTH
x1=basin_char[3] #SCALING PARAMETER FOR BASE LENGTH
wb=basin_char[4] #BASE WIDTH FOR SATURATED FLOW
por=basin_char[5] #DRAINABLE POROSITY
ws=wb/2

#write('BASIN CHARACTERISTICS: AREA Lb X1 Wb POR','')
#write(basin_char,'')

#GROUNDWATER HYDRAULIC PARAMETERS:
alpha=gw_hyd[1] #SURFACE HYDRAULIC GRADIENT
beta=gw_hyd[2] #EXPONENT FOR Z = f(X/x1)
ks=gw_hyd[3] #HYDRAULIC CONDUCTIVITY FOR SURFACE
kb=gw_hyd[4] #HORIZONTAL HYDRAULIC CONDUCTIVITY FOR BASE DISCHARGE
kz=gw_hyd[5] #VERTICAL HYDRAULIC CONDUCTIVTIY FOR RECHARGE

#write('GROUNDWATER HYDRAULICS: ALPHA BETA Ks Kb Kz','')
#write(gw_hyd,'')

#FLOW METRICS
qthresh=flow[1] #BASE FLOW USED TO CALIBRATE kb AND INITIALIZE BASE STORAGE WHEN DATA ARE MISSING
rs=flow[2] #EXPONENTIAL RATE CONSTANT [1/T] FOR STORM FLOW
rb1=flow[3] #EXPONENTIAL RATE CONSTANT [1/T] FOR HIGH BASE FLOW
rb2=flow[4] #EXPONENTIAL RATE CONSTANT [1/T] FOR LOW BASE FLOW
prec=flow[5] #PRECISION OF STREAMFLOW VALUES
fr4rise=flow[6] #INCREASE IN STREAMFLOW AS FRACTION FOR DESIGNATING TIME STEPS WITH IMPULSES
qmean=mean(qin,na.rm=T)

#write('FLOW METRICS: Qthresh Rs Rb1 Rb2 Prec Frac4Rise','')
#write(flow,'')

#LOOKUP TABLE FOR BASE STORAGE-DISCHARGE RELATION WITH S,Q,X,Z,dz/dx
SBT=base_table(lb,x1,wb,beta,kb,qin,por)

#ERROR TOLERANCE USED SEQUENTIALLY TO REFINE IMPULSE
ifact=c(2,1.1)

#NUMBER OF TIME STEPS
p=length(qin)

#LOGICAL FOR BAD PARAMETER
bad_params=FALSE
##################################
#DESIGNATE RECESSIONAL TIME STEPS
##################################
#TIME SERIES OF FRACTIONAL DAILY CHANGE IN STREAMFLOW (NEGATIVE VALUE FOR RECESSION)
if(timestep=='day'){dq=c(0,qin[2:p]-qin[1:(p-1)])}
if(timestep=='hour'){dq=rep(0,p)
for(y in 25:p){dq[y]=qin[y]-max(qin[(y-24):(y-1)],na.rm=T)}} #USE MAXIMUM OF PREVIOUS 24 HOURS TO IDENTIFY RECESSIONAL TIME STEPS FOR HOURLY DATA
dqfr=dq/qin
dqfr[(dq==0)&(qin==0)]=0
dqfr[(dq<0)&(qin==0)]=1

#RISE IS LOGICAL VECTOR INDICATING WHEN CHANGE FROM PREVIOUS TIME STEP MEETS CRITERIA:
#GREATER THAN SPECIFICED FRACTIONAL CHANGE AND PRECISION
rise=(dqfr>fr4rise & dq>prec)
rise[is.na(rise)]=FALSE

#RECESS IS LOGICAL VECTOR INDICATING RECESSIONAL AND RELATIVELY CONSTANT DISCHARGE PERIODS (INCLUDING PERIODS WHEN STREAMFLOW RISES LESS THAN CRITERIA FOR "RISE"
recess=c(dqfr<=fr4rise | dq<prec)
recess[is.na(recess)]=FALSE

#CALCULATE CONSECUTIVE DAYS OF RECESSION, USE FOR WEIGHTING ERROR
recess_day=cumsum(recess)-cummax((!recess)*cumsum(recess))
###################################################################
#OUTPUT VARIABLES
X=rep(NA,p) #LONGITUDINAL LOCATION OF BASE WATER LEVEL INTERSECTION WITH SURFACE
qcomp=array(NA,dim=c(p,3)) #THREE FLOW COMPONENTS: surface flow; base flow; direct runoff from saturated areas
ETA=rep(NA, p) #STATE DISTURBANCES (POSITIVE VALUES REPRESENT INPUTS) [L3]
I=rep(NA,p) #PRECIPITATION CALCULATED FROM eta
Z=array(NA, dim=c(p,2)) #WATER SURFACE ELEVATION OF SURFACE (CHANNEL IS DATUM) AND BASE (BASIN OUTLET IS DATUM)
ST=array(NA, dim=c(p,2)) #STORAGE
EXC=array(NA,dim=c(p,2)) #EXCHANGES - INFILTRATION AND RECHARGE
#####################################################################
#CHECK PARAMETERS, END PROCESS IF PARAMETERS ARE BAD
if(any(c(lb,x1,wb,alpha,beta,ks,kb,kz,por,qthresh,-rs,-rb1,-rb2,prec,fr4rise)<0)) {bad_params=TRUE
    	write('Negative parameter(s)','')}
if(any(c(lb,x1,wb,alpha,beta,ks,kb,kz,por,qthresh,-rs,-rb1,-rb2,prec,fr4rise)==Inf)) {bad_params=TRUE
    	write('Infinite parameter(s)','')}
if(lb*wb > area) {bad_params=TRUE
    	write('lb x wb > area','')}
if(any(is.na(SBT))) {bad_params=TRUE
    	write('Cannot calculate discharge for base parameters')}

ts=1 #INITIAL TIME STEP
stts=ts #STARTING TIME STEP, stts, FOR ERROR CALCULATION

#ADVANCE TIME STEP TO FIRST STREAMFLOW OBSERVATION
while(is.na(qin[ts])){ts=ts+1; stts=ts}

#INITIALIZE VARIABLES
ts_ini=TRUE
qb_in=min(qin[ts],qmean) #INITIAL BASEFLOW IS SMALLER OF OBSERVED AND MEAN STREAMFLOW
qb_en=NA
xb_in=SBT$Xb[sum(SBT$Q<=qb_in)]
zb_in=SBT$Z[sum(SBT$Q<=qb_in)]
sb_in=SBT$S[sum(SBT$Q<=qb_in)]
sba=max(SBT$S)-sb_in #AVAILABLE BASE STORAGE

qs_in=max(0,qin[ts]-qb_in,na.rm=T)
zs_in=min(qs_in/(2*lb*ks*alpha),ws*alpha, na.rm=T)
ss_in=sur_store(lb,alpha,ws,por,zs_in)
ssa=sur_store(lb,alpha,ws,por,ws*alpha)-ss_in #AVAILABLE SURFACE STORAGE

infil_in=0
rech_in=recharge(lb,xb_in,ws,kz,zs_in,por)

#LOOP THROUGH TIME STEPS
while(ts<=p) {

#INITIALIZE TIME STEP USING STATE VARIABLES FOR PREVIOUS TIME STEP IF AVAILABLE
  if(!ts_ini){
    xb_in=xb_in=X[ts-1]
    zb_in=Z[ts-1,2]
    sb_in=ST[ts-1,2]
    qb_in=SBT$Q[sum(SBT$Xb<=xb_in)]

    zs_in=Z[ts-1,1]
    ss_in=ST[ts-1,1]
    qs_in=sur_q(lb,alpha,ks,zs_in)

#STORAGE CAPACITY AVAILABLE
    ssa=sur_store(lb,alpha,ws,por,ws*alpha)-ss_in
    sba=max(SBT$S)-sb_in #BASE ZONE
    rech_in=min(recharge(lb,xb_in,ws,kz,zs_in,por),sba+qb_in) #INITIAL RECHARGE LIMITED TO AVAILABLE BASE STORAGE CAPACITY + BASE FLOW
    qd=0
    infil_in=0

#IMPULSE (PPT) NEEDED TO GENERATE OBSERVED STREAMFLOW
#CONDITION FOR ALLOWING IMPULSE DURING A TIME STEP
    I[ts]=0 #SET IMPULSE TO ZERO
    etaest=max(c(0,qin[ts]-qb_in-qs_in),na.rm=T)

#INITIAL ESTIMATE OF IMPULSE REQUIRED FOR ADDITIONAL SURFACE FLOW
    if((ts>1)&(etaest>0)) {if(rise[ts]|rise[ts-1]) {I[ts]=etaest/(2*lb*ws)
      zs=zs_in
      qs=qs_in

#LOOP TO CALCULATE ADDITIONAL IMPULSE NEEDED TO REDUCE ETA
#USE PROGRESSIVELY SMALLER INCREMENTAL CHANGES IN IMPULSE (ifact) FOR ITERATIONS
#STOP WHEN ETA IS LESS THAN PRECISION OR 1% OF Q
      for(x in ifact) {i=I[ts];eta=etaest
        while((eta>max(prec,qin[ts]/100)) & (i>0)) {I[ts]=i;etaest=eta
          i=x*i
          infil=min(infiltration(lb,ws,ks,alpha,(zs_in+zs)/2,i),ssa) #LIMIT INFILTRATION TO AVAILABLE STORAGE
          ss=max(ss_in+infil-rech_in-qs,0) #UPDATE SURFACE STORAGE
          zs=sur_z(lb,alpha,ws,por,ss)
          qs=sur_q(lb,alpha,ks,zs) #SURFACE DISCHARGE
          qd=dir_q(lb,alpha,zs_in,i)+dir_q(lb,alpha,(zs-zs_in),i/2)+max(2*lb*(ws-zs_in/alpha)*(I[ts]-ks),0)
          eta=qin[ts]-qs-qd-qb_in}}}} #CLOSE RISE, X-FACTOR, AND ETA CONDITIONS

  infil_in=min(infiltration(lb,ws,ks,alpha,zs_in,I[ts]),ssa)} #CLOSE INITIAL CALCULATIONS WHEN STREAMFLOW RECORD IS AVAILABLE (NOT PROJECTION)

#END OF TIME STEP CALCULATIONS
  ss_en=max(ss_in+infil_in-rech_in-qs_in,0)
  zs_en=sur_z(lb,alpha,ws,por,ss_en)
  qs_en=sur_q(lb,alpha,ks,zs_en)
  infil_en=min(infiltration(lb,ws,ks,alpha,zs_en,I[ts]),ssa,na.rm=T)
  rech_en=min(recharge(lb,xb_in,ws,kz,zs_en,por),sba+qb_in)
  sb_en=max(sb_in+rech_en-qb_in,0)
  xb_en=SBT$Xb[max(sum(SBT$S<sb_en),1)]
  zb_en=SBT$Z[max(sum(SBT$S<sb_en),1)]
  qb_en=SBT$Q[max(sum(SBT$S<sb_en),1)]

############################################################
#FINAL CALCS FOR TIME STEP
  qcomp[ts,1]=(qs_in+qs_en)/2 #SURFACE FLOW
  qcomp[ts,2]=(qb_in+qb_en)/2 #BASE FLOW

  EXC[ts,1]=(infil_in+infil_en)/2
  EXC[ts,2]=(rech_in+rech_en)/2

#FOR INITIAL TIME STEP
  if(ts_ini){ST[ts,]=c(ss_en,sb_en);Z[ts,]=c(zs_en,zb_en)}

#FOR TIME STEPS WHEN STATES ARE AVAILABLE FOR PREVIOUS TIME STEP
  if(!ts_ini) {
    ST[ts,1]=max(ST[ts-1,1]+EXC[ts,1]-qcomp[ts,1]-EXC[ts,2],0)
    ST[ts,1]=min(ST[ts,1],sur_store(lb,alpha,ws,por,ws*alpha))
    Z[ts,1]=sur_z(lb,alpha,ws,por,ST[ts,1])

    ST[ts,2]=max(ST[ts-1,2]+EXC[ts,2]-qcomp[ts,2],0)
    ST[ts,2]=min(ST[ts,2],max(SBT$S))
    Z[ts,2]=SBT$Z[max(sum(SBT$S<=ST[ts,2]),1)]

#DIRECT RUNOFF INCLUDEs ADDITIONAL SATURATED AREA X HALF OF RAINFALL(EXCESS AFTER INFILTRATION) AND ANY PRECIPITATION THAT EXCEEDS INFILTRATION RATE
    qcomp[ts,3]=dir_q(lb,alpha,zs_in,I[ts]) + dir_q(lb,alpha,(Z[ts,1]-zs_in),I[ts]/2) + max(2*lb*(ws-zs_in/alpha)*(I[ts]-ks),0)} #DIRECT RUNOFF

  ETA[ts]=qin[ts]-sum(qcomp[ts,1:3]) #STREAMFLOW RESIDUAL
  X[ts]=SBT$Xb[max(sum(SBT$S<=ST[ts,2]),1)]
  ts=ts+1
  ts_ini=FALSE

  } #CLOSE CONDITION ts<p
########
#OUTPUT
########
if(error_basis=='base') {q4er=qcomp[,2]}
if(error_basis=='total') {q4er=rowSums(qcomp)}

#ADJUSTED PERCENT ERROR
APE=(q4er-qin)/(qin+prec)
APE[qin==0 & q4er==0]=0
Weight=1-exp(rs*recess_day) #WEIGHT VARIES FROM O TO 1 WITH INCREASING LENGTH OF RECESSION
Weight[APE>0]=1 #WEIGHT OF 1 IS ASSIGNED TO OVER PREDICTION

tmp=data.frame(dys,qin,rowSums(qcomp),qcomp,ETA,ST[,1:2],I,Z,EXC,recess_day,APE,Weight,stringsAsFactors=FALSE)
tmp[,2:16]=signif(tmp[,2:16],6)

dimnames(tmp)[[2]]=c('Date','Qob.L3','Qsim.L3','SurfaceFlow.L3','Baseflow.L3','DirectRunoff.L3','Eta.L3','StSur.L3','StBase.L3','Impulse.L','Zs.L','Zb.L','Infil.L3','Rech.L3','RecessCount.T','AdjPctEr','Weight')

tmp.ci=bf_ci(tmp)
tmp=cbind(tmp,signif(tmp.ci[[2]],3))

assign('bfs_out',tmp, envir=.GlobalEnv)

#RETURN MEAN ABSOLUTE WEIGHTED APE
signif(sum(abs(APE*Weight),na.rm=T)/sum(Weight,na.rm=T),3)}

