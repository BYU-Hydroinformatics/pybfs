#########################################################################
#UTILITY FUNCTIONS NEEDED FOR BASE FLOW SEPARATION MODEL:
#base_table; sur_store; sur_z; sur_q; dir_q; infiltration; recharge; bf_ci; flow_metrics
#########################################################################
#DESCRIPTION OF VARIABLES USED AS ARGUMENTS:
#a - HYDRAUlIC GRADIENT FOR SURFACE (CONSTANT)
#b - EXPONENT FOR SATURATED HEIGHT AS A FUNCTION OF LENGTH:Z=[(xb-X0)/x1]^BETA
#lb - BASIN LENGTH
#x1 - LONGITUDINAL SCALING PARAMETER
#xb - LONGITUDINAL LOCATION OF BASE WATER LEVEL INTERSECTION WITH SURFACE
#ws - WIDTH OF SURACE ZONE (ONE SIDE OF CHANNEL)
#wb - WIDTH OF BASE ZONE
#kb - HORIZONAL HYDRAULIC CONDUCTIVITY OF BASE
#ks - HYDRAULIC CONDUCTIVITY OF SURFACE
#kz - VERTICAL HYDRAULIC CONDUCTIVITY
#por - DRAINABLE POROSITY
#zb - BASE WATER SURFACE ELEVATION AT INTERSECTION WITH SURFACE
#zs - SATURATED THICKNESS OF SURFACE
############################################################################################
#LOOKUP TABLE FOR BASE VARIABLES: SATURATED HEIGHT, INTERSECTION, HYDRAULIC GRADIENT, STORAGE, AND DISCHARGE
base_table=function(lb,x1,wb,b,kb,qin,por) {

tmp.q=qin[!is.na(qin)]
tmp.range=log(c(min(tmp.q[tmp.q>0])/10,max(tmp.q[tmp.q>0])),10)
qq=10^(seq(tmp.range[1],tmp.range[2],(tmp.range[2]-tmp.range[1])/999))
qq=c(0,qq)

#LONG PROFILE OF BASE WATER SURFACE:  Z IS ELEVATION, X IS LONGITUDINAL DISTANCE FROM DOWNSTREAM END 
#z = (x/x1)^b
#x = x1 * z^1/b

#DERIVATION OF HYDRAULIC GRADIENT AS A FUNCTION OF Z
#dzdx = b / (x1^b) * x^(b-1)
#dzdx = b / (x1^b) * (x1 * z^(1/b))^(b-1)
#dzdx = b * (x1^(b-1) / x1^b) *  z^((b-1)/b)
#dzdx = b / x1 *  z^((b-1)/b)

#DIFFERENTIAL SOLUTION WHERE I[] IS INTEGRAL 
#q = wb * kb * I[dzdx dz]
#q = wb * kb * I[b/x1 * z*((b-1)/b dz]
#q = wb * kb * b/x1 * I[z^((b-1)/b) dz]

#I[z^((b-1)/b) dz] = b / (2*b - 1) * z^((2*b-1)/b)
#IF b == 0.5, I[z^((b-1)/b) dz]= ln(z)-1

#q = wb * kb * b^2 / (x1 * (2 * b - 1)) * z^((2*b-1)/b)
#q * (x1 * (2 * b - 1))  = wb * kb * b^2 * z^((2*b-1)/b)
#z = (q * x1 * (2 * b - 1) / (wb * kb * b^2))^(b/(2*b-1))

#q = wb * kb * b^2 / (x1 * (2 * b - 1)) * (x/x1)^(2*b-1)
#q * x1 * (2 * b - 1) / (wb * kb * b^2) = (x/x1)^(2*b-1)
#x = x1 * (q * x1 * (2 * b - 1) / (wb * kb * b^2))^(1/(2*b-1))
#x = x1^(2*b/(2*b-1)) * (q * (2 * b - 1) / (wb * kb * b^2))^(1/(2*b-1))

#CALCULATE Z, X, AND STORAGE FOR QUANTILES OF Q RANGE
if(b!=0.5){z = (qq * x1 * (2 * b - 1) / (wb * kb * b^2))^(b/(2*b-1))}
if(b==0.5){z = exp(qq * 2* x1/(wb * kb))}

x = x1 * z^(1/b)

#s = wb * por * [I(z dx) +(lb - x) * z]
s = wb * por * (1/x1^b * 1/(b+1) * x^(b+1) + (lb - x) * z)

BT=data.frame(signif(x,5),signif(z,5),signif(s,5),signif(qq,5))
dimnames(BT)[[2]]=c('Xb','Z','S','Q')

BT[BT$Xb<lb,]}
#################################################################
#SURFACE STORAGE IS CALCULATED AS THE DIFFERENCE BETWEEN 
#THE TOTAL SURFACE ZONE VOLUME AND 
#THE UNSATURATED PORTION BELOW Z (SURFACE SATURATION THICKNESS)
#BOTH ARE TRIANGULAR PRISMS
sur_store=function(lb,a,ws,por,zs) {
z=min(c(ws*a,zs,0),na.rm=T)
lb*(2*ws*zs-zs^2/a)*por}
##################################################################
#SATURATED THICKNESS IN SURFACE ZONE 
sur_z=function(lb,a,ws,por,ss){
a1=1/(2*a); b1=-2*ws; c1=ss/(lb*por)
if((b1^2-4*a1*c1)<0) {ws*a} else {(-b1-sqrt(b1^2-4*a1*c1))/(2*a1)}}
###################################################################
#DISCHARGE FOR SURFACE ZONE
sur_q=function(lb,a,ks,z) {2*lb*z*a*ks} #Z is WATER SURFACE ELEVATION AT SURFACE OF SURFACE ZONE (DATUM IS CHANNEL)
#########################################################################
#DIRECT RUNOFF
dir_q=function(lb,a,z,i){2*lb*z/a*i} #z/b IS THE SATURATED SURFACE WIDTH
##########################################################################
#INFILTRATION BASED ON UNSATURATED SURFACE AREA
infiltration=function(lb,ws,ks,a,zs,i) {(2*lb*(ws-zs/a)*min(i,ks))} #Z IS WATER SURFACE ELEVATION AT SURFACE OF SURFACE ZONE (DATUM IS CHANNEL), i IS IMPULSE DEPTH
#####################################################################
#RECHARGE BASED ON SATURATED SURFACE DEPTH AND UNSATURATED BASE AREA
recharge=function(lb,xb,ws,kz,zs,por) {(lb-xb)*2*ws*min(zs*por,kz)}
#####################################################################
#CREDIBLE PREDICTION INTERVALS FOR STREAMFLOW
bf_ci=function(bfs_out) {
tmp.error=(bfs_out$Qsim.L3-bfs_out$Qob.L3)/bfs_out$Qsim.L3
tmp.error[abs(tmp.error==-Inf)]=1
tmp.error[bfs_out$DirectRunoff.L3>0]=NA
tmp.q=bfs_out$Qsim.L3
tmp.q[bfs_out$DirectRunoff.L3>0]=NA
qnts=quantile(tmp.q[tmp.q>0],p=seq(0.05,0.95,0.05),na.rm=TRUE)

#STAND ALONE TABLE
ci_table=data.frame(array(dim=c(0,5)))
for(x in 2:18) {y=(tmp.q>qnts[x-1]) & (tmp.q<qnts[x+1])
ci_table=rbind(ci_table,c(qnts[x],quantile(tmp.error[y],p=c(0.05,0.5,0.95),na.rm=TRUE)))}
dimnames(ci_table)[[2]]=c('Qsim.L3.T','FrEr0.05','FrEr0.50','FrEr0.95')

#CI FOR EACH DAILY VALUE
qnt=rep(NA,length(tmp.q))
for(t in 1:length(tmp.q)){if(bfs_out$DirectRunoff[t] %in% 0) {qnt[t]=match(TRUE,bfs_out$Qsim.L3[t]<(qnts[1:18]+qnts[2:19])/2)}}
ci=data.frame(bfs_out$Qsim.L3*(1-ci_table$FrEr0.95[qnt]),bfs_out$Qsim.L3*(1-ci_table$FrEr0.05[qnt]))
dimnames(ci)[[2]]=c('CB0.05','CB0.95')

list(ci_table,ci)}
######################################################################################################
#FUNCTION TO CALCULATE STREAMFLOW METRICS USED FOR BASE-FLOW SEPARATION FROM A STREAMFLOW TIME SERIES
flow_metrics=function(qin,timestep,fr4rise) {

#ARGUMENTS:
#qin IS A NUMERIC VECTOR WITH STREAMFLOW TIME SERIES;
#timestep IS EITHER "day" or "hour";

#OUTPUTS:
#Qthresh IS A LOWER THRESHOLD ON BASEFLOW SPECIFIED BY THE LARGER OF MINIMUM STREAMFLOW (NOT INCLUDING NO FLOW) AND FLOW AT WHICH RECESSION RATES INCREASE
#Rs IS EXPONENTIAL MODEL COEFFICIENT FOR RAPID RECESSION (STORM FLOW): 95TH PERCENTILE OF STANDARDIZED (PERCENT CHANGE) 2-DAY RECESSIONS WHERE Q0 IS RISE IN STREAMFLOW FOR DAY 0 [1/TIMESTEP], USING THE RISE IN STREAMFLOW ACCOMMODATES RAPID RECESSION AFTER SMALL EVENTS DURING PERIODS OF LOW BASE FLOW
#Rb1 IS EXPONENTIAL MODEL COEFFICIENT FOR SLOW RECESSION (BASE FLOW): 95TH PERCENTILE OF STANDARDIZED (PERCENT CHANGE) 10-DAY RECESSIONS WHERE Q0 IS STREAMFLOW FOR DAY 0 [1/TIMESTEP]
#Rb2 IS EXPONENTIAL MODEL COEFFICIENT FOR INTERMEDIATE RECESSION (HIGH BASE FLOW) : 50TH PERCENTILE OF STANDARDIZED (PERCENT CHANGE) 10-DAY RECESSIONS WHERE Q0 IS STREAMFLOW FOR DAY 0 [1/TIMESTEP]
#prec IS THE PRECISION BASED ON THE GREATER OF MEASUREMENT RESOLUTION, MEDIAN DAILY RANGE OF HOURLY STREAMFLOW (DOES NOT APPLY TO DAILY STREAMFLOW)  OR 10% MEASUREMENT ERROR ON THE 1% FLOW: errors less than prec are ignored, baseflow less than prec is zero, mean percent error for zero flow is calculated using prec as the basis

library(quantreg)

xx=length(qin)

#SET NEGATIVE FLOWS TO NA
qin[qin<0]=NA

#PRECISION OF MEASURED BASEFLOW: 
Q01=quantile(qin,0.01,na.rm=T)
tmp=qin-Q01  
prec=min(tmp[tmp>0], na.rm=T) #SMALLEST POSITIVE DIFFERENCE BETWEEN STREAMFLOW AND Q01 IS USED TO INDICATE LOW-FLOW MEASUREMENT RESOLUTION 
prec=max(prec,0.1*Q01,na.rm=T) #INCREASE PRECISION IF 10% OF Q01 (NOMINAL LOW-FLOW MEASUREMENT ERROR) IS LARGER

#INCREASE PRECISION TO MEDIAN DAILY RANGE (PRESUMED TO BE NOISE) FOR HOURLY DATA
if(timestep=='hour') {tmp.x=cbind(c(rep(1,24),c(2:(xx-23))),c(1:xx))
tmp.range=c()

for (t in 1:xx) {tmp.range=c(tmp.range,max(qin[tmp.x[t,1]:tmp.x[t,2]],na.rm=T)-min(qin[tmp.x[t,1]:tmp.x[t,2]],na.rm=T))} #DAILY RANGE
prec=max(prec,median(tmp.range,na.rm=T),na.rm=T)}

if(timestep=='day') {x=2}
if(timestep=='hour') {x=48}

#INTERPOLATE STREAMFLOW  IF MISSING PERIOD IS LESS THAN 2 DAYS
tmp1=c(1:xx)
tmp2=tmp1-cummax(tmp1*(!is.na(qin)))
tmp3=tmp1-cummax(tmp1*(!is.na(qin[xx:1])))
tmp3=tmp3[xx:1]
tmp=apply(cbind(tmp2,tmp3),1,max)

qint=cbind(tmp,qin,qin)
for(y in 2:(xx-1)) {if(is.na(qint[y,2])){qint[y,2]=qint[y-1,2]}
if(is.na(qint[xx-y+1,3])){qint[xx-y+1,3]=qint[xx-y+2,3]}}
qin[is.na(qin)&(qint[,1]<x)]=(qint[is.na(qin)&(qint[,1]<x),2]+qint[is.na(qin)&(qint[,1]<x),3])/2

#IDENTIFY RECESSION DAYS
#USING WINDOW OF PREVIOUS 2 DAYS
rp_indices=array(c(1:xx),dim=c(xx,1)) #DAILY INDICES FOR THE TIME STEPS USED TO CALCULATE RECESSION: INDICES FOR TS=Y ARE Y,Y-1,...,Y-X
for(y in 1:x){rp_indices=cbind(rp_indices,c((1-y):(xx-y)))}

rec=rep(FALSE,xx) #LOGICAL VECTOR INDICATING RECESSION TIME STEPS
rise=rep(FALSE,xx) #LOGICAL VECTOR INDICATING RISE TIME STEPS
qrise=rep(0,xx) #RISE OVER WINDOW USED AS BASIS TO CALCULATE RECESSION RATE
qmax=rep(FALSE,xx) #TIME STEP HAS MAX Q FOR WINDOW

#CALCULATE THE CHANGE IN STREAMFLOW OVER PREVIOUS 2 DAYS FOR EACH TIME-STEP
for(y in (x+1):xx){qp=qin[rp_indices[y,]]
if(all(!is.na(qp))) {dq=qp[1:x]-qp[2:(x+1)] #CHANGE IN STREAMFLOW FOR EACH TIME STEP IN WINDOW OF X TIME STEPS
q4rise=max(fr4rise*qp[1],prec,na.rm=T) #THRESHOLD STREAMFLOW FOR RISE
qrise[y]=sum(dq,na.rm=T) #INCREASE FOR WINDOW
rise[y]=(qrise[y]>q4rise)
if(qrise[y]<q4rise){qrise[y]=NA}

rec[y]=(all(dq<prec) & (sum(dq,na.rm=T)<0)) #RECESSION TIME STEPS: CHANGES FOR ALL TIME STEPS IN WINDOW MUST BE LESS THAN PREC AND CUMULATIVE CHANGE MUST BE LESS THAN ZERO; CREATES A BUFFER FOR X TIME STEPS AFTER RISE THAT WILL NOT BE CONSIDERED RECESSION TIME STEPS
qmax[y]=(qp[1]==max(qp,na.rm=T))}} #TIME STEP HAS MAXIMUM STREAMFLOW FOR WINDOW

#2-DAY RECESSION RATES FOLLOWING A RISE
dq=c(rep(NA,x),(qin[(x+1):xx]-qin[1:(xx-x)])/qrise[1:(xx-x)]) #CHANGE IN STREAMFLOW AS FRACTION OF RISE
dq[!rec]=NA #LIMIT dq TO RECESSION PERIODS
dq[!c(rep(FALSE,x),qmax[1:(xx-x)])]=NA #SET dq TO NA FOR TIME STEP WITH MAX STREAMFLOW AT START OF WINDOW
dq=quantile(dq,0.95,na.rm=T)
if(!is.na(dq) & (dq<(-1))){dq=(-0.9)} #REDUCE 100% RECESSION TO 90% TO PREVENT ERROR FOR RECESSION COEFFICIENT (NEXT LINE)
Rs=log(1+dq)/x #RECESSION COEFFICIENT FOR EXP MODEL

#10-DAY RECESSIONS
if(timestep=='day') {x=10}
if(timestep=='hour') {x=240}
rec=rep(FALSE,xx) #LOGICAL VECTOR INDICATING RECESSIONAL TIME STEPS

#CALCULATE THE CHANGE IN STREAMFLOW OVER PREVIOUS 10 DAYS FOR EACH TIME-STEP
rp_indices=array(c(1:xx),dim=c(xx,1)) #DAILY INDICES FOR THE TIME STEPS USED TO CALCULATE RECESSION: INDICES FOR TS=Y ARE Y,Y-1,...,Y-X
for(y in 1:x){rp_indices=cbind(rp_indices,c((1-y):(xx-y)))}

for(y in (x+1):xx){qp=qin[rp_indices[y,]] #STREAMFLOW FOR PREVIOUS X TIME STEPS (WINDOW)
if(any(!is.na(qp))){dq=qp[1:x]-qp[2:(x+1)] #CHANGE IN STREAMFLOW FOR EACH TIME STEP IN WINDOW
q4rise=fr4rise*max(c(qp[1],prec),na.rm=T) #THRESHOLD STREAMFLOW FOR RISE
rec[y]=(sum(dq,na.rm=T)<=0) & (sum(dq[dq>0])<q4rise) #RECESSION PERIOD: CUMULATIVE CHANGE MUST BE LESS THAN ZERO AND NO TIME STEP IN WINDOW CAN HAVE A CHANGE GREATER THAN RISE (CREATES A BUFFER FOR X TIME STEPS AFTER RISE)
rec[y]=(qp[1]==min(qp,na.rm=T))} #TIME STEP MUST HAVE MINIMUM STREAMFLOW FOR WINDOW 
} #CLOSE Y LOOP 

#10-DAY RECESSION RATES
dq=c(rep(NA,x),(qin[(x+1):xx]-qin[1:(xx-x)])/qin[1:(xx-x)]) #CHANGE IN STREAMFLOW FOR TIME STEP FROM PREVIOUS X TIME STEP AS FRACTION OF STREAMFLOW FOR PREVIOUS X TIME STEP
dq[!rec]=NA #LIMIT dq TO RECESSION PERIODS

#LIMIT RECESSION CALCULATIONS TO DAYS WITH POSITIVE STREAMFLOW WHEN STREAM DID NOT DRY COMPLETELY (dq > -1)
tmp.d=qin>0
tmp.d[is.na(qin)]=FALSE
tmp.d[dq==-1]=FALSE
logq=log(qin[tmp.d],10)
rb=log(1+dq[tmp.d],10)/x

#SET Qthresh TO MINIMUM STREAMFLOW EXCLUDING NO FLOW PERIODS
Qthresh=min(qin[qin>0],na.rm=T)
qmean=mean(qin,na.rm=T)

#QUANTILES OF STREAMFLOW, EXCLUDING NO-FLOW PERIODS, UP TO MEAN
tmp.qq=quantile(qin[(qin>0) & (qin<=mean(qin,na.rm=T))],p=seq(0,1,0.01),na.rm=T) 

#CALCULATE MEDIAN RECESSION RATES FOR STREAMFLOW QUANTILES
tmp.r=c()
for(qq in tmp.qq) {tmp.r=c(tmp.r,quantile(dq[qin<qq],p=0.95,na.rm=T))} #RECESSION RATE FOR Q<QQ

#SET QTHRESH TO STREAMFLOW WITH MAXIMUM (LEAST NEGATIVE) RECESSION RATE, USED TO FILTER OUT LOW FLOWS WHEN RECESSION INCREASES BECAUSE OF STREAM DRYING 
tmp.qqthresh=match(max(tmp.r,na.rm=T),tmp.r)
tmp.qthresh=tmp.qq[tmp.qqthresh]
Qthresh=max(c(Qthresh,tmp.qthresh),na.rm=T)

#LOWER LIMIT ON FASTEST RECESSION RATES AS A FUNCTION OF LOG STREAMFLOW
logq=logq[!is.na(rb)]
rb=rb[!is.na(rb)]

if((length(rb)<9)| (min(logq)==max(logq))){rb10=c(NA,NA)} else {

  tmp=summary(rq(rb~logq,tau=0.1))
  rb10=tmp$coefficients[1:2] #INTERCEPT AND SLOPE
  names(rb10)=c('RbIntcpt','RbSlope')}

assign('rb10',rb10,envir=.GlobalEnv)

#RECESSION RATE FOR MEAN STREAMFLOW
Rb1=rb10[1]+rb10[2]*log(qmean,10)
Rb1=min(Rb1,-0.002) #LIMIT TO -0.002

#RECESSION RATE FOR LOW BASEFLOW (10TH PERCENTILE OF FILTERED STREAMFLOW)
tmp.logq10=quantile(logq,p=0.1,na.rm=T)
Rb2=rb10[1]+rb10[2]*tmp.logq10
Rb2=min(Rb2,-0.001) #LIMIT TO -0.001

out=c(Qthresh,Rs,Rb1,Rb2,prec,fr4rise)
names(out)=c('Qthresh','Rs','Rb1','Rb2','Prec','Fr4Rise')
out}
