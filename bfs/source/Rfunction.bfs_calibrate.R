#CALIBRATE A SITE
bfs_calibrate=function(tmp.site,tmp.area,tmp.q,dys) {

#RUN FLOW METRICS FUNCITON TO GENERATED VARIABLES IN flow VECTOR
  flow=flow_metrics(tmp.q,timestep='day',fr4rise=0.05)

#ONLY PROCEED WITH CALIBRATION IF flow_metrics RETURNS VALUES FOR ALL PARAMETER
  if(any(c(is.na(flow),flow==-Inf))) {bf_params=NA; bff=NA; ci_table=NA; bfs_out=NA} else {
    Qthresh=flow[1] #BASE FLOW USED TO CALIBRATE Kb AND INITIALIZE BASE STORAGE WHEN DATA ARE MISSING
    Rs=flow[2] #STORMFLOW RECESSION COEFFICIENT BASED ON 2-DAY RECESSION, EXCLUDES TIME STEPS USED FOR CALIBRATION
    Rb1=flow[3] #BASEFLOW RECESSION COEFFICIENT AT Qmean, BASED ON 10-DAY RECESSION
    Rb2=flow[4] #BASEFLOW RECESSION RATE AT Qthresh
    Prec=flow[5] #PRECISION OF LOW FLOW MEASUREMENTS
    Frac4Rise=flow[6] #THRESHOLD FOR RISE IN STREAMFLOW AS FRACTION

    Qmean=mean(tmp.q[tmp.q>=0],na.rm=TRUE) #LIMIT MEAN STREAMFLOW TO ZERO OR POSITIVE FLOWS (FOR TIDAL SYSTEMS)

    RbI=rb10[1] #INTERCEPT OF 10th QUANTILE OF Rb ~ LOG Q
    RbS=rb10[2] #SLOPE

#INITIALIZE PARAMETERS WITH NOMINAL ESTIMATES
    Lb=2*(tmp.area/2)^0.5 #BASIN LENGTH [L]
    Wb=tmp.area/Lb/10 #WIDTH OF BASE
    Ws=Wb/2
    POR=0.15 #DRAINABLE POROSITY

    ALPHA=0.01 #SURFACE HYDRAULIC GRADIENT
    BETA=1 #BASE SURFACE EXPONENT
    X1=1/ALPHA #BASE GRADIENT IS EQUAL TO SURFACE GRADIENT

    tmp=ini_params(tmp.area,Lb,X1,Wb,POR,BETA,Rb1,tmp.q)
    Lb=tmp$Lb
    Wb=tmp$Wb
    Kb=tmp$Kb

    Ks=(1-exp(Rs))*POR*(3/4*Ws)/ALPHA #SURFACE HYDRAULIC CONDUCTIVITY
    Kz=10*Qmean/(Lb*Wb) #VERTICAL HYDRAULIC CONDUCTIVITY OF BASE, MAXIMUM RECHARGE RATE

    basin_char=c(tmp.area,Lb,X1,Wb,POR)
    gw_hyd=c(ALPHA,BETA,Ks,Kb,Kz)

    Error=bfs(tmp.q,dys,timestep='day',error_basis='base',basin_char,gw_hyd,flow)
    BFF=sum(bfs_out$Baseflow.L3,na.rm=T)/sum(bfs_out$Qob.L3,na.rm=T)
    bf_params=data.frame(tmp.site,tmp.area,Lb,X1,Wb,POR,ALPHA,BETA,Ks,Kb,Kz,Qthresh,Rs,Rb1,Rb2,Prec,Frac4Rise,Error,BFF)
#########################################################################################################
#STEP 1. CALIBRATE ASSUMING BETA = 1
#########################################################################################################
    X=c(Lb,Wb,ALPHA,Ks,Kb,Kz)

    LOGX=log(X,10)

    tmp=optim(LOGX,cal_initial,tmp.q=tmp.q,dys=dys,timestep='day',error_basis='base',basin_char=basin_char,gw_hyd=gw_hyd,flow=flow,qmean=Qmean,
      control=list(maxit=1000,parscale=LOGX,reltol=0.01))

#UPDATE PARAMETERS AND PARAMETER SETS
    if(all(is.finite(tmp$par))){
      Lb=10^tmp$par[1];Wb=10^tmp$par[2];ALPHA=10^tmp$par[3];Ks=10^tmp$par[4];Kb=10^tmp$par[5];Kz=10^tmp$par[6]
      basin_char=c(tmp.area,Lb,X1,Wb,POR)
      gw_hyd=c(ALPHA,BETA,Ks,Kb,Kz)
      
      X=c(Lb,Wb,ALPHA,Ks)

      LOGX=log(X,10)

      tmp=optim(LOGX,cal_surface,tmp.q=tmp.q,dys=dys,timestep='day',error_basis='base',basin_char=basin_char,gw_hyd=gw_hyd,flow=flow,
        control=list(maxit=1000,parscale=LOGX,reltol=0.01))}

    if(all(is.finite(tmp$par))){
      Lb=10^tmp$par[1];Wb=10^tmp$par[2];ALPHA=10^tmp$par[3];Ks=10^tmp$par[4]
      basin_char=c(tmp.area,Lb,X1,Wb,POR)
      gw_hyd=c(ALPHA,BETA,Ks,Kb,Kz)

      Error=bfs(tmp.q,dys,timestep='day',error_basis='base',basin_char,gw_hyd,flow)
      tmp.bf=bfs_out$Baseflow.L3
      tmp.ov=bfs_out$Baseflow.L3>bfs_out$Qob.L3
      tmp.ov[is.na(tmp.ov)]=FALSE
      tmp.bf[tmp.ov]=bfs_out$Qob.L3[tmp.ov]
      BFF=sum(tmp.bf[!is.na(bfs_out$Qob.L3)])/sum(bfs_out$Qob.L3[!is.na(bfs_out$Qob.L3)])}

    bf_params=rbind(bf_params,data.frame(tmp.site,tmp.area,Lb,X1,Wb,POR,ALPHA,BETA,Ks,Kb,Kz,Qthresh,Rs,Rb1,Rb2,Prec,Frac4Rise,Error,BFF))
##################################################################################################################################
#STEP 2. CALIBRATE NON-LINEAR BASEFLOW FUNCTION FOR RECESSION RATES AT BETA = 1 TO 20, PICK BEST BETA AND RE-CALIBRATE BF_SEP
##################################################################################################################################
    tmp.out=array(dim=c(0,6))

    n=0
    b=0.5
    continue=TRUE
    
    while(continue) {b=b+0.1

#DERIVATION OF X1 GIVEN x and q
#Q = Wb * Kb I[dz/dx dz]
#Q = Wb * Kb I[b/X1 z^((b-1)/b) dz]
#Q = Wb * Kb * b / X1 * (b/(2* b-1)) * z^((2* b-1)/b)
#Q = Wb * Kb * b^2 / (X1 * (2* b - 1)) * ((x/X1)^b)^((2 * b-1)/b)
#Q = Wb * Kb * b^2 / (X1 * (2* b - 1)) * (x/X1)^(2 * b-1)
#X1^(2*b) = Wb * Kb * b^2 / ((2*b - 1) * Q) * x^(2 * b-1)
#X1 = (Wb * Kb * b^2 / ((2*b - 1) * Q) * x^(2 * b-1))^(1/(2*b))

#X1 SO THAT Qb=QMEAN @ Xb=Lb/2 SEE DERIVATION ABOVE
      X1=(Wb * Kb * b^2 / ((2*b - 1) * Qmean) * (Lb/2)^(2 * b-1))^(1/(2*b))
      
      if(is.finite(X1)) {X=c(Lb,X1,Wb,Kb)

        tmp=optim(X,cal_basetable,b=b,params=c(tmp.area,POR,Qmean,Qthresh,RbI,RbS),tmp.q=tmp.q,
          control=list(maxit=1000,parscale=X,reltol=0.01))

        if(all(is.finite(tmp$par))) {n=n+1
          tmp.out=rbind(tmp.out,c(tmp$par[1],tmp$par[2],tmp$par[3],b,tmp$par[4],tmp$value))}}

      if(b>10) {if(abs((tmp$value-tmp.out[n,6])/tmp$value)<0.001) {continue=FALSE}}
      if(b==20) {continue=FALSE}}

    X=tmp.out[match(min(tmp.out[,6]),tmp.out[,6]),1:5]

    Lb=X[1];X1=X[2];Wb=X[3]
    basin_char=c(tmp.area,X[1],X[2],Wb,POR)

    BETA=X[4];Kb=X[5]
    gw_hyd=c(ALPHA,BETA,Ks,Kb,Kz)

    X=c(Lb,Wb,Kb,Kz)

    tmp=optim(X,cal_base,tmp.q=tmp.q,dys=dys,timestep='day',error_basis='base',basin_char=basin_char,gw_hyd=gw_hyd,flow=flow,
      control=list(maxit=1000,parscale=X,reltol=0.01))
    tmp$par[tmp$par<0]=NA

#UPDATE PARAMETERS AND PARAMETER SETS
    if(all(is.finite(tmp$par))){
      Lb=tmp$par[1];Wb=tmp$par[2];Kb=tmp$par[3];Kz=tmp$par[4]
      
      X=c(Lb,Wb,ALPHA,Ks)
      LOGX=log(X,10)
      tmp=optim(LOGX,cal_surface,tmp.q=tmp.q,dys=dys,timestep='day',error_basis='base',basin_char=basin_char,gw_hyd=gw_hyd,flow=flow,
        control=list(maxit=1000,parscale=LOGX,reltol=0.01))

      if(all(is.finite(tmp$par))){
        Lb=10^tmp$par[1];Wb=10^tmp$par[2];ALPHA=10^tmp$par[3];Ks=10^tmp$par[4]
        basin_char=c(tmp.area,Lb,X1,Wb,POR)
        gw_hyd=c(ALPHA,BETA,Ks,Kb,Kz)}

      Error=bfs(tmp.q,dys,timestep='day',error_basis='base',basin_char,gw_hyd,flow)
      tmp.bf=bfs_out$Baseflow.L3
      tmp.ov=bfs_out$Baseflow.L3>bfs_out$Qob.L3
      tmp.ov[is.na(tmp.ov)]=FALSE
      tmp.bf[tmp.ov]=bfs_out$Qob.L3[tmp.ov]
      BFF=sum(tmp.bf[!is.na(bfs_out$Qob.L3)])/sum(bfs_out$Qob.L3[!is.na(bfs_out$Qob.L3)])

      bf_params=rbind(bf_params,data.frame(tmp.site,tmp.area,Lb,X1,Wb,POR,ALPHA,BETA,Ks,Kb,Kz,Qthresh,Rs,Rb1,Rb2,Prec,Frac4Rise,Error,BFF))}
#######################################################################################################
#STEP 3. SELECT BEST PARAMETERS TO MAXIMIZE BASEFLOW AND RE-CALIBRATE SURFACE PARAMETERS
#######################################################################################################
    g=match(max(bf_params$BFF[-1],na.rm=T),bf_params$BFF)
    Lb=bf_params$Lb[g]
    X1=bf_params$X1[g]

    Wb=bf_params$Wb[g]
    basin_char=c(tmp.area,Lb,X1,Wb,POR)

    BETA=bf_params$BETA[g]
    Kb=bf_params$Kb[g]
    Ks=bf_params$Ks[g]
    Kz=bf_params$Kz[g]
    gw_hyd=c(ALPHA,BETA,Ks,Kb,Kz)

    X=c(Lb,Wb,ALPHA,Ks,Kb,Kz)
    LOGX=log(X,10)
    tmp=optim(LOGX,cal_initial,tmp.q=tmp.q,dys=dys,timestep='day',error_basis='base',basin_char=basin_char,gw_hyd=gw_hyd,flow=flow, control=list(maxit=1000,parscale=LOGX,reltol=0.01))
    
    if(all(is.finite(tmp$par))) {
      Lb=10^tmp$par[1];Wb=10^tmp$par[2];ALPHA=10^tmp$par[3];Ks=10^tmp$par[4];Kb=10^tmp$par[5];Kz=10^tmp$par[6]
      basin_char=c(tmp.area,Lb,X1,Wb,POR)
      gw_hyd=c(ALPHA,BETA,Ks,Kb,Kz)}
#######################################################################################################
#FINAL RUN
#######################################################################################################
    Error=bfs(tmp.q,dys,timestep='day',error_basis='base',basin_char,gw_hyd,flow)
    tmp.bf=bfs_out$Baseflow.L3
    tmp.ov=bfs_out$Baseflow.L3>bfs_out$Qob.L3
    tmp.ov[is.na(tmp.ov)]=FALSE
    tmp.bf[tmp.ov]=bfs_out$Qob.L3[tmp.ov]
    BFF=sum(tmp.bf[!is.na(bfs_out$Qob.L3)])/sum(bfs_out$Qob.L3[!is.na(bfs_out$Qob.L3)])

    bf_params=data.frame(tmp.site,tmp.area,Lb,X1,Wb,POR,ALPHA,BETA,Ks,Kb,Kz,Qthresh,Rs,Rb1,Rb2,Prec,Frac4Rise,Error,BFF)
    bf_params[,-1]=signif(bf_params[,-1],6)

#GENERATE CREDIBLE INTERVAL TABLE    
    tmp.ci=bf_ci(bfs_out)
    ci_table=tmp.ci[[1]]
#######################################################################################################
#CALCULATE THE COMPONENTS AS FRACTIONS OF STREAMFLOW
    bff=data.frame(array(NA,dim=c(1,6)))
    dimnames(bff)[[2]]=c('tmp.site','Qmean','BFF','SFF','DRF','Error')
    bff$tmp.site=tmp.site

    tmp=(bfs_out$Qob.L3>0)
    tmp[is.na(tmp)]=FALSE
    tmp[is.na(bfs_out$Qsim.L3)]=FALSE #LIMIT FLOW COMPONENT CALCULATIONS TO TIME STEPS WITH POSITIVE FLOW (FOR TIDAL SYSTEMS)
    bfs_out=bfs_out[tmp,]

    bff$Qmean=signif(mean(bfs_out$Qob.L3),6)

#DAILY BASEFLOW AS FRACTION OF OBSERVED STREAMFLOW
    tmp.bf=bfs_out$Baseflow.L3/bfs_out$Qob.L3
    tmp.bf[tmp.bf>1]=1

#DAILY SURFACE FLOW AS FRACTION OF OBSERVED STREAMFLOW, LIMITED TO 1 - baseflow fraction
    tmp.sf=bfs_out$SurfaceFlow.L3/bfs_out$Qob.L3
    tmp.sf[(tmp.sf+tmp.bf)>0]=1-tmp.bf[(tmp.sf+tmp.bf)>0]

#MEAN BASEFLOW, SURFACE FLOW, AND DIRECT RUNOFF FRACTIONS
    bff$BFF=round(sum(tmp.bf*bfs_out$Qob.L3)/sum(bfs_out$Qob.L3),3) #BASEFLOW FRACTION
    bff$SFF=round(sum(tmp.sf*bfs_out$Qob.L3)/sum(bfs_out$Qob.L3),3) #BASEFLOW
    bff$DRF=round(sum(bfs_out$DirectRunoff.L3)/sum(bfs_out$Qob.L3),3) #DIRECT-RUNOFF FRACTION
    bff$Error=Error} #CLOSE FLOW METRIC CHECK
#######################################################################################################

  list(bf_params,bff,ci_table,bfs_out)}