
#dir_home='C:/cpk/'
#dir_home='~/Projects-active/'

#WORKING DIRECTORY FOR CALIBRATION
#setwd(paste(dir_home,'Baseflow/calibration_20250314',sep=''))
dir_out=paste(getwd(),'/out/',sep='')


#COMPILE PARAMETERS AND BFF
bf_params=data.frame(array(dim=c(0,19)))

bff=data.frame(array(,dim=c(0,6)))

tmp.files=dir(paste(dir_out,'site_sum',sep=''))
  
sites=substr(dir(paste(dir_out,'hydrographs',sep='')),5,25)
  
sites=gsub('.pdf','',sites)

ns=length(sites)  

for (s in 1:ns) {sitefile=paste('bf_params_',sites[s],'.csv',sep='')

  if(sitefile %in% tmp.files) {

    tmp=read.csv(paste(dir_out,'site_sum/',sitefile,sep=''),colClasses=c('character',rep('numeric',18)))
    bf_params=rbind(bf_params,tmp)

    tmp=read.csv(paste(dir_out,'site_sum/bff_',sites[s],'.csv',sep=''),colClasses=c('character',rep('numeric',5)))
    bff=rbind(bff,tmp)}}

dimnames(bf_params)[[2]]=c('SiteID','AREA','Lb','X1','Wb','POR','ALPHA','BETA','Ks','Kb','Kz','Qthresh','Rs','Rb1','Rb2','Prec','Frac4Rise','Error','BFF')
write.csv(bf_params,paste(dir_out,'bf_params.csv',sep=''),row.names=F)
  
dimnames(bff)[[2]]=c('SiteID','Qmean','BFF','SFF','DRF','Error.Total')
write.csv(bff,paste(dir_out,'bff.csv',sep=''),row.names=F)


