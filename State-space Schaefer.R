
setwd( "C:/Users/James.Thorson/Desktop/Project_git/state_space_production_model/" )
library( TMB )

# Parameters
r = 0.5 
k = 1000
sigmap = 0.05
sigmam = 0.1
n_years = 100
n_rep = 100
overdispersion = 2  # 1=No overdispersion

#
Date = Sys.Date()
  DateDir = paste0(getwd(),"/",Date,"/")
  dir.create(DateDir)

# Compile results
Results = array(NA, dim=c(2,n_rep,3), dimnames=list(c("without","with"),paste0("Rep_",1:n_rep),c("k","RE_final_depletion","sigmam")))

# Loop
for(repI in 1:n_rep){
  set.seed(repI)

  # exploitation fraction
  exploit_max = r * 0.6
  exploit_t = c( rep(0.01,ceiling(n_years/4)),seq(0.01,exploit_max,length=ceiling(n_years/4)),rep(exploit_max,floor(n_years/4)),seq(exploit_max,r/4,length=floor(n_years/4)) )
  
  # Simulate
  catch_t = n_t = rep(NA,n_years-1)
  n_t[1] = rnorm(1, k, sigmap)
  for(t in 2:n_years){
    catch_t[t-1] = n_t[t-1] * exploit_t[t-1]
    n_t[t] = (n_t[t-1]-catch_t[t-1]) * (1 + r*(1-(n_t[t-1]-catch_t[t-1])/k)) + rnorm(1, mean=0, sd=sigmap*n_t[t-1])
    if( n_t[t]< k/100 ) n_t[t] = k/100
  }
  nobs_t = n_t + rnorm(n_years, mean=0, sd=sigmam*n_t)
  
  # Visualize
  png( paste0(DateDir,"Rep_",repI,"_true_dynamics.png"), width=10, height=7.5, res=200, units="in")
    par(mfrow=c(1,2))
    matplot( cbind(n_t,nobs_t), type="l", ylim=c(0,k*1.5))
    plot(catch_t)
  dev.off()
  
  #################
  # TMB prep
  #################
  
  # Compile
  Version_statespace = "statespace_v1"
  compile( paste0(Version_statespace,".cpp") )
  
  # Fit the model
  Data = list("y_t"=nobs_t, "c_t"=catch_t, "ysd_t"=rep(sigmam/overdispersion,n_years), "penalties_z"=c(log(0.1),log(1)) )
  Params = list("log_r"=log(r), "log_k"=log(max(nobs_t)), "log_q"=log(1), "log_sigmap"=log(1), "log_sigmam"=log(1e-10), "log_sigmac"=log(0.01), "log_x_t"=log(nobs_t), "logit_exploit_t"=rep(0,n_years-1) )
  Random = c("log_x_t", "logit_exploit_t")
  
  # Save plot
  Plot_Fn = function( report, sdsummary, dir, name, ... ){
    png( paste0(dir,"/",name), width=10, height=7.5, res=200, units="in")
      par( mfrow=c(2,2), mar=c(3,3,2,0), mgp=c(2,0.5,0), tck=-0.02)
      # Biomass
      matplot( cbind(n_t,report$x_t), type="l", ylim=c(0,k*1.5), ylab="", xlab="year", main="Biomass")
      Mat = cbind( report$x_t, sdsummary[which(rownames(sdsummary)=="x_t"),"Std. Error"])
      polygon( x=c(1:n_years,n_years:1), y=c(Mat[,1]+Mat[,2],rev(Mat[,1]-Mat[,2])), col=rgb(1,0,0,0.2), border=NA)
      # Depletion
      matplot( cbind(n_t/k,report$Depletion_t), type="l", ylim=c(0,1.5), ylab="", xlab="year", main="Relative biomass")
      Mat = cbind( report$Depletion_t, sdsummary[which(rownames(sdsummary)=="Depletion_t"),"Std. Error"])
      polygon( x=c(1:n_years,n_years:1), y=c(Mat[,1]+Mat[,2],rev(Mat[,1]-Mat[,2])), col=rgb(1,0,0,0.2), border=NA)
      # Catch
      matplot( cbind(catch_t,report$cpred_t[-n_years]), type="l", log="y", ylab="", xlab="year", main="Catch")
      # Exploitation rate
      matplot( cbind(exploit_t[-n_years],report$exploit_t), type="l", log="y", ylab="", xlab="year", main="Exploitation fraction")
    dev.off()  
  }
  
  ######
  # Run without estimating overdispersion
  ######
  Map = list()
  Map[["log_r"]] = factor(NA)
  Map[["log_sigmac"]] = factor(NA)
  Map[["log_sigmam"]] = factor(NA)
  
  # Compile
  dyn.load( dynlib(Version_statespace) )
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map)
  
  # Optimize
  Opt = nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr, control=list("trace"=1, "eval.max"=1e4, "iter.max"=1e4))
  Opt[["diagnostics"]] = data.frame( "Est"=Opt$par, "final_gradient"=Obj$gr(Opt$par) )
  Report = Obj$report()
  SD = try(sdreport( Obj ))
  
  # visualize and save results
  if( all(abs(Opt$diagnostics$final_gradient)<0.01) ){
    Results["without",repI,c("k","RE_final_depletion","sigmam")] = c(Report$k, Report$Depletion_t[n_years]/(n_t[n_years]/k), exp(Report$log_sigmam))
    if( !inherits(SD,"try-error")) Plot_Fn( report=Report, sdsummary=summary(SD), dir=DateDir, name=paste0("Rep_",repI,"_Without_overdispersion.png"))
  }
  
  ######
  # Run with estimating overdispersion
  ######
  Map = list()
  Map[["log_r"]] = factor(NA)
  Map[["log_sigmac"]] = factor(NA)
  Params[["log_sigmam"]] = log(1)
    
  # Compile
  dyn.load( dynlib(Version_statespace) )
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map)
  
  # Optimize
  Opt = nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr, control=list("trace"=1, "eval.max"=1e4, "iter.max"=1e4))
  Opt[["diagnostics"]] = data.frame( "Est"=Opt$par, "final_gradient"=Obj$gr(Opt$par) )
  Report = Obj$report()
  SD = try(sdreport( Obj ))
  
  # visualize and save results
  if( all(abs(Opt$diagnostics$final_gradient)<0.01) ){
    Results["with",repI,c("k","RE_final_depletion","sigmam")] = c(Report$k, Report$Depletion_t[n_years]/(n_t[n_years]/k), exp(Report$log_sigmam))
    if( !inherits(SD,"try-error")) Plot_Fn( report=Report, sdsummary=summary(SD), dir=DateDir, name=paste0("Rep_",repI,"_With_overdispersion.png"))
  }
}

#####################
# Compile results
#####################

HistFn = function(x0, x1, trueval, nbreaks=10, ...){
  xlim = range(pretty(range( c(trueval,x0,x1),na.rm=TRUE )))
  Hist0 = hist(x0, plot=FALSE, breaks=seq(xlim[1],xlim[2],length=nbreaks))
  Hist1 = hist(x1, plot=FALSE, breaks=seq(xlim[1],xlim[2],length=nbreaks))
  ylim = c(0, max(c(Hist0$counts,Hist1$counts)) )
  plot(Hist0, xlim=xlim, col=rgb(0,0,0,0.2), ylim=ylim, ...)
  plot(Hist1, xlim=xlim, col=rgb(1,0,0,0.2), ylim=ylim, add=TRUE)
  abline(v=c(mean(x0,na.rm=TRUE),mean(x1,na.rm=TRUE),trueval), lwd=3, col=c("black","red","blue"))    
}
png( paste0(DateDir,"_Boxplot_summary.png"), width=10, height=5, res=200, units="in")
  par( mfrow=c(1,3), mar=c(3,3,4,0), mgp=c(2,0.5,0), tck=-0.02, yaxs="i")
  HistFn( x0=Results["with",,"k"], x1=Results["without",,"k"], trueval=k, xlab="", ylab="", main="Carrying capacity", yaxt="n" )
  HistFn( x0=Results["with",,"RE_final_depletion"], x1=Results["without",,"RE_final_depletion"], trueval=1, xlab="", ylab="", main="Error in final\nrelative abundance", yaxt="n" )
  HistFn( x0=Results["with",,"sigmam"], x1=Results["without",,"sigmam"], trueval=sigmam, xlab="", ylab="", main="Overdispersion", yaxt="n" )
dev.off()

