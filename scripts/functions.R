##########################################
# Function running the breakpoint models #
##########################################
# Function to  estimated an optimal breakpoint in a linear breakpoint model 
# using segmented package 
breakpoint.con<-function(y=community,x=time){
require(segmented)
seg.mod<-segmented(lm(y~x), ~x, psi=median(x)) # estimates the breakpoint
bp<-confint.segmented(seg.mod)[1] # estimated breakpoint 
bp95CI<-confint.segmented(seg.mod) # confidence intervals
mod<-lm(y~x+I((x-bp)*ifelse(x>bp,1,0))) # alternative implementation
list(mod=mod,bp=bp,seg.mod=seg.mod,bp95CI=bp95CI)
}
#####
# Function to analyse change in pollen community turnover with human arrival
# Pollen community data are summarised in the first axis of a DCA
# and then related to time using a linear breakpoint model with human arrival as breakpoint
# see manuscript for more details
break.point<-function(community, time, bp,name,plot="no"){
require("vegan")
require("MuMIn")
require("grDevices")
#
# Extract fist axis of a DCA
scale <-scores(decorana(community))[,"DCA1"]
#
# change axis so that the change towars present is always positive (irrelevant for analysis, but important for visualisation):
if(mean(scale[time> bp])>mean(scale[time<bp])){scale <-scale*-1} 
if(name %in% c("Mauritius","Flores")){scale <-scale*-1} # to do
#
models<-list()
b1 <- function(x, bp) ifelse(x < bp, bp - x, 0)
b2 <- function(x, bp) ifelse(x < bp, 0, x - bp)
models[["null"]] <- lm(scale ~ 1)
models[["linear"]] <- lm(scale ~ time)
models[["continuous_BP"]] <- lm(scale ~ b1(time, bp) + b2(time, bp))
models[["continuous_BP_flexible"]] <-breakpoint.con(scale,time)
#
# Extract model information:
res<-data.frame(name=name,
bp=bp,
bp.est=models[["continuous_BP_flexible"]]$bp,
bp.95CI.low=models[["continuous_BP_flexible"]]$bp95CI[2],
bp.95CI.up=models[["continuous_BP_flexible"]]$bp95CI[3],
slope_before_BP=NA,
slope_after_BP=NA,
AICc_null=AICc(models[["null"]]),
AICc_linear=AICc(models[["linear"]]),
AICc_continuous_BP=AICc(models[["continuous_BP"]]),
AICc_continuous_BP_flexible =AICc(models[["continuous_BP_flexible"]]$seg.mod),
R2_linear=summary(models[["linear"]])$adj.,
R2_continuous_BP=summary(models[["continuous_BP"]])$adj., 
R2_continuous_BP_flexible=summary(models[["continuous_BP_flexible"]]$mod)$adj.)
res[,c("slope_before_BP","slope_after_BP")]<-diff(predict(models[["continuous_BP"]],newdata=data.frame(time=seq(bp+200,bp-200,length.out=5)[-3])))[-2]/100
#
# Visualistaion (not needed)
if(plot=="yes"){
quartz(,3.5,3.1)
par(mfrow=c(1,1),cex=1) # c(bottom, left, top, right)
plot(scale ~ time,ylab="System state",xlab="",type="p",pch=16,cex=0.8,xlim=c(5000,0),ylim=c(min(scale)-diff(range(scale))*0.1,max(scale)+diff(range(scale))*0.1),main=name)
prd<-predict(models$continuous_BP_flexible$mod, interval="confidence", level = 0.95)
polygon(c(rev(time), time), c(rev(prd[ ,3]), prd[ ,2]), col = 'orange', border = NA)
lines(time,prd[,2],col="steelblue4",lty=2)
lines(time,prd[,3],col="steelblue4",lty=2) 
abline(v=models$continuous_BP_flexible$bp,lty=2,lwd=2,col="orange")
prd<-predict(models[["continuous_BP"]], interval="confidence", level = 0.95)
polygon(c(rev(time), time), c(rev(prd[ ,3]), prd[ ,2]), col = adjustcolor( "lightblue", alpha.f = 0.85), border = NA)
lines(time,prd[,2],col="steelblue4",lty=2)
lines(time,prd[,3],col="steelblue4",lty=2) 
points(prd[,1]~ time,col="steelblue4",type="l",lwd=2)
points(scale ~ time,type="p",pch=16,cex=0.8,col="black")
abline(v=bp,lty=2,lwd=2)}
res} # end of function
##########################

# Function to run breakpoint models for figure 1
run.models<-function(scale,time,bp){
models<-list()
b1 <- function(x, bp) ifelse(x < bp, bp - x, 0)
b2 <- function(x, bp) ifelse(x < bp, 0, x - bp)
models[["continuous_BP"]] <- lm(scale ~ b1(time, bp) + b2(time, bp))
seg.mod<-segmented(lm(scale ~ time), ~ time, psi=median(time)) # estimates the breakpoint
bp<-confint.segmented(seg.mod)[1] # estimated breakpointz 
bp95CI<-confint.segmented(seg.mod) # confidence intervals
mod<-lm(scale ~ time +I((time-bp)*ifelse(time>bp,1,0))) # alternative implementation
models[["continuous_BP_flexible"]] <-list(mod=mod,bp=bp,seg.mod=seg.mod,bp95CI=bp95CI)
models}

########
# Functions to provide "*" dependent on p-values
star<-function(x){
a<-""
if(!is.na(x)){
if(x<0.05){a<-"*"}	
if(x<0.01){a<-"**"}
if(x<0.001){a<-"***"}	
}
a
}
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}
########

# Function tio calculate distance between scores of the first DCA axis before and after a breakpoint
dists<-function(spec, time, bp,name){
require("vegan")
spec <-spec[time<5000,]
spec<-spec[,colSums(spec)>0]
time<-time[time<5000]
dca<-decorana(spec)
dca1.iner<-dca $ evals[1]
scale <-scores(dca)[,"DCA1"]
rm(dca)
if(mean(scale[time>env[i,"human.arrival"]])>mean(scale[time<env[i,"human.arrival"]])){scale <-scale*-1} 
if(env$island[i] %in% c("Mauritius","Maui")){scale <-scale*-1}
res<-c(
mean(abs(dist(scale[time>bp]))), # before human arrival
mean(abs(as.matrix(dist(scale))[time<bp,time>bp])), #before to after
mean(abs(dist(scale[time<bp]))) # after human arrival
)
names(res)<-c("before","between","after")
res
}