rm(list=ls())
####################################################
#############
# load data 
setwd("~/Documents/biodiversity-changes-on-islands")
######
# Information per island. See Supplemental material for more explanations on the used variables.
env<-read.csv2("data/island_variables/env.csv")
#######
#######
# Pollen cores for each island #
	# Pollen community data are stored in the list "specs"
	# Age estimates are stored in the list "ages"
specs<-list()
ages<-list()
for(i in env$island){
	temp<-read.csv2(paste("data/pollen/",i,".csv",sep=""))
	ages[[i]]<-temp[,2]
	specs[[i]]<-temp[,3:ncol(temp)]
	specs[[i]]<-specs[[i]][,colSums(specs[[i]])>0]
rm(temp)}
rm(i)
#######
#######
# Digital elevation model as backround for figures (fig. 1 and S1)
# load raster data for backround map
require("raster")
elevation<-raster("data/elevation/elevation.small.tif")
proj4string(elevation) <- CRS("+proj=longlat +datum=WGS84")
elevation<-extend(elevation,c(-250,215,-80,90))
#######
# end load data
#############
####################################################



####################################################
#############
# load functions
source('scripts/functions.R', chdir = TRUE)
#######
#######
# Test breakpoint models for any island
test<-"Tenerife"
break.point(community=specs[[test]], time=ages[[test]],bp=env[env$island==test,"human.arrival"],name=test,plot="yes")	
#############
####################################################


####################################################
#############
# Run breakpoint models for all islands
# results are stored in the file "modelinfo"
modelinfo<-list()
for(i in 1:nrow(env)){	
spec <-specs[[i]]
time<-ages[[i]]
modelinfo[[as.character(env$island[i])]]<-break.point(spec,time,env[i,"human.arrival"],env[i,"island"])
 # add   plot="yes"	 for visualisation
# dev.off()
}
modelinfo <-do.call(rbind.data.frame, modelinfo)
modelinfo$slope.diff<-apply(abs(modelinfo[,c("slope_before_BP", "slope_after_BP")]),1,diff)
modelinfo$Island<-env$island # name of the island
modelinfo$lat<-env$lat # latitude
modelinfo$lon<-env$lon # longitude
modelinfo$ele<-env$ele # elevation
modelinfo$area<-env$area # Island area
modelinfo$MED <-env$MED # median sea level stand over the previous nine glacial-interglacial cycles (Norder et al. 2019).
modelinfo$Dist <-env$Dist # Distance to next island (Weigelt et al.)
modelinfo$SLMP <-env$SLMP  # Surrounding landmass area (Weigelt et al.)
# write.csv2(modelinfo,"modelinfo.csv")
#############
####################################################


####################################################
##################
# Results #
sum(0<modelinfo$slope.diff) # Number of islands where human arrival caused an increase in turnover rate
#####
#####
# Slope before and after human arrival
apply(abs(modelinfo[,c("slope_before_BP", "slope_after_BP")]),2,median) # median, mean, sd
# paired t-test for systematic difference in slope #
t.test(abs(modelinfo[,c("slope_before_BP")]),abs(modelinfo[,c("slope_after_BP")]),paired=T)
#####
#####
# Rates of pollen compositional turnover increase following human arrival by:
median(abs(modelinfo$slope_after_BP)/abs(modelinfo$slope_before_BP))
mean(abs(modelinfo$slope_after_BP)/abs(modelinfo$slope_before_BP))
sd(abs(modelinfo$slope_after_BP)/abs(modelinfo$slope_before_BP))
#####
#####
# R2-values for different models:
mean(modelinfo$R2_continuous_BP) # with human arrival as breakpoint
mean(modelinfo$R2_linear) # simple linear model
mean(modelinfo$R2_continuous_BP_flexible) # with flexible breakpoint
#
sd(modelinfo$R2_continuous_BP)
sd(modelinfo$R2_linear)
sd(modelinfo$R2_continuous_BP_flexible)
#
median(modelinfo$R2_continuous_BP)
median(modelinfo$R2_linear)
median(modelinfo$R2_continuous_BP_flexible)
#
median(modelinfo$R2_continuous_BP/modelinfo$R2_linear)
#####
#####
# AICc difference larger than 2
sum(((modelinfo$AICc_continuous_BP - modelinfo$AICc_linear) +2)<=0) # in so many cases, the BP model is better than the linear one (AIC diff 2 or larger)
#####
#####
#Number of islands were human arrival falls within the 95CI of optimised breakpoints
sum(modelinfo$bp>modelinfo$bp.95CI.low & modelinfo$bp<modelinfo$bp.95CI.up) # 41%
#
# Number of islands were difference between human arrival estimate and optimised breakpoint is smaller than 500 years
sum(abs(modelinfo$bp-modelinfo$bp.est)<500) # 70.4%
median(abs(modelinfo$bp-modelinfo$bp.est)) # median difference
dif<-modelinfo$bp-modelinfo$bp.est
t.test(dif,mu=0) # no tendency of the human arrival to be earlier or later than estimated breakpoint
#####
#####

####################################################
############
# Figure 1 #
require("sp")
require("colorspace") 
require("lattice")
require("gridBase")
require("grid")
require("vegan")
require("grDevices")
require("segmented")

##########################################################################################
# coordinates defining were the plot for each island is placed on the map
env[env$island =="Hispaniola",c("x","y")]<-c(0.12,0.995) # Hispaniola
env[env$island =="San Cristobal",c("x","y")]<-c(0.09,0.83) # Galapagos
env[env$island =="Maui",c("x","y")]<-c(0.06,0.665) # Maui
env[env$island =="Upolu",c("x","y")]<-c(0.03,0.50) # Upolu
env[env$island =="Raivavae",c("x","y")]<-c(0.03,0.33) # Raivavae (ganz oben)
env[env$island =="Moorea",c("x","y")]<-c(0.03,0.17) # Moorea (ganz links)
env[env$island =="Rimatara",c("x","y")]<-c(0.165,0.17) # Rimatara (mitte)
env[env$island =="Rapa Iti",c("x","y")]<-c(0.20,0.33) # Rapa Iti (ganz unten)
env[env$island =="Mauritius",c("x","y")]<-c(0.64,0.27) #00.49,0.32 
env[env$island =="Iceland",c("x","y")]<-c(0.49,0.94) # Iceland
env[env$island =="Tenerife",c("x","y")]<-c(0.50,0.745) # 
env[env$island =="Gran Canaria",c("x","y")]<-c(0.68,0.91) # 
env[env$island =="La Gomera",c("x","y")]<-c(0.625,0.745) # 
env[env$island =="Pico",c("x","y")]<-c(0.30,0.97) # 
env[env$island =="Flores",c("x","y")]<-c(0.325,0.805) #
env[env$island =="Great Mercury",c("x","y")]<-c(0.88,0.335) # Mercury Island
env[env$island =="Tawhiti Rahi",c("x","y")]<-c(0.80,0.17) # Tawhiti Ravi
env[env$island =="New Caledonia",c("x","y")]<-c(0.717,0.43) #0.65,0.49
env[env$island =="Taveuni",c("x","y")]<-c(0.89,0.53) #
env[env$island =="Viti Levu",c("x","y")]<-c(0.88,0.695) #
env[env$island =="Yacata",c("x","y")]<-c(0.825,0.86) # Yacata
env[env$island =="Santo Antao",c("x","y")]<-c(0.475,0.52) # Cap Verde
env[env$island =="Sao Nicolau",c("x","y")]<-c(0.60,0.57)#0.485,0.375         
env[env$island =="Robinson Crusoe",c("x","y")]<-c(0.35,0.33) # 0.19,0.26
env[env$island =="Alexander Selkirk",c("x","y")]<-c(0.31,0.17) #
env[env$island =="Tristan da Cunha",c("x","y")]<-c(0.50,0.33) #
env[env$island =="Nightingale Island",c("x","y")]<-c(0.46,0.17) # 

# Plot the map
quartz(, width = 12, height =6.5)
par(mar=c(0, 0, 0, 0),bty="n")# par(mar=c(5, 4, 4, 2))
plot(elevation,col=rev(gray.colors(20,alpha=1)),legend=F,xaxt="n",yaxt="n")
lines(x=c(-25.2,-14),y=c(17.1,5),lwd=2,lty=3)
lines(x=c(-25.2,44),y=c(17.2,21),lwd=2,lty=3)
lines(x=c(-70.924,-135),y=c(19.028,105),lwd=2,lty=3)
lines(x=c(-89.48004,-150),y=c(-0.895295,60),lwd=2,lty=3)
lines(x=c(-20.6,13),y=c(64.1,88),lwd=2,lty=3) 
lines(x=c(-30,-34),y=c(39.4,75),lwd=2,lty=3)
lines(x=c(-16.4,0),y=c(28.5,39),lwd=2,lty=3)
lines(x=c(174.7,152),y=c(-35.5,-82.5),lwd=2,lty=3)
lines(x=c(175,190),y=c(-37,-40),lwd=2,lty=3)
lines(x=c(-13,1),y=c(-37,-40),lwd=2,lty=3)
lines(x=c(-13,-15),y=c(-37,-83),lwd=2,lty=3)
lines(x=c(-79,-70),y=c(-33,-40),lwd=2,lty=3)
lines(x=c(-80,-93),y=c(-33,-83),lwd=2,lty=3)
lines(x=c(58.5,70),y=c(-21,-59),lwd=2,lty=3)
lines(x=c(-157,-166),y=c(21,30),lwd=2,lty=3)
lines(x=c(-174,-182),y=c(-13,-12),lwd=2,lty=3)
lines(x=c(-144,-138),y=c(-29,-40),lwd=2,lty=3)
lines(x=c(-148,-167),y=c(-25,-83),lwd=2,lty=3)
lines(x=c(-153,-183),y=c(-24,-83),lwd=2,lty=3)
lines(x=c(-150,-180),y=c(-17,-42),lwd=2,lty=3)
lines(x=c(177.5,164),y=c(-15,65),lwd=2,lty=3)
lines(x=c(180,190),y=c(-15,-12),lwd=2,lty=3)
points(lat~ long,data=env,col="orange",bg="black",pch=21,cex=1.1,lwd=2)

##########################################################################################
# show<-"breakpoint and human" will build figure 1
# show<-"boxplot" will build figure S1
show<-"breakpoint and human" #  "boxplot" ; "breakpoint and human" ;
if(show=="boxplot"){
for(i in 1:nrow(env)){
spec<-specs[[i]]
time<-ages[[i]]
spec <-spec[time<5000 & time>-50,]
time<-time[time<5000 & time>-50]
scale <-scores(decorana(spec))[,"DCA1"]	
if(mean(scale[time>env[i,"human.arrival"]])>mean(scale[time<env[i,"human.arrival"]])){scale <-scale*-1} 
if(env$island[i] %in% c("Mauritius","Maui")){scale <-scale*-1}
before.human<-scale[time>env[i,"human.arrival"]]
after.human<-scale[time<env[i,"human.arrival"]]
pushViewport(viewport(x=env[i,"x"],y=env[i,"y"],width=.08,height=.12,just=c("left","top")))
grid.rect()
par(plt = gridPLT(), new=TRUE)	
boxplot(before.human,after.human,col=c("snow3","firebrick"),names=c("",""),range=0,ylim=c(min(scale),max(scale)+diff(range(scale))/4))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = rgb(1, 1, 1, 0.5, names = NULL, maxColorValue = 1)) 
boxplot(before.human,after.human,col=c("snow3","firebrick"),add=T,names=c("",""),range=0,ylim=c(min(scale),max(scale)+diff(range(scale))/4))
mtext(text=env[i,"island"],side=3,line=-0.9,cex=1,adj=0,at=0.5)
popViewport( )
}}

if(show=="breakpoint and human"){
for(i in 1:nrow(env)){
spec<-specs[[i]]
time<-ages[[i]]
spec <-spec[time<5000 & time>-50,]
time<-time[time<5000 & time>-50]
scale <-scores(decorana(spec))[,"DCA1"]
# change axis so that the change towars present is always positive:
if(mean(scale[time> env[i,"human.arrival"]])>mean(scale[time<env[i,"human.arrival"]])){scale <-scale*-1} 
if(env$island[i] %in% c("Mauritius","Maui")){scale <-scale*-1}
models<-run.models(scale=scale,time=time,bp=env[i,"human.arrival"])
pushViewport(viewport(x=env[i,"x"],y=env[i,"y"],width=.105,height=.12,just=c("left","top")))
grid.rect()
par(plt = gridPLT(), new=TRUE)	
plot(scale ~ time,ylab="",xlab="",type="p",pch=16,cex=0.5,xlim=c(5000,0),ylim=c(min(scale)-diff(range(scale))*0.1,max(scale)+diff(range(scale))*0.1),mgp = c(1.6, 0.3, 0),cex.axis=0.7, tcl =-0.3,col="black")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = rgb(1, 1, 1, 0.5, names = NULL, maxColorValue = 1)) 
abline(v=models$continuous_BP_flexible$bp,lty=3,lwd=1,col="black")
abline(v=env[i,"human.arrival"],lty=1,lwd=2,col="orange")
abline(v=env[i,"human.arrival"],lty=3)
prd<-predict(models[["continuous_BP"]], interval="confidence", level = 0.95)
polygon(c(rev(time), time), c(rev(prd[ ,3]), prd[ ,2]), col = adjustcolor( "lightblue", alpha.f = 0.85), border = NA)
lines(time,prd[,2],col="steelblue4",lty=2)
lines(time,prd[,3],col="steelblue4",lty=2) 
points(prd[,1]~ time,col="steelblue4",type="l",lwd=2)
points(scale~ time,type="p",pch=16,cex=0.5,col="black")
if (env$ID[i] %in% c("Nightingale Island Ljung 2007")){
mtext(text=env[i,"island"],side=1,line=-1,cex=0.7,adj=0,at=4800)	
}
else
{mtext(text=env[i,"island"],side=3,line=-1,cex=0.7,adj=0,at=4800)}
popViewport( )	
}}
############
####################################################

####################################################
############
# Figure 2 #
quartz(,3,4 )
boxplot(abs(modelinfo[,c("slope_before_BP", "slope_after_BP")]),range=0,ylab="Rate of compositional turnover",names=c("pre-human ","human"),col=c("steelblue4","orange"))
############
####################################################

####################################################
############
# Figure 3 #
quartz(,14,7)
par(mfrow=c(2,4),cex=1)
plot(slope.diff~bp,dat= modelinfo,ylab="Difference in turnover rate with human arrival",xlab="Human arrival [years BP]",pch=21,col="black",bg="orange")
mod<-lm((slope.diff)~(bp),dat= modelinfo)
mod<-lm((slope.diff)~log(bp),dat= modelinfo)
summary(mod)
curve(coef(mod)[1]+coef(mod)[2]*log(x),from=100,to=3500,add=T,lwd=2)
mtext(text=bquote(italic(R)^2 == .(round(summary(mod)$r.squared, 2)) * .(star(lmp(mod)))),line=-1.5,cex=1,side=3,at=2200)
###
# Rate of change preceding human arrival
plot(slope.diff~slope_before_BP,dat= modelinfo,ylab="",xlab="Turnover rate before humans",pch=21,col="black",bg="orange")
mod<-lm((slope.diff)~slope_before_BP,dat= modelinfo)
summary(mod)
mtext(text=bquote(italic(R)^2 == .(round(summary(mod)$adj.r, 2)) * .(star(lmp(mod)))),line=-1.5,cex=1,side=3,at=0)
curve(coef(mod)[1]+coef(mod)[2]*x,from=-0.003,to=0.002,add=T,lwd=1,lty=2)
###
plot(slope.diff~abs(lat),dat= modelinfo,ylab="",xlab="Latitude [abs °]",pch=21,col="black",bg="orange")
mod<-lm((slope.diff)~abs(lat),dat= modelinfo)
curve(coef(mod)[1]+coef(mod)[2]*x,from=0,to=600,add=T,lwd=1,lty=2)
summary(mod)
mtext(text=bquote(italic(R)^2 == .(round(summary(mod)$r.squared, 2)) * .(star(lmp(mod)))),line=-1.5,cex=1,side=3,at=50)
###
plot(slope.diff~ele,dat= modelinfo,ylab="",xlab="Elevation [m]",pch=21,col="black",bg="orange")
mod<-lm((slope.diff)~ele,dat= modelinfo)
mtext(text=bquote(italic(R)^2 == .(round(summary(mod)$r.squared, 2)) * .(star(lmp(mod)))),line=-1.5,cex=1,side=3,at=2000)
curve(coef(mod)[1]+coef(mod)[2]*x,from=0,to=2500,add=T,lwd=1,lty=2)
###
plot(slope.diff~log10(area),dat= modelinfo,ylab="",xlab="Area (log10 scale)[sqkm]",pch=21,col="black",bg="orange")
mod<-lm((slope.diff)~log10(area) ,dat= modelinfo)
summary(mod)
mtext(text=bquote(italic(R)^2 == .(round(summary(mod)$r.squared, 2)) * .(star(lmp(mod)))),line=-1.5,cex=1,side=3,at=3)
curve(coef(mod)[1]+coef(mod)[2]*x,from=0,to=2500,add=T,lwd=1,lty=2)
###
plot(slope.diff~log10(MED),dat= modelinfo,ylab="",xlab="Glacial-interglacial area (lig10 scale) [sqkm]",pch=21,col="black",bg="orange")
mod<-lm((slope.diff)~ log10(MED),dat= modelinfo)
summary(mod)
curve(coef(mod)[1]+coef(mod)[2]*x,from=0,to=12,add=T,lwd=1,lty=2)
mtext(text=bquote(italic(R)^2 == .(round(summary(mod)$adj.r, 2)) * .(star(lmp(mod)))),line=-1.5,cex=1,side=3,at=4)
###
plot(slope.diff~Dist,dat= modelinfo,ylab="",xlab="Distance mainland [km]",pch=21,col="black",bg="orange")
mod<-lm((slope.diff)~ Dist,dat= modelinfo)
summary(mod)
mtext(text=bquote(italic(R)^2 == .(round(summary(mod)$adj.r, 2)) * .(star(lmp(mod)))),line=-1.5,cex=1,side=3,at=4000)
curve(coef(mod)[1]+coef(mod)[2]*x,from=0,to=6000,add=T,lwd=1,lty=2)
###
plot(slope.diff~SLMP,dat= modelinfo,ylab="",xlab="Sourrounding landmass [sqkm]",pch=21,col="black",bg="orange")
mod<-lm((slope.diff)~ SLMP,dat= modelinfo)
summary(mod)
mtext(text=bquote(italic(R)^2 == .(round(summary(mod)$adj.r, 2)) * .(star(lmp(mod)))),line=-1.5,cex=1,side=3,at=0.6)
curve(coef(mod)[1]+coef(mod)[2]*x,from=0,to=1.4,add=T,lwd=1,lty=2)
############
####################################################

####################################################
#############
# Figure S2 #
# Absolute change in pollen composition with human arrival
# Calculate distance in axis scores of the first DCA axis
# comparing distance before with distance after human arrival
dist.res<-list()
for(i in 1:nrow(env)){	
spec <-specs[[i]]
time<-ages[[i]]
dist.res[[as.character(env$island[i])]]<-dists(spec,time,bp=env[i,"human.arrival"],env[i,"island"])
}
dist.res <-do.call(cbind.data.frame, dist.res)
boxplot(t(dist.res)) # figures S2
round(apply(dist.res,1,median),2)
round(apply(dist.res,1,mean),2)
round(apply(dist.res,1,sd),2)
t.test(t(dist.res)[,2],t(dist.res)[,3])
############
####################################################

###################################################
######################
# Randomisation test #
# Repeat breakpoint models on randomised data
# Community data are randomised, while human arrival and time is left unchanged
b1 <- function(x, bp) ifelse(x < bp, bp - x, 0)
b2 <- function(x, bp) ifelse(x < bp, 0, x - bp)
runs<-list()
for(r in 1:1000){
randomodel<-list() # estimated slopes with random data
for(i in 1:nrow(env)){	
spec <-specs[[i]]
time<-ages[[i]]
spec<-spec[,colSums(spec)>0]
spec<-spec[sample(1:length(time),length(time)),] # random order
scale <-scores(decorana(spec))[,"DCA1"]
mod <- lm(scale ~ b1(time, env[i,"human.arrival"]) + b2(time, env[i,"human.arrival"]))
randomodel[[as.character(env$island[i])]]<-diff(predict(mod,newdata=data.frame(time=seq(env[i,"human.arrival"]+200,env[i,"human.arrival"]-200,length.out=5)[-3])))[-2]/100
names(randomodel[[as.character(env$island[i])]])<-c("slope_before_BP","slope_after_BP")
}
runs[[r]] <-as.data.frame(t(as.data.frame(randomodel)))
rm(randomodel)
cat("\r", r, "of", 1000) 
flush.console()  
}
rm(b1,b2)
######################
# Number of cases in which the true measured distance in turnover rate from before to after human arrival is higher than that in any of the 1000 runs with randomized community data
sum(unlist(lapply(runs,function(x){mean(abs(x$slope_after_BP)-abs(x$slope_before_BP))}))<mean(apply(abs(modelinfo[,c("slope_before_BP", "slope_after_BP")]),1,diff)))
######################
###################################################