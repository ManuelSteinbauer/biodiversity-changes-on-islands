# The script simulates pollen community data with a long term gradient in composition as well as random short term nois.
# The data are use to show that Ordination based analyses are able to detect this long term gradient


#############
# Functions #
#############
rate.of.change.with.dca<-function(community,timesteps){
dca<-scores(decorana(community[timesteps,]))[,"DCA1"]
if(dca[length(dca)]<dca[1]){dca <-dca*-1}
coef(lm(dca~ timesteps))[2]}

rate.of.change.with.similarity <-function(community,timesteps){
results<-data.frame(Rate.of.change=rep(NA, length(timesteps)-1),Mean.time.between.samples=NA)
community.s<-decostand(community,method="chi.square")
for(i in 1:(length(timesteps)-1)){
results[i,"Rate.of.change"]<-vegdist(rbind(community.s[timesteps[i],], community.s[timesteps[i+1],]),method="euclidean")/(timesteps[i+1]-timesteps[i])
results[i,"Mean.time.between.samples"]<-i+(timesteps[i+1]-timesteps[i])}
mean(results$Rate.of.change)}
#############
#############

############################
# Simualted community data #
############################
specs.with.variability<-data.frame(spec1=1:100,spec2=100:1,spec3=50+sample(-50:+50,100)) 
# three species with:
# one changing abundance from 1 to 100
# one changing abundance from 100 to 1
# one with abundance of 50 + Variability

specs.without.variability<-data.frame(spec1=1:100,spec2=100:1,spec3=50) 
# three species with:
# one changing abundance from 1 to 100
# one changing abundance from 100 to 1
# one with abundance of 50

############################
############################
library(vegan)
# Rate-of-change quantified with dissimilarity between neighbouring temporal units
rate.of.change.with.similarity(community = specs.without.variability, timesteps=1:100)

# even sampling:
rate.of.change.with.similarity(community = specs.without.variability, timesteps=seq(1,100,20)) # here for distance 20 (n=5 samples)


# Rate-of-change with variability: 
result<-rep(NA,100) # 100 simulations with independent data
for(r in 1:100){
specs.with.variability<-data.frame(spec1=1:100,spec2=100:1,spec3=50+sample(-50:+50,100)) 
result[r]<-rate.of.change.with.similarity(community = specs.with.variability, timesteps=1:100)
}
mean(result)
sd(result)

for(r in 1:100){
specs.with.variability<-data.frame(spec1=1:100,spec2=100:1,spec3=50+sample(-50:+50,100)) 
result[r]<-rate.of.change.with.similarity(community = specs.with.variability, timesteps=seq(1,100,20)) # here for distance 20 (n=5 samples)
}
mean(result)
sd(result)

# uneven sampling:
result<-rep(NA,1000)
for(r in 1:1000){
result[r]<-rate.of.change.with.similarity(community = specs.without.variability, timesteps=sort(sample(1:100,5))) # here for 5 samples
}
mean(result)
sd(result)

for(r in 1:1000){
specs.with.variability<-data.frame(spec1=1:100,spec2=100:1,spec3=50+sample(-50:+50,100)) 
result[r]<-rate.of.change.with.similarity(community = specs.with.variability, sort(sample(1:100,100))) # here for 5 samples

}
mean(result)
sd(result)



#########
# DCA
# Rate-of-change measured in DCA units
rate.of.change.with.dca(community = specs.without.variability, timesteps=1:100)
rate.of.change.with.dca(community = specs.without.variability, timesteps=seq(1,100,20))

# with random varibility
result<-rep(NA,100)
for(r in 1:100){
specs.with.variability<-data.frame(spec1=1:100,spec2=100:1,spec3=50+sample(-50:+50,100)) 
result[r]<-rate.of.change.with.dca(community = specs.with.variability, timesteps=seq(1,100,20)) # even sampling with distance 20 (n=5)
}
mean(result)
sd(result)

result<-rep(NA,100)
for(r in 1:100){
specs.with.variability<-data.frame(spec1=1:100,spec2=100:1,spec3=50+sample(-50:+50,100)) 
result[r]<-rate.of.change.with.dca(community = specs.with.variability, timesteps=sort(sample(1:100,20)))# uneven sampling of 20 samples
}
}
mean(result)
sd(result)

# Visualisation of DCA results #

par(mfrow=c(2,2))
# DCA sampling all time intervals:
dca<-scores(decorana(specs.without.variability))[,"DCA1"]
if(dca[length(dca)]<dca[1]){dca <-dca*-1}
plot(dca~time)
mod<-lm(dca~time)
abline(mod)
coef(mod)[2] # predicted slope, modelled slope without error is 0.01402

# DCA with uneven sampling
sampled.sites<-sort(sample(1:100,10))
dca<-scores(decorana(specs.without.variability[sampled.sites,]))[,"DCA1"]
if(dca[length(dca)]<dca[1]){dca <-dca*-1}
plot(dca~ sampled.sites,xlim=c(0,100))
mod<-lm(dca~ sampled.sites)
abline(mod)
coef(mod)[2] # predicted slope, modelled slope without error is 0.01402
# -> DCA can handle uneven sampling
###########

# DCA sampling all time intervals:
dca<-scores(decorana(specs.with.variability))[,"DCA1"]
if(dca[length(dca)]<dca[1]){dca <-dca*-1}
plot(dca~time)
mod<-lm(dca~time)
abline(mod)
coef(mod)[2] # predicted slope, modelled slope without error is 0.01402

# DCA with uneven sampling
sampled.sites<-sort(sample(1:100,10))
dca<-scores(decorana(specs.with.variability[sampled.sites,]))[,"DCA1"]
if(dca[length(dca)]<dca[1]){dca <-dca*-1}
plot(dca~ sampled.sites,xlim=c(0,100))
mod<-lm(dca~ sampled.sites)
abline(mod)
coef(mod)[2] # predicted slope, modelled slope without error is 0.01402
# -> DCA can handle uneven sampling
###########
