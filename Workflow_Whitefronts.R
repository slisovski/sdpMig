## Install spdMig package (only once)
# library(devtools)
# install_github("slisovski/sdpMig", force = TRUE)
 
 
library(spdMig)

###############################################
### Parameters ################################
###############################################

parms <- list(
  MaxT   = 100,  ## Maximum time
  MaxX   = 100,  ## Maximum body condition
  NSites = 8,    ## Excl. Breeding site
  
  ### Species specific Parameters ##
  
  B0 = 3,       ## Future reproductive success
  w  = 0.028,   ## Parameters for sigmoidal TR function
  xc = 55,      ## 
  max_u = 1.0,  ## maximum foraging intensity
  f = 1.0,      ## proportion of x that cannot be devoted to flying
  
  
  ## Flying capacities
  c     = 14776,
  speed = 1440,
  
  ## Wind
  WindAssist = 0,
  WindProb   = 1,
  
  ## Decision Error
  decError = 8000,
  
  ## Terminal Reward
  xFTReward = c(0, 86, 87, 97, 98, 100),
  yFTReward = c(0, 0,   2,  2,  0,   0),
  
  
  ### Site specific Parameters ###
  path = "g://Dropbox//Modelling_Migration_SDP//Packages//sdpMig//Parameters//WFGparameter_ex4.csv",

  pred_a1  = 2,
  pred_a2  = 2,
  
  ### Accuracy
  ZStdNorm = c(-2.5, -2.0, -1.5, -1.0, -0.5,  0.0,  0.5,  1.0,  1.5,  2.0,  2.5),
  PStdNorm = c(0.0092, 0.0279, 0.0655, 0.1210, 0.1747, 0.2034, 0.1747, 0.1210, 0.0655, 0.0279, 0.0092)
  
) ## End parameter ####


#######################
### Simulation ########
#######################


setwd("g:/Dropbox/Modelling_Migration_SDP/packages/sdpMig/Output/")
write("Scenario,mort,stag1,stag2,stag3,stag4,stag5,stag6,stag7,stag8,stag9,meanX1,meanX2,meanX3,meanX4,meanX5,meanX6,meanX7,meanX8,meanX9,fuelrate1,fuelrate2,fuelrate3,fuelrate4,fuelrate5,fuelrate6,fuelrate7,fuelrate8,fuelrate9",
      file="test_whitefronts.csv", sep=",", append=F)

sims <- c("default","change_B0","change_max_u","change_f","change_flightspeed","change_expend","change_sitefuel")

for (scenario in 1:length(sims)) {

    sdp  <- makeSDPmig(parms, "Whitefronted Geese")
    
    switch(sims[scenario],
           "change_B0" =           {sdp@Species$B0           <- sdp@Species$B0 * 0.9 },   # 2
           "change_max_u" =        {sdp@Species$max_u        <- sdp@Species$max_u * 0.9}, # 3
           "change_f" =            {sdp@Species$f            <- sdp@Species$f * 0.9},     # 4
           "change_flightspeed" =  {sdp@Species$speed        <- sdp@Species$speed * 0.9 },# 5
           "change_expend" =       {sdp@Sites$expend         <- sdp@Sites$expend * 1.1 }, # 6
           "change_sitefuel" =     {sdp@Sites$gain$gain.y    <- sdp@Sites$gain$gain.y * 0.9}  # 7
           )

    sdpM <- bwdIteration(sdp, pbar = TRUE)
    simu <- MigSim(sdpM, 100, 1, 1, c(25, 35)) #MigSim(object, NrInd, starttime, startsite, start_x(min, max))
    smP  <- simuPlot(simu, sdpM, fun = "mean")
    decisionPlot(sdpM)
    
    
    # Calculate three major output measures: mortality, staging times per site, body condition per site, fuelling rates per site 
    # Mortality
    mort <- sum(apply(smP, 1, function(x) any(is.na(x))))/nrow(smP); mort 
    
    # Staging times per site
    smP # individual staging times
    smP[apply(smP, 1, function(x) any(is.na(x))),] <- cbind(rep(NA, ncol(smP)))
    mean.staging <- apply(smP, 2, mean, na.rm = T); mean.staging
    
    # Mean body condition per site over time
    meanX <- matrix(ncol = sdp@Init$MaxT, nrow = sdp@Init$NSites+1)
    for(i in 1:ncol(meanX)) {
      t1 <- aggregate(simu[,3,i], by = list(simu[,2,i]), FUN = function(x) mean(x, na.rm = T))
      meanX[t1$Group.1+1, i]  <- t1$x
    }
    t(meanX)
    mean.bodycond <- apply(t(meanX), 2, weighted.mean, w= , na.rm=T); mean.bodycond
    
    ### Fuelling rates per site
    cond <- array(dim = c(sdp@Init$NSites+1, 4, dim(simu)[1]))
    for(i in 1:dim(simu)[1]) {
      if(all(simu[i,5,]==0)) {
        tmp  <- simu[i,,-which(is.na(simu[i,4,]) & simu[i,2,]<sdp@Init$NSites)]
        tmp2 <- lapply(split(as.data.frame(t(tmp)), f = tmp[2,]), function(x) cbind(nrow(x), x[1,3], x[nrow(x),3]))
        out <- cbind(0:sdp@Init$NSites, NA, NA, NA)
        out[match(as.numeric(names(tmp2)), out[,1]),-1] <- do.call("rbind", tmp2)
        cond[,,i] <- out
      }
    }
    indRate <- t(apply(cond, c(1,3), function(x) (x[4]-x[3])/x[2]))
      indRate[,sdp@Init$NSites+1] <-  apply(cond, 3, function(x) x[nrow(x),3])
    mean.fuelrates <- apply(indRate, 2, mean, na.rm=T); mean.fuelrates
    
    #  plot(apply(indRate, 2, median, na.rm = T)[-9]); points(apply(indRate, 2, mean, na.rm = T)[-9], pch = 16)  
      
    write(cbind(sims[scenario], mort, t(mean.staging), t(mean.bodycond), t(mean.fuelrates) ), 
          file="test_whitefronts.csv", sep=",", append=T, ncolumns = (sdp@Init$NSites+1) * 3 + 2)
}












matplot(t(simu[,3,]), type = "l", col = "grey90", lty = 1)
matplot(t(meanX), type = "l", col = rainbow(9), lwd = 2, lty = 1, add = T)












#apply(simu, 1, function(x) aggregate(x[2,], by = list(x[,3], FUN = apply(x[2:3,], 1, sum))
#merge(data.frame(Var1 = 0:sdp@Init$NSites), as.data.frame(table(x[2])), all.x = T)[,2]}) 



######### Calculate individual staging times per site

res <- matrix(NA, ncol = sdpM@Init$NSites+1, nrow = dim(simu)[1])

for(i in 1:nrow(res)) {
  tmp <- simu[i,,]
  if(sum(tmp[5,])>0) tmp <- tmp[,1:(min(which(tmp[5,]==1)))]
  res[i,] <- apply(cbind(0:sdpM@Init$NSites), 1, function(x) {
    ifelse(all(tmp[2,]<x) , NA, sum(tmp[2,]==x, na.rm = T))
  })
}
Staging <- res
Staging[,1] # staging on site 1

plot(  Staging[,]) 




######### Plot Individual body condition over time

plot(NA, xlim = c(0, sdpM@Init$MaxT+1), ylim = c(1, sdpM@Init$MaxX), bty = "n", xlab = "time", ylab = "body condition", 
     las = 1, cex.lab  = 1.2)

for(i in 1:dim(simu)[1]) {
  tmp <- simu[i,,]
  tmp <- tmp[,tmp[5,]!=1]
  if(!is.null(nrow(tmp))) lines(tmp[1,], tmp[3,], col = adjustcolor("grey60", alpha.f = 0.6))
}
for(i in 1:dim(simu)[1]) {
  tmp <- simu[i,,]
  tmp <- tmp[,!is.na(tmp[2,]) & tmp[5,]!=1]
  if(!is.null(nrow(tmp))) points(tmp[1,], tmp[3,], col = c(rainbow(sdpM@Init$NSites+1, alpha = 0.6)[tmp[2,]+1]), pch = 16) else {
    points(tmp[1], tmp[3], col = c(rainbow(sdpM@Init$NSites, alpha = 0.6)[tmp[2]]), pch = 16) 
  }
}
legend("topleft", paste("site", 1:sdpM@Init$NSites), pch = 16, col = rainbow(sdpM@Init$NSites), bty="n",ncol=2)

#########


cls <- rainbow(sdp@Init$NSites+1)
matplot(t(simu[,3,]), type = "p", pch = 16, cex = .5, col = "grey90")
matplot(t(simu[,2,]), type = "p", pch = 16, cex = .5, col = "grey90")




