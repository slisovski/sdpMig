## Install spdMig package (only once)
# library(devtools)
# install_github("slisovski/sdpMig", force = TRUE)


library(spdMig)

###############################################
### Parameters ################################
###############################################

parms <- list(
  MaxT   = 180,  ## Maximum time
  MaxX   = 100,  ## Maximum body condition
  NSites = 6,    ## Excl. Breeding site
  
  ### Species specific Parameters ##
  
  B0 = 3,       ## Future reproductive success
  w  = 0.028,   ## Parameters for sigmoidal TR function
  xc = 35,      ## 
  max_u = 1.0,  ## maximum foraging intensity
  f = 1.0,      ## proportion of x that cannot be devoted to flying
  
  
  ## Flying capacities
  c     = 11677,
  speed = 1440,
  
  ## Wind
  WindAssist = 0,
  WindProb   = 1,
  
  ## Decision Error
  decError = 5000,
  
  ## Terminal Reward
  xFTReward = c(0, 114, 139, 151, 180),
  yFTReward = c(0, 0,   2,  0,  0),
  
  
  ### Site specific Parameters ###
  path = "g://Dropbox//Modelling_Migration_SDP//Packages//sdpMig//Parameters//RedKnots.csv",
  
  pred_a1  = 2,
  pred_a2  = 2,
  
  ### Accuracy
  ZStdNorm = c(-2.5, -2.0, -1.5, -1.0, -0.5,  0.0,  0.5,  1.0,  1.5,  2.0,  2.5),
  PStdNorm = c(0.0092, 0.0279, 0.0655, 0.1210, 0.1747, 0.2034, 0.1747, 0.1210, 0.0655, 0.0279, 0.0092)
  
) ## End parameter ####


#######################
### Simulation ########
#######################

## Einzelrun
# 
# sdp  <- makeSDPmig(parms, "RedKnots")
# sdpM <- bwdIteration(sdp, pbar = TRUE)
# simu <- MigSim(sdpM, 100, 1, 2, c(15, 25)) #MigSim(object, NrInd, starttime, startsite, start_x(min, max))
# smP  <- simuPlot(simu, sdpM, fun = "mean")
# smP[apply(smP, 1, function(x) any(is.na(x))),] <- cbind(rep(NA, ncol(smP)))
# mean.staging <- apply(smP, 2, mean, na.rm = T); mean.staging


## Multiple runs in scenarios

setwd("g:/Dropbox/Modelling_Migration_SDP/packages/sdpMig/Output/")
Out.file <- "test_redknots.csv"

col.heads <- array(0, dim=(parms$NSites+1)*2+2)
col.heads[1] <- "Scenario"; col.heads[2] <- "mort"
for (c in 1:(parms$NSites+1)){
  col.heads[c+2] <- paste("stag",c)
  col.heads[c+(parms$NSites+1)+2] <- paste("meanX",c)
  col.heads[c+2*(parms$NSites+1)+2] <- paste("fuelrate",c)
}
write.table(rbind(col.heads), file=Out.file, sep=",", append=F, col.names = F, row.names = F )


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
  simu <- MigSim(sdpM, 100, 1, 2, c(10, 15)) #MigSim(object, NrInd, starttime, startsite, start_x(min, max))
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
  
  write(cbind(sims[scenario], mort, t(mean.staging), t(mean.bodycond), t(mean.fuelrates) ), 
        file=Out.file, sep=",", append=T, ncolumns = (sdp@Init$NSites+1) * 3 + 2)
}
 
