## Install spdMig package (only once)
# library(devtools)
# install_github("slisovski/sdpMig", force = TRUE)
 
library(sdpMig)
  library(zoo)
  library(maptools)
  data(wrld_simpl)

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
  path = "Parameters/WFGParameter_ex4.csv",

  pred_a1  = 2,
  pred_a2  = 2,
  
  ### Accuracy
  ZStdNorm = c(-2.5, -2.0, -1.5, -1.0, -0.5,  0.0,  0.5,  1.0,  1.5,  2.0,  2.5),
  PStdNorm = c(0.0092, 0.0279, 0.0655, 0.1210, 0.1747, 0.2034, 0.1747, 0.1210, 0.0655, 0.0279, 0.0092)
  
) ## End parameter ####


## Site specific paramters exploration
siteTab <- read.csv(parms$path, skip = 2, sep = ",", dec = ".")
head(siteTab)


cls <-  topo.colors(nrow(siteTab))

### Sites
opar <- par(mar = c(6,6,1,1))
plot(siteTab$Lon, siteTab$Lat, type = "n", las = 1, xlab = "Longitude", ylab = "Latitude", bty = "n", mgp = c(4,2.7,2))
plot(wrld_simpl, col = "grey90", border = "grey50", add = T)
points(siteTab$Lon, siteTab$Lat, type = "b", cex = 5, pch= 21, lwd = 2, bg = adjustcolor(cls, alpha.f = 0.7))
text(siteTab$Lon, siteTab$Lat, 1:nrow(siteTab), cex = 1.1)
par(opar)

### Quality (energy intake)
time <- siteTab[,substring(names(siteTab), 1,1)=="x"]
dei  <- siteTab[,substring(names(siteTab), 1,1)=="y"]

plot(NA, xlim = c(range(time, na.rm = T)), ylim = range(dei, na.rm = T), xlab = "time", ylab = "MEI (x)")
for(i in 1:nrow(siteTab)) lines(time[i,], dei[i,], type = "o", pch = 16, lwd = 2, col = cls[i])

par(new = T)
plot(parms$xFTReward, parms$yFTReward, xlim = c(range(time, na.rm = T)), ylim = range(parms$yFTReward, na.rm = T),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", type = "l", lty = 2)



#############################
### Backward Simulation #####
#############################


sdp  <- makeSDPmig(parms, "Whitefronted Geese")
sdpM <- bwdIteration(sdp, pbar = TRUE)


#############################
### Forward Simulation ######
#############################

simu <- MigSim(sdpM, 100, 1, 1, c(25, 35))

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
    legend("topleft", paste("site", 1:sdpM@Init$NSites), pch = 16, col = rainbow(sdpM@Init$NSites), bty="n",ncol=2)
    #########




