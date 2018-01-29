require(plyr)
require(data.table)
#require(Cairo)
#CairoWin()
# Uncomment the previous two lines to use the Cairo package. This will prompt a new window
# where the plot, antialiased, will be generated, as RStudio for Windows does not antialise plots.

# Parameters
params <- list(
  a = 1.1, 
  b = 0.015, 
  t0 = 11,
  decay = 0.5, # I changed the sign since it gave a more approximate result, default is 0.5
  tau = 20,
  subjects = 50, 
  numBlocks = 2,
  trials = 90,
  groups = c("a", "b", "c")
)

# Check if value is even, returns 1 if TRUE, returns 2 if FALSE:
is.even <- function(number){
  if(number %% 2 == 0){
    1
  }
  else{
    2
  }
}

# Act-R noise
actr.noise <- function(sd,n=1) {
  rand <- runif(n,min=0.0001,max=0.9999)
  sd * log((1 - rand ) / rand)
}

# Miliseconds to pulses
msecToPulses <- function(t) {
  tick <- params$t0
  timeCounter <- tick
  tickCounter <- 0
  while(timeCounter <= t){
    sd <- params$b * params$a * tick
    tick <- params$a * tick + actr.noise(sd)
    timeCounter <- timeCounter + tick
    tickCounter <- tickCounter + 1
  }
  tickCounter
}

# Pulses to miliseconds
pulsesToMsec <- function(p) {
  timeCounter <- 0
  tick <- params$t0
  while(p > 0){
    sd <- params$b * params$a * tick
    tick <- params$a * tick + actr.noise(sd)  
    timeCounter <- timeCounter + tick
    p <- p - 1
  }
  timeCounter
}

# Setting up the basics of the experiment  
setupExperiment <- function(currentSubject){
  time <- 0 #Real time on the experiment
  D <- matrix(NA, 35, 900) # Declarative memory
  block <- 1
  shortValuesB1 <- c(630, 750, 900)
  longValuesB1 <- c(900, 1080, 1300)
  shortValuesB2 <- c(630, 750, 690)
  longValuesB2 <- c(1300, 1080, 1190)
  anchorValueB2 <- 900
  fixation1 <- runif(1, 1500, 2000)
  fixation2 <- runif(1, 2000, 3000)
  group <- sample(params$groups, 1, replace = TRUE)
  shortSampleB1 <- sample(shortValuesB1, params$trials, replace = TRUE)
  longSampleB1 <- sample(longValuesB1, params$trials, replace = TRUE)
  shortSampleB2 <- sample(c(shortValuesB2[which(letters == group)], anchorValueB2),
                          params$trials, replace = TRUE)
  longSampleB2 <- sample(c(longValuesB2[which(letters == group)], anchorValueB2),
                         params$trials, replace = TRUE)
  list(
    time = time,
    subjectID = currentSubject,
    block = block,
    group = group,
    D = D,
    shortSampleB1 = shortSampleB1,
    longSampleB1 = longSampleB1,
    anchorValueB2 = anchorValueB2,
    shortSampleB2 = shortSampleB2,
    longSampleB2 = longSampleB2,
    fixation1 = fixation1,
    fixation2 = fixation2
  )
  
}

# Adding encounter to DM
addEncounter <- function(subjectInfo, pulse){
  tmp <- sum(!is.na(subjectInfo$D[pulse, ])) + 1
  subjectInfo$D[pulse, tmp] <- subjectInfo$time
  subjectInfo$time <- subjectInfo$time + subjectInfo$fixation1
  subjectInfo
}

# Retrieving encounter from DM
getEncounter <- function(subjectInfo, row){
  tmp <- subjectInfo$D[row,]
  if(length(tmp[!is.na(tmp)]) == 0){
    return(NA)
  }
  else{
    return(tmp[!is.na(tmp)])
  }
}

# Chunk activation level and decay
chunkActivation <- function(encounters, currentTime){
  if(currentTime < min(encounters)){
    return(NA)
  }
  else{
    sum((currentTime - encounters[encounters<currentTime])^params$decay)
  }
}

# Probability of retrieving chunk based on its activation
calculateProbability <- function(currentActivation, totalActivation){
  prob <- (sum(currentActivation) / params$tau) / ( totalActivation / params$tau )
}

# Implementing blending
blending <- function(subjectInfo){
  blendedResponse <- 0
  totalActivation <- 0
  
  for(row in 1:nrow(subjectInfo$D)){
    numberOfPresentations <- sum(!is.na(subjectInfo$D[row, ]))
    if(numberOfPresentations > 0){
      encounters <- getEncounter(subjectInfo, row)
      currentActivation <- chunkActivation(encounters, subjectInfo$time)
      totalActivation <- totalActivation + currentActivation 
    }
  }
  for(row in 1:nrow(subjectInfo$D)){
    numberOfPresentations <- sum(!is.na(subjectInfo$D[row, ]))
    if(numberOfPresentations > 0){
      encounters <- getEncounter(subjectInfo, row)
      currentActivation <- chunkActivation(encounters, subjectInfo$time)
      prob <- calculateProbability(currentActivation, totalActivation)
      blendedResponse <- blendedResponse + (prob * row)
    }
  }
  blendedResponse
 
}

# Setting up an individual trial
experimentTrial <- function(currentCondition, currentBlock, subjectInfo){
  currentSample <- switch(currentBlock,
                          switch(currentCondition, 
                                 subjectInfo$shortSampleB1, 
                                 subjectInfo$longSampleB1),
                          switch(currentCondition, 
                                 subjectInfo$shortSampleB2, 
                                 subjectInfo$longSampleB2))
  
  for(i in 1:length(currentSample)){
    ts <- currentSample[i]
    subjectInfo$time <- subjectInfo$time + subjectInfo$fixation1
    pulse <- msecToPulses(ts)
    subjectInfo <- addEncounter(subjectInfo, pulse)
    blendedResponse <- blending(subjectInfo)
    
    # Modifying the response based on the previous one.
    
    if(i>1 && blendedResponse != (0 || msecToPulses((-500 * 
                                                     as.numeric(subjectInfo$subjectData[i - 1, 7])) /
                                                    tp ))){ 
      blendedResponse <- blendedResponse  - 
        msecToPulses((-500 * as.numeric(subjectInfo$subjectData[i - 1, 7])) / tp ) 
      }

    if(i>1 && blendedResponse != (0 || msecToPulses((-100 * 
                                                     as.numeric(subjectInfo$subjectData[i - 1, 7])) / 
                                                    tp ))){ 
      blendedResponse <- blendedResponse  -  
        msecToPulses((-100 * as.numeric(subjectInfo$subjectData[i - 1, 7])) / tp ) 
    }

    if(currentSample[i] == min(currentSample)) {
      blendedResponse <- blendedResponse * 0.98
      }
    else if(currentSample[i] == max(currentSample)) {
          blendedResponse <- blendedResponse * 1.02
    }
    else if(currentSample[i] != (max(currentSample) || min(currentSample))){
      switch(currentCondition, blendedResponse <- blendedResponse * 0.998, 
             blendedResponse <- blendedResponse * 1.002)
    }
    
    tp <- pulsesToMsec(blendedResponse)

    
    subjectInfo$time <- subjectInfo$time + tp + subjectInfo$fixation2
    subjectInfo$subjectData <- rbind(subjectInfo$subjectData, c(subjectInfo$subjectID, 
                                                                subjectInfo$group, 
                                                                currentBlock,
                                                                currentCondition, i, ts, tp)) 
  }
  
  subjectInfo
}

# Training the model with the time stimulus
train <- function(currentCondition, currentBlock, subjectInfo){
  currentSample <- switch(currentBlock,
                          switch(currentCondition, 
                                 subjectInfo$shortSampleB1, 
                                 subjectInfo$longSampleB1),
                          switch(currentCondition, 
                                 subjectInfo$shortSampleB2, 
                                 subjectInfo$longSampleB2))
  for(i in 1:params$trials){
    ts <-  currentSample[i]
    subjectInfo$time <- subjectInfo$time + ts
    pulse <- msecToPulses(ts)
    subjectInfo <- addEncounter(subjectInfo, pulse)
  }
  subjectInfo
}

# Combining it all in the framework of the experiment
experiment <- function(){
  simulatedData <- data.frame(matrix(ncol=6, nrow=0))
  
  
  # Loop per subject
  for(subjectNumber in 1:params$subjects){
    for(currentBlock in 1:params$numBlocks){
      subjectInfo <- setupExperiment(subjectNumber)
      subjectInfo$D <- matrix(NA, 35, 900)
      subjectInfo <- train(is.even(subjectNumber), currentBlock, subjectInfo)
      subjectInfo <- experimentTrial(is.even(subjectNumber), currentBlock, subjectInfo)
      simulatedData <- rbind(simulatedData, subjectInfo$subjectData)
    }
    cat("Processed subject", subjectNumber, "of", params$subjects, "\n")
  }

  colnames(simulatedData) <- c("subject ID", "Group", "Block", "Cond", "Trial", "Ts", "Tp")
  simulatedData

}

# Run experiment

run <- function(){
  experiment <- experiment()
  setDT(experiment)
  experiment[,Tp := as.numeric(as.character(Tp))]
  experiment[,Ts := as.numeric(as.character(Ts))]
  experiment[,plot(density(Tp),lwd=2,xlim=c(0,10000), ylim=c(0,0.004))]
  experiment[Cond==1,lines(density(Tp),col="red")]
  experiment[Cond==2,lines(density(Tp),col="darkgreen")]
  ssDev <- experiment[Block==1 & Tp>250 & Tp<1750,list(SsDev=mean(Tp)-mean(Ts)),by=list(`subject ID`)]
  experiment <- merge(experiment,ssDev,by="subject ID")
  abline(v=250)
  abline(v=1750)
  
  ## Create means per participant:
  plotdat <- experiment[,list(mTp=mean(Tp),nTp=length(Tp),sdTp=sd(Tp)),by=list(Ts,Group, Block, Cond, Trial, `subject ID`)]
  
  plotdat <- plotdat[,list(mTp=mean(mTp),sdTp=sd(mTp)),by=list(Cond, Ts, Block, Group)]
  
  plotdatBlock1 <- plotdat[Block==1,list(mTp=mean(mTp),sdTp=sd(mTp)),by=list(Block,Cond,Ts)]
  
  setkey(plotdatBlock1,Ts)
  plotdatBlock1[Cond==1,plot(Ts,mTp,type="b",pch=20,xlim=c(500,1500),ylim=c(500,1500), col="red")]
  plotdatBlock1[Cond==2,lines(Ts,mTp,type="b",pch=20,xlim=c(500,1500),ylim=c(500,1500),col="darkgreen")]
  
  abline(a=0,b=1,lty=2)
  
  plotdatBlock2 <- plotdat[Block==2,]
  setkey(plotdatBlock2,Ts)
  
  for (group in c("a","b","c")) {
    colourID <- 1
    for (lengthC in c(1,2)) {
      plotdatBlock2[Group==group & Cond==lengthC,
                    lines(Ts,mTp,type="b",pch=toupper(group), col=colourID)]
    }
    colourID <- colourID + 1 
  }
}
