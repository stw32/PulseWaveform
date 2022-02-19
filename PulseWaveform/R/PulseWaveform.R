# PPG function list:

# 1. Preproc  (note for bioradio data only)
# 2. find_w   
# 3. find_u_v
# 4. find_o
# 5. preclean_wuv
# 6. Baseline
# 7. Clean_wuv
# 8. sep_beats
# 9. find_average
# 10. find_sd
# 11. diast_pk     
# 12. osnd_of_average
# 13. feature_extract


preproc <- function(data){
  dat<-data[!(data$PPG.PulseOx1=='NaN'),]
  
  #Downsample
  # The BioRadio device provides 250 samples per second, but the PPG is only sampled 75 times per second, 
  # so we typically have repeated values in a pattern of 3-3-4 repeating. DownSample tries to retrieve the 
  # unique values, with an effort to be robust against variation in the repeat pattern and also against 
  # genuine repeated values.
  
  list<-rle(dat$PPG.PulseOx1)
  ID <- rep(1:length(list$values), times = list$lengths)
  data_downsampled <- c()
  
  nSrc <- nrow(dat)
  iDst <- 1
  iSrc <- 1
  iVal <- 1
  print("Deduplicating data...")
  if(TRUE){ # Mode switch in case new method is worse or somehow broken
    startClock <- Sys.time()
    reportClock <- 10
    
    data_downsampled <- cbind(dat, ID)
    
    while (iSrc <= nSrc){
      if (as.double(difftime(Sys.time(),startClock,unit="secs"))>reportClock){
        min <- as.integer(reportClock/60)
        sec <- reportClock %% 60
        time <- paste("[",min,":",sep="")
        if (sec == 0){ time <- paste(time,"00]",sep="")} else { time <- paste(time,sec,"]",sep="")}
        print(paste(time,iSrc,"of",nSrc))
        reportClock <- reportClock + 10
      }
      
      data_downsampled[iDst,] <- data_downsampled[iSrc,]
      iDst <- iDst + 1
      if (list$lengths[iVal] <= 4){
        iSrc <- iSrc + list$lengths[iVal]
        iVal <- iVal + 1
      } else {
        iSrc <- iSrc + 4
        list$lengths[iVal] <- list$lengths[iVal] - 4
      }
    }
    
    # Trim downsampled data to size
    data_downsampled <- data_downsampled[1:(iDst-1),]
    
  }else{ # Old method, using rbind which gets slow
    
    data2 <- cbind(dat, ID)
    data_downsampled <-c()
    
    for (i in 1:max(ID)){
      sub.data <- dplyr::filter(data2, ID == i)
      if(nrow(sub.data) <= 4){
        data_downsampled <- rbind(data_downsampled, sub.data[1,])
      }else if(nrow(sub.data) > 4 ){data_downsampled <- rbind(data_downsampled, sub.data[1,], sub.data[5,])}
    }
    
  } # End of method selector
  
  print("Removing DC blocker...")
  #Undetrend
  # Analysis of device output indicates that the PPG signal is detrended by application of the following
  # formula: OUT[i] = 80 + (OUT[i-1]-80) * 0.96875 + (IN[i] - [IN[i-1]), where the constant 0.96875 is
  # an approximation fitted to the data.
  # Individual pulse events are more comprehensible if the detrending is not used, so this function 
  #removes it by inverting the above function. 
  undetrended <-replicate(length(data_downsampled$PPG.PulseOx1)-1, 0) 
  undetrended<-c(data_downsampled$PPG.PulseOx1[1],undetrended) #add first detrended value to vector
  for (i in 2:length(data_downsampled$PPG.PulseOx1)){
    undetrended[i]<-((data_downsampled$PPG.PulseOx1[i]-80) - ((data_downsampled$PPG.PulseOx1[i-1]-80) * 0.96875) + (undetrended[i-1]))
  }
  print("Done")
  return(undetrended)
}



find_w <- function(d1p, deriv1, sp, sr){
  ########################################################################################################################################
  # FindW identifies peaks in the first derivative of the PPG time series (denoted w on the original PPG pulse wave (Elgendi et al, 2018)).  
  # A rolling window relative to heart rate is applied to identify beats and artefacts. 
  # Peaks identified within a window are confirmed as genuine / artefactual with a series of checks against thresholds derived from beats 
  # local to the window and across the entire time series. 
  
  # Inputs: 
  # d1p (polynomial spline of first derivative of the PPG time series)
  # deriv1 (first derivative in discrete form)
  # sp (polynomial spline of original time series)
  # sr (sample rate)
  
  # Outputs:
  # wX (x axis coordinates of all w points) 
  # wY (y axis coordinates of all w points on original PPG trace)
  # wYD1 (y axis coordinates of all w points on 1st derivative trace)
  ########################################################################################################################################
  
  d1InflxX <- solve(d1p, b = 0, deriv = 1)                 # Find inflection points on 1st deriv
  d1InflxY <- predict(d1p, d1InflxX)
  wX <- c()                                                # Create vectors for w values to be stored
  wY <- c()
  window <- list()                                         # Define window within which to identify peaks 
  prevPkDist <- c()                                        # Create vector to store peak to peak distances
  
  ############# Identify first two peaks: #############
  
  a <- 2                                                                                                                                        # Start with a very small window (spanning three (1 + a) inflection points), 
  while(length(wX) < 2){                                                                                                                        # and widen until two peaks are identified within the window
    
    firstWindow <- data.frame(d1InflxX[1]:d1InflxX[1+a], deriv1[d1InflxX[1]:d1InflxX[1+a]])                                                     # Create window
    windowInflxY <- d1InflxY[which(d1InflxX > firstWindow[1, 1] & d1InflxX < firstWindow[, 1][length(firstWindow[, 1])])]                       # Find the inflection points which fall within the window 
    windowInflxX <- d1InflxX[which(d1InflxX > firstWindow[1, 1]  & d1InflxX < firstWindow[, 1][length(firstWindow[, 1])])]                     
    threshold <- quantile(firstWindow[, 2], probs=c(.95))                                                                                       # Define a threshold for finding peaks, of 0.95 
    windowPks <- which(windowInflxY > threshold)                                                                                                # Identify inflection points above the threshold as peaks
    
    if(length(windowPks) == 2 & windowInflxX[windowPks[2]] - windowInflxX[windowPks[1]] > sr){                                                                                      # If the first two peaks indicate a heart rate < 60 bpm, 
      confirmWindow <- data.frame(d1InflxX[windowPks[1]+2]:d1InflxX[windowPks[2]-2],   deriv1[d1InflxX[windowPks[1]+2]:d1InflxX[windowPks[2]-2]])       # ensure a peak between them has not been missed with the current threshold
      inflxConformY <- d1InflxY[which(d1InflxX > confirmWindow[1, 1] & d1InflxX < confirmWindow[, 1][length(confirmWindow[, 1])])]
      inflxConformX <- d1InflxX[which(d1InflxX > confirmWindow[1, 1]  & d1InflxX < confirmWindow[, 1][length(confirmWindow[, 1])])]                     # Find inflection points within the window (defined as between peaks 1 and 2)
      threshold <- max(quantile(confirmWindow[, 2], probs=c(.95)), windowInflxY[windowPks[1]]/2)     
      missedPks <- inflxConformX[which(inflxConformY > threshold)]                                                                                      # Identify peaks above threshold (and greater than 50% of the height of the first peak)
      missedPks <- which(windowInflxX == missedPks[1])                                                                                                  # Identify the inflection point in the original window which corresponds to the missed peak 
      if(length(missedPks) > 0){windowPks[2] <- missedPks[1]}                                                                                                                       # Assign the missed peak as the second peak
    }
    
    if(length(windowPks) == 2 & windowInflxX[windowPks[2]] - windowInflxX[windowPks[1]] < (sr/10)){                                                                         # Rarely, the first peak contains two inflection points above threshold, 
      windowPks <- windowPks[-2]                                                                                                                # if peaks are too close (< 1/10th of a second), the second peak is discounted
    }
    
    windowPksY <- windowInflxY[windowPks]                                                                                                       # Identify the inflection points corresponding to peaks
    windowPks <- windowInflxX[windowPks]     
    
    if(length(windowPks) > 2 & (windowPks[3] - windowPks[1] < 1)){                                                                              # Very rarely, the first peak contains three inflection points above threshold. Discount two. 
      windowPks <- windowPks[-c(1, 3)]
    }
    
    m <- mean(d1InflxY[1:(1+a)])                                                                                                                # Calculate mean and SD of inflection points within the window
    std <- sd(d1InflxY[1:(1+a)])
    
    if(length(windowPks) == 2){                                                                                                                 # If two peaks are found, confirm that they are both significantly higher than surrounding                                                                                 
      if(windowPksY[1] > (m + (1.5*std)) & windowPksY[1] > (median(deriv1[1:100])+std(deriv1[1:100]))){                                         # inflection points, and that the second peak is greater than half the height of the first     
        wX[1] <-  windowPks[1]                                                                                                                  # peak (unless the first peak is an artefact)
        if(windowPksY[2] > (m + (1.5*std)) & windowPksY[2] > (windowPksY[1]/2) | windowPksY[1] > (mean(deriv1) + (5*std(deriv1)))){  
          wX[2] <-  windowPks[2]                                                                                                                          
        }
      }
    }   
    
    if(length(windowPks) > 3){                                                                                                                  # If more than 3 peaks identified, assume the time series begins with an artefact and skip 
      d1InflxX <- d1InflxX[-c(1:(sr*(100/75)))]                                                                                                 # forward, resetting a. 
      d1InflxY <- d1InflxY[-c(1:(sr*(100/75)))]
      a <- -3  
      print("skipped")
    }
    
    a <- a + 5                                                                                                                                  # Increases until two legitimate peaks are found
  }
  
  ############# Find remaining peaks (with dynamic moving window dependent on peak to peak distance): #############
  
  artefacts <- c()                                                                                                                              # Create vector to store artefactual peaks
  m <- mean(d1InflxY)                                                                                                                           # Calculate mean and SD of inflection points across the time series
  std <- sd(d1InflxY)
  for(i in 3:length(d1InflxX)){  
    prevPkDist[i] <- wX[length(wX)] - wX[length(wX)-1]                                                                                          # Peak to peak distance defined as distance between two previous peaks
    if(prevPkDist[i] > (mean(prevPkDist[!is.na(prevPkDist)])*2)){                                                                               # Identify abnormally large peak to peak distances (this can be because 
      prevPkDist[i] <- prevPkDist[3]                                                                                                            # there was a large gap between an artefact and the next legitimate peak) 
    }                                                                                                                                           # and correct them by resetting the value (this avoids leapfrogging)     
    if((wX[length(wX)] + (4*prevPkDist[i]) > length(deriv1))){                            
      break                                                                                                                                     # If the next window goes beyond the length of the data, break the loop
    }                                      
    windowPks <- c()
    windowExtnd <- 1.35                                                                                                                         # Define parameters by which to extend the window or move it forward
    windowStart <- 0.5
    printed <- NA
    if(length(artefacts) > 0){                                                                                                                  # Recalculate mean and standard deviation after removing any artefacts identified 
      remove <- c()
      for(j in artefacts[which(artefacts < length(wX))]){
        remove[j] <- which(abs(d1InflxX - wX[j]) == min(abs(d1InflxX - wX[j])))
      }
      remove <- remove[!is.na(remove)]
      newRem <- c()
      for(j in 1:length(remove)){
        newRem[(length(newRem)+1):(length(newRem)+21)] <- (remove[j] -10): (remove[j] + 10)                                                     # Remove inflection points around artefact peaks
      }
      if(sum(newRem < 1) > 0){  
        newRem <- newRem[-(which(newRem < 1))] 
      }    
      m <- mean(d1InflxY[-newRem])       
      std <- sd(d1InflxY[-newRem])  
    }
    
    while(length(windowPks) < 1){                                                                                                               # Each subsequent window will adjust (if required) until the next peak is detected
      
      if(windowExtnd > 10){                  
        windowStart <- 2
        windowExtnd <- 2.5
      }
      
      window[[i]] <- data.frame((wX[length(wX)] + windowStart*prevPkDist[i]):(wX[length(wX)] + windowExtnd*prevPkDist[i]),                      # Define the window initially as from (the previous peak + (peak to peak distance / 2)) 
                                deriv1[(wX[length(wX)] + windowStart*prevPkDist[i]):(wX[length(wX)] + windowExtnd*prevPkDist[i])])                         # to (previous peak + (peak to peak distance * 1.35))
      windowInflxY <- d1InflxY[which(d1InflxX > window[[i]][1, 1] & d1InflxX < window[[i]][length(window[[i]][, 1]), 1])]        
      windowInflxX <- d1InflxX[which(d1InflxX > window[[i]][1, 1] & d1InflxX < window[[i]][length(window[[i]][, 1]), 1])]                       # Identify which inflection points are within the window       
      threshold <- quantile(window[[i]][, 2], probs=c(0.95))                                                                                    # Define threshold                                 
      windowPks <- which(windowInflxY > threshold)                                                                                              # Identify peaks
      windowPksY <- windowInflxY[windowPks]
      windowPks <- windowInflxX[windowPks]     
      
      # plot(window[[i]], type = "l")                                                                                                           # For debugging purposes, the window can be plotted
      # points(windowPks, windowPksY, pch = 19)
      
      ############# Window Peak Confirmation Criteria: #############
      
      if(length(windowPks) > 2){                                                                                                                # If more than 2 peaks are identified within the window, this can be due to four reasons:
        if(max(windowPksY) > (m+(2*std)) |  max(windowPksY) > (mean(deriv1[windowPks[1]:(windowPks[1]                                           # 1. a window does not include a genuine peak, and there are multiple secondary ones of similar height 
                                                                                         + (3*prevPkDist[i]))]) + 2*std(deriv1[windowPks[1]:(windowPks[1] + (3*prevPkDist[i]))]))){                           # 2. there is an artefact with more than two peaks
          wX[length(wX)+1] <- windowPks[which(windowPksY == max(windowPksY))]                                          # 3. there are significantly large secondary peaks that also exceed the threshold for identification
          lowPks <- windowPksY[order(windowPksY)[1:2]]                                                                 # 4. a genuine peak has multiple inflection points   
          if(lowPks[1] > m+(2*std) & lowPks[2] > m+(2*std)                                                                                      
             & (windowPks[3] - windowPks[1]) > (prevPkDist[i]/10)){                                                                             # Check if the maximum peak is above a threshold relative to the time series (no 1)    
            cat('\n','Potential artefact',  ', plot(', (wX[i-1]-100), ':', (wX[i-1]+300), ', deriv1[', (wX[i-1]-100), ':', (wX[i-1]+300),    # If it is, check if both lower peaks also exceed the threshold, 
                '], type = "l") ,', 'wave', i, '+/- 2 removed because the non-max peaks were high')                                              # if so mark them as artefactual (no 2), if not assume they are secondary peaks (no 3)
            artefacts[length(artefacts) + c(1, 2, 3, 4, 5)] <- c(i-2, i-1, i, i+1, i+2)                                                      # An additional condition of the above is that the 'peaks' are not too close together so as to be inflection points (no 4)
          }
        }else{                                                                                                                                  
          windowExtnd <- windowExtnd + 0.5                                                                                                      # If no 1 is the case, extend the window to look for peaks again
          windowPks <- c()
        }
      }
      
      if(length(windowPks) == 2){                                                                                                               # If two peaks are identified, they should be confirmed as genuine:
        if(max(windowPksY) > (m+(2*std)) | max(windowPksY) > (mean(deriv1[windowPks[1]:(windowPks[1]                                            # Check if the maximum peak exceeds the global or local threshold 
                                                                                        + (3*prevPkDist[i]))]) + 2*std(deriv1[windowPks[1]:(windowPks[1] + (3*prevPkDist[i]))]))){                                           # (local defined relative to peak to peak distance)
          if((windowPks[2] - windowPks[1]) < (prevPkDist[i]/3)){                                                                                # If it does, check if the two peaks are close togetber in time 
            wX[length(wX)+1] <- windowPks[which(windowPksY == max(windowPksY))]                                                                 # (within 1/3rd of the peak to peak distance) 
            cat('\n','Potential artefact',  ', plot(', (wX[i-1]-100), ':', (wX[i-1]+300),                                                       # If they are, mark the highest peak as artefactual
                ', deriv1[', (wX[i-1]-100), ':', (wX[i-1]+300), '], type = "l") ,', 'wave',
                i, '+/- 2 removed because two peaks found too close together') 
          }else{                                                                                                                                # If they are not, go through both peaks in turn to confirm if genuine:
            if((windowPksY[1] > (m+(2*std))) | windowPksY[1] > (mean(deriv1[windowPks[1]:(windowPks[1]                                          # Check if the first peak can be marked as genuine by:
                                                                                          + (3*prevPkDist[i]))]) + 2*std(deriv1[windowPks[1]:(windowPks[1] + (3*prevPkDist[i]))])) &                                        # 1. comparing to time series threshold 
               windowPksY[1] > (windowPksY[2]/2)){                                                                                               # 2. comparing to local threshold 
              wX[length(wX)+1] <- windowPks[1]                                                                                                # 3. comparing to height of the other peak (should exceed it's half maximum)
              if(windowPksY[2] > (m+(2*std)) | windowPksY[2] > (mean(deriv1[windowPks[1]:(windowPks[1]                                        # If it can, the second peak can be marked as genuine with the same criteria
                                                                                          + (3*prevPkDist[i]))]) + 2*std(deriv1[windowPks[1]:(windowPks[1] + (3*prevPkDist[i]))])) &                                      
                 windowPksY[2] > (windowPksY[1]/2)){   
                wX[length(wX)+1] <- windowPks[2]  
              }
            }else{                                                                                                                              # If the first peak is not genuine, check the second and identify only the 
              if(windowPksY[2] > (m+(2*std)) | windowPksY[2] > (mean(deriv1[windowPks[1]:(windowPks[1] +                                        # second peak as genuine if appropriate. 
                                                                                          (3*prevPkDist[i]))]) + 2*std(deriv1[windowPks[1]:(windowPks[1] + (3*prevPkDist[i]))])) & 
                 windowPksY[2] > (windowPksY[1]/2)){
                wX[length(wX)+1] <- windowPks[2] 
              }
            }
          }
        }else{                                                                                                                                  # If the maximum peak did not exceed threshold, extend the window and look again. 
          windowExtnd <- windowExtnd + 0.5
          windowPks <- c()
        }   
      }
      
      if(length(windowPks) == 1){                                                                                                               # If one peak is identified, confirm it is genuine by: 
        if(windowPksY > (m+(2*std)) | windowPksY > (mean(deriv1[windowPks[1]:(windowPks[1] +                                                    # 1. comparing to time series threshold
                                                                              (3*prevPkDist[i]))]) + 2*std(deriv1[windowPks[1]:(windowPks[1] + (3*prevPkDist[i]))])) |                                              # 2. comparing to local threshold
           windowPksY > (predict(d1p, wX[i-1])*0.9)){                                                                                            # 3. comparing it to the height of the previous peak                               
          wX[length(wX)+1] <- windowPks 
        }else{
          if((i-1) %in% artefacts){                                                                                                             # If the above criteria are not met, check if the previous peak was artefactual
            wX[length(wX)+1] <- windowPks                                                                                           
          }else{                                                                                                                                # If it was not, extend the window
            windowExtnd <- windowExtnd + 0.5
            windowPks <- c()
          }
        }
      }
      
      #############Window Artefact Identification: #############
      
      if(windowExtnd > 2 & is.na(printed)){                                                                                                     # If the window has been extended more than twice to find a peak, 
        cat('\n','Potential artefact',  ', plot(', (wX[i-1]-100), ':', (wX[i-1]+300), ', deriv1[',                                              # it is likely that that peak may be an artefact; label it as such
            (wX[i-1]-100), ':', (wX[i-1]+300), '], type = "l") ,', 'wave', i, '+/- 2 removed because window needed extending')  
        printed <- 1 
        artefacts[length(artefacts) + c(1, 2, 3, 4, 5)] <- c(i-2, i-1, i, i+1, i+2)                                                             # i-1 and subsequent adjustments are used because wX[i] does not yet exist
      }
      
      if(min(window[[i]]) < (m-(3*std)) & is.na(printed)){                                                                                      # If a window contains a value that drops considerably below the mean, 
        cat('\n','Potential artefact',  ', plot(', (wX[i-1]-100), ':', (wX[i-1]+300), ', deriv1[',                                              # label the peak found in that window as an artefact 
            (wX[i-1]-100), ':', (wX[i-1]+300), '], type = "l") ,', 'wave', i, '+/- 2 removed because the window contains a very low value') 
        artefacts[length(artefacts) + c(1, 2, 3, 4, 5)] <- c(i-2, i-1, i, i+1, i+2)
      }
      
      windowExtnd <- windowExtnd + 0.1                                                                                                          # If no peaks are found in a window (or only spurious ones found), 
    }                                                                                                                                           # extend the window and look again
  }       
  
  ############# Post Peak Identification Cleaning: #############
  
  d1wY <- predict(d1p, wX)                                                                                                                      # Find summary statistics of identified peak y-values
  m <- mean(d1wY) 
  std <- sd(d1wY) 
  if(length(artefacts) > 0){                                                                                                                    # Remove identified artefacts (should this be before calculating mean??)
    wX <- wX[-artefacts]
    d1wY <- d1wY[-artefacts]
  }
  
  artefacts <- c()                                                                                                                              # Check for unusually tall peaks, and label them as artefacts if:
  for(i in 2:(length(wX)-1)){     
    if(d1wY[i] > 1.5*d1wY[i+1] & d1wY[i] > 1.5*d1wY[i-1] | d1wY[i] > (m+(5*std))){                                                              # 1. the ith wave height exceeds local (comparison with i-1 and i+1) or time series thresholds
      cat('\n','Potential artefact',  ', plot(', (wX[i-1]-100), ':', (wX[i-1]+300), ', deriv1[',
          (wX[i-1]-100), ':', (wX[i-1]+300), '], type = "l") ,', 'wave', i, '+/- 2 removed because of a very tall peak')
      artefacts[length(artefacts) + c(1, 2, 3, 4, 5)] <- c(i-2, i-1, i, i+1, i+2)
    }
    if(i > 5 & i < (length(wX)-5)){
      if(d1wY[i] > (1.5*mean(d1wY[c((i-5):(i+5))]))){
        cat('\n','Potential artefact',  ', plot(', (wX[i-1]-100), ':', (wX[i-1]+300), ', deriv1[',                                              # 2. the ith wave height exceeds another local threshold (derived from i-5 to i+5)  
            (wX[i-1]-100), ':', (wX[i-1]+300), '], type = "l") ,', 'wave', i, '+/- 2 removed because of a very tall peak')                        
        artefacts[length(artefacts) + c(1, 2, 3, 4, 5)] <- c(i-2, i-1, i, i+1, i+2)
      }
    }
  }       
  
  ############# Output: #############
  
  wY <- predict(sp, wX)                                    # Find W y_values on original trace
  Ws <- data.frame(wX, wY, d1wY)                           # Create output dataframe 
  if(length(artefacts) > 0){                               # Remove any further artefacts identified
    Ws <- Ws[-artefacts, ]
  }
  Ws <- Ws[-1, ]                                           # Remove the first and last waves
  Ws <-Ws[-nrow(Ws), ]
  
  colnames(Ws) <- c("wX", "wY", "wYD1")                    # Outputs as described above
  return(Ws)
}


find_u_v <- function(wx, wy, d1, d1p, spline, sr = samplingRate, plot = FALSE){
  ########################################################################################################################################
  # FindUV identifies the points on the systolic upstroke that correspond to the two half maximum values on the corresponding first  
  # derivative peak. The first, before w, is denoted U. The second, after w, is denoted V. 
  
  # Inputs: 
  # wx (A vector of x-coordinates corresponding to w points on the PPG time series)
  # wy (A vector of y-coordinates corresponding to w points on the PPG time series)
  # d1 (The first derivative time series, discrete form)
  # d1p (The first derivative time series, polynomial spline form)
  # spline (The original PPG time series, polynomial spline form)
  # sr (samplingRate)
  # plot (Logical, will plot the first derivative time series with U and V values as points if set to true)
  
  # Outputs:
  # uX (A vector of x-xoordinates corresponding to U points on the PPG time series)
  # uY (A vector of y-xoordinates corresponding to U points      "       "        )
  # vX (A vector of x-xoordinates corresponding to V points      "       "        )
  # vY (A vector of y-xoordinates corresponding to V points      "       "        )
  ########################################################################################################################################
  
  wHalfHeight <- predict(d1p, wx)/2                                                           # For each w point, identify the half maximum on first derivative
  halfHeightX <- c()                                                                          
  halfHeightY <- c()
  for(i in 1:length(wHalfHeight)){                                                            # For each 1st derivative peak, create a polynomial spline of the peak only
    d1PeakSub <- CubicInterpSplineAsPiecePoly((round(wx[i])-(sr/8)):(round(wx[i])+(sr/8)),
                                              d1[(round(wx[i])-(sr/8)):(round(wx[i])+(sr/8))], "natural")                            
    preHalfHeights <- solve(d1PeakSub, b = wHalfHeight[i])                                    # Identify the x-coordinates of the new spline when the y-value = the half maximum
    
    if(length(preHalfHeights) < 2){                                                           # If only one coordinate is found, extend the length of the spline, 
      d1PeakSub <- CubicInterpSplineAsPiecePoly((round(wx[i])-(sr/4)):(round(wx[i])+(sr/4)),          # and search again
                                                d1[(round(wx[i])-(sr/4)):(round(wx[i])+(sr/4))], "natural")                        
      preHalfHeights <- solve(d1PeakSub, b = wHalfHeight[i]) 
    }
    
    if(length(preHalfHeights) > 2){                                                           # More than one x-coordinate is sometimes found if the gradient of the upstroke 
      a <- preHalfHeights - wx[i]                                                             # is particularly variable
      b <- c()                                                                                # In such cases, find the time (x-axis) difference between each detected half height and w
      for(j in 1:(length(a)-1)){                                                              # Select the two values with the most similar difference before and after w
        b[j] <- a[j] - a[j+1]
      } 
      b[length(b) + 1] <- a[1] - a[length(a)]
      c <- which(abs(b) == min(abs(b)))
      if(c == length(b)){
        preHalfHeights <- c(preHalfHeights[c], preHalfHeights[1])
      }else{
        preHalfHeights <- c(preHalfHeights[c], preHalfHeights[c+1])
      }
    }
    halfHeightX[c((2*(i)-1), (2*(i)))] <- preHalfHeights                                      # Add the two values per peak to a vector of all in the time series
    halfHeightY[c((2*(i)-1), (2*(i)))] <- predict(d1PeakSub,                                  # And identify the corresponding y-coordinates for the (first derivative) times series 
                                                  halfHeightX[c((2*(i)-1), (2*(i)))])  
  }
  
  if(plot){                                                                                   # Check results
    plot(d1p)
    points(halfHeightX, halfHeightY, pch = 19)
  }
  
  uX <- halfHeightX[seq_along(halfHeightX) %%2 != 0]                                          # Split the time series values into vectors for U and V
  vX <- halfHeightX[seq_along(halfHeightX) %%2 == 0]
  
  uY <- c()                                                                                   # Find the y-coordinates for U and V on the original PPG time series
  for(i in 1:length(uX)){
    uY[i] <- predict(spline, uX[i])
  }
  vY <- c()                                                                 
  for(i in 1:length(vX)){
    vY[i] <- predict(spline, vX[i])
  }
  
  df <- data.frame(uX, uY, vX, vY)                                                            # Organize outputs (see summary box above)
  return(df)
}


find_o <- function(wx, inx, iny, d1p, sp){
  ########################################################################################################################################
  # FindO identifies the origin of the systolic peak for each beat. In the absence of a clear inflection point at the origin, O points are
  # derived from inflection points in the first derivative. 
  
  # Inputs: 
  # wx (A vector of x-coordinates corresponding to w points on the PPG time series)
  # inx (A vector of all inflection point x-coordinates in the PPG time series)
  # iny (A vector of all inflection point y-coordinates in the PPG time series)
  # d1p (The first derivative time series, polynomial spline form)
  # sp (The original PPG time series, polynomial spline form)
  
  # Outputs:
  # o (vector of O points)
  # inx (vector of updated inflection point x-coordinates)
  # iny (vector of updated inflection point y-coordinates)
  ########################################################################################################################################
  
  o <- c()                                                                                    # Create vector to store O values
  for(i in 1:length(wx)){
    o[i] <- max(which(inx < wx[i]))                                                           # Identify O points as the those inflection points that 
  }                                                                                           # immediately precede w points. 
  
  inflexD1 <- solve(d1p, b = 0, deriv = 1)                                                    # Find inflection points on the first derivative
  inflexD1y <- predict(d1p, inflexD1)                                                         
  o2 <- c()                                                                                   # For each beat, find the inflection point before w that is 
  for(i in 1:length(wx)){                                                                     # also below 0 (points above 0 tend to be spurious)
    o2[i] <- max(which(inflexD1 < wx[i] & inflexD1y < 0))
  }
  
  for(i in 1:length(wx)){                                                                     # For each beat, replace the originally identified O with 
    if((inx[o][i] - inflexD1[o2][i]) < 0){                                                    # the O derived from the first derivative if the latter 
      inx[o][i] <- inflexD1[o2][i]                                                            # is closer to w in time (this occurs when there is no 
      iny[o][i] <- predict(sp, inx[o][i])                                                     # inflection point at the origin of the original PPG signal)
    }
  }
  
  tmp <- list(inx, iny, o)                                                                    # return the o points, as well as the inflection points 
  return(tmp)                                                                                 # which have been repositioned due to the above replacement
}


preclean_wuv <- function(w, uv, o, samp, sp, q = FALSE){
  ########################################################################################################################################
  # PreClean takes u, v and w values and assesses the x-axis (time) difference between them for each beat. For normal beats, time from
  # u to w is around half (50%) of the time from u to v. Beats with abnormal / artefactual systolic upstrokes tend to have outlying 
  # values for this measure. Thus PreClean can identify and remove them. 
  
  # Inputs: 
  # w (vector of all x-axis coordinates corresponding to w points)
  # uv (a dataframe containing x and y coordinates for all u and v values)
  # o (vector of O points)
  # samp (sample rate)
  # sp (The original PPG time series, discrete form)
  # q (Logical, will pause function and give the option of plotting rejected beats)
  
  # Outputs:
  # w (as above, with aberrant beats removed)
  # uv (              "       "             )
  # o (               "       "             )
  ########################################################################################################################################
  
  vDist <- c()                                                                                # For each wave, find the timing of w relative  
  for(i in 1:length(w$wX)){                                                                   # to u and v, expressed as a percentage 
    vDist[i] <- (w$wX[i] - uv$uX[i]) / (uv$vX[i] - uv$uX[i])*100
  }
  # plot(w$wX, vDist, type = "l")                                                             # This value may be of interest in itself (see plot)
  
  sdpdtv <- sd(vDist)                                                                         # Make a vector of all beats where this value is abnormal, 
  pdtvWaves <- c(which(vDist > (sdpdtv + median(vDist)) & vDist > 70),                        # either by being greater than one standard deviation 
                 which(vDist < (median(vDist) - sdpdtv) & vDist < 30))                        # from the median, or having an absolute value < 30 or > 70
  
  if(length(pdtvWaves) > 0){
    cat("\n", length(pdtvWaves), 
        "waves removed due to abnormal distances between u, v and w" )                        
    if(q == TRUE){                                                                            # If q = T, enter debug / plot routine
      plotyyy <- 0
      while(plotyyy == 0){
        plotyy <- readline(prompt = "Would you like to view? (enter yes or no)")
        if(plotyy == "yes"){
          for(i in 1:length(pdtvWaves)){
            if(pdtvWaves[i] == 1){
              plot(1:(samp*10), sp[1:(samp*10)], type = "l")
              points(w$wX, w$wY)
              points(w$wX[pdtvWaves[i]], w$wY[pdtvWaves[i]], pch =19, col = 'red')
              points(uv$uX, uv$uY)
              points(uv$vX, uv$vY)
            }else{
              plot((w$wX[pdtvWaves[i]]-samp*5):(w$wX[pdtvWaves[i]]+samp*5), 
                   sp[(w$wX[pdtvWaves[i]]-samp*5):(w$wX[pdtvWaves[i]]+samp*5)], type = "l")
              points(w$wX, w$wY)
              points(w$wX[pdtvWaves[i]], w$wY[pdtvWaves[i]], pch =19, col = 'red')
              points(uv$uX, uv$uY)
              points(uv$vX, uv$vY)
            }
          }
          plotyyy <- 1
        }
        if(plotyy == "no"){
          cat("\n", "ok") 
          plotyyy <- 1
        }
        if(plotyy != "yes" & plotyy != "no"){cat("\n", "please enter 'yes' or 'no'")}
      }
    }
    w <- w[-pdtvWaves, ]                                                                      # Remove beats with abnormal values from all relevant vectors
    uv <- uv[-pdtvWaves, ]
    o <- o[-pdtvWaves]       
  }
  return(list(w, uv, o))
}


baseline <- function(inx, iny, o, dat, sp, plot = FALSE){
  ########################################################################################################################################
  # RemoveBaselineDrift generates a spline to fit the wandering baseline of the PPG time series and subtracts it from the time series. 
  # Given baseline drift is indicative of vasomotion and respiratory modulation, it may be of interest to extract the spline.
  
  # Inputs: 
  # inx (A vector of all inflection point x-coordinates in the PPG time series)
  # iny (A vector of all inflection point y-coordinates in the PPG time series)
  # o (A vector of all O points)
  # dat (The original PPG time series (preprocessed))
  # sp (The original PPG time series, polynomial spline form)
  # plot (Logical, plots baseline spline and the original PPG time series)
  
  # Outputs:
  # baseCor (time series with baseline wander removed)
  ########################################################################################################################################
  
  sfunction2 <- splinefun(inx[o], iny[o], method = "natural")                                 # Create a spline from the O points through
  splineBase <- sfunction2(seq(1, length(dat)), deriv = 0)                                    # the time series
  
  if(plot){                                                                                   # Plot the baseline spline and time series
    plot(sp)
    points(inx[o], iny[o], pch = 19)
    lines(splineBase)
  }
  
  baseCor <- dat - splineBase                                                                 # Subtract the baseline spline from the time 
  if(plot){                                                                                   # series, thereby removing baseline drift
    plot(baseCor, type = "l")                                                                  
    lines(1:length(dat), seq(from = 0, to = 0, length.out = length(dat)))                     # Plot new baseline (y = 0)     
  }
  
  return(baseCor)
}



clean_wuv <- function(wuv, sp, inx, o, samp, bc, q = FALSE){
  ########################################################################################################################################
  # CleanWuv first identifies erroneous beats and removes them, by identifying abnormal values of U and V, and the y-axis difference 
  # between them. It then calculates O-O intervals and inter-beat intervals, and removes beats preceding large intervals. 
  
  # Inputs: 
  # wuv (a dataframe of all w, u and v values corresponding to detected and pre-cleaned beats in the time series)
  # sp (the baseline-corrected PPG time series, polynomial spline form)
  # inx (a vector of all inflection point x-coordinates in the PPG time series)
  # o (a vector of all O points)
  # samp (sampling rate)
  # bc (the baseline-corrected PPG time series, discrete form)
  # q (logical, will pause function and give the option of plotting rejected beats)
  
  # Outputs:
  # d (a dataframe combining the following three structures):
  # wuv (the inputted dataframe of w, u and v values, with rows corresponding to erroneous beats removed)
  # diffVU (a vector of y-axis differences between u and v for each beat - for scaling purposes)
  # o2 (a vector of O points, with those corresponding to erroneous beats removed)
  # ibi (a vector of all inter-beat intervals (x-axis differences between sucessive W points))
  # oDiff (a vector of all O-O intervals (x-axis differences between successive O points))
  ########################################################################################################################################
  
  o2 <- o                                                                                     # Define a second vector for O from which 
  # rejected waves are to be removed
  
  falseU <- c()                                                                               # Identify U values that are implausibly far 
  for(i in 1:(nrow(wuv)-1)){                                                                  # from baseline    
    if(wuv$uY[i] > (median(wuv$uY) + 2*std(wuv$uY)) & wuv$uY[i] > 1 &                         
       wuv$uY[i] > wuv$uY[i+1]*2 | wuv$uY[i] < -50){
      falseU[i] <- i
    }
  }
  if(length(falseU) > 0){
    falseU <- falseU[!is.na(falseU)]
    cat("\n", length(falseU), "/", nrow(wuv), 
        "waves removed due to U having an abnormally high y-value relative to baseline")
    if(q == TRUE){                                                                            # Optional debug mode for plotting waves  
      plotyyy <- 0                                                                            # with erroneous U values
      while(plotyyy == 0){
        plotyy <- readline(prompt = "Would you like to view? (enter yes or no)")
        if(plotyy == "yes"){
          for(i in 1:length(falseU)){
            plot((wuv$wX[falseU[i]]-samp*2):(wuv$wX[falseU[i]]+samp*2), 
                 bc[(wuv$wX[falseU[i]]-samp*2):(wuv$wX[falseU[i]]+samp*2)], type = "l")
            points(wuv$uX[falseU[i]], wuv$uY[falseU[i]], pch = 19)
            points(wuv$uX[falseU[i]-1], wuv$uY[falseU[i]-1])
            points(wuv$uX[falseU[i]+1], wuv$uY[falseU[i]+1])
          }
          plotyyy <- 1
        }
        if(plotyy == "no"){
          cat("\n", "ok") 
          plotyyy <- 1
        }
        if(plotyy != "yes" & plotyy != "no"){cat("\n", "please enter 'yes' or 'no'")}
      }
    }
    if(length(falseU) > 0){
      o2 <- o2[-falseU]                                                                       # Remove erroneous U values
      wuv <- wuv[-falseU, ]
    }
  }
  
  diffVU <- wuv$vY - wuv$uY                                                                   # Find y-axis v-u difference (scale factor) for each wave 
  if(length(which(diffVU == 0)) > 0){                                                         # Remove waves where the above difference is 
    dup <- which(diffVU == 0)                                                                 # 0 (indicating an artefact)
    wuv <- wuv[-dup, ]
    diffVU <- diffVU[-dup]
    cat(length(which(diffVU == 0)), "/", nrow(wuv), 
        "waves removed due to scale factors of 0 (u and v incorrectly identified)")
  }
  
  falseScale <- which(diffVU < (median(diffVU) - 5*(std(diffVU))))                            # Identify beats where the v-u difference is abnormally small  
  if(length(falseScale) > 0){                                                                 # and remove them (with optional debug mode as above)
    cat("/n", length(falseScale), "/", nrow(wuv), 
        "waves removed due to scale factors of 0 (u and v incorrectly identified)")
    if(q == TRUE){
      plotyyy <- 0
      while(plotyyy == 0){
        plotyy <- readline(prompt = "Would you like to view? (enter yes or no)")
        if(plotyy == "yes"){
          for(i in 1:length(falseScale)){
            plot((wuv$wX[falseScale[i]]-samp*2):(wuv$wX[falseScale[i]]+samp*2), 
                 bc[(wuv$wX[falseScale[i]]-samp*2):(wuv$wX[falseScale[i]]+samp*2)], type = "l")
            points(wuv$uX[falseScale[i]], wuv$uY[falseScale[i]])
            points(wuv$vX[falseScale[i]], wuv$vY[falseScale[i]])
          }
          plotyyy <- 1
        }
        if(plotyy == "no"){
          cat("\n", "ok") 
          plotyyy <- 1
        }
        if(plotyy != "yes" & plotyy != "no"){cat("\n", "please enter 'yes' or 'no'")}
      }
    }
    if(length(falseScale) > 0){
      o2 <- o2[-falseScale]
      wuv <- wuv[-falseScale, ]
      diffVU <- diffVU[-falseScale]
    }
  }
  
  oY <- predict(sp, inx[o])                                                                   # Identify O y-values on the original PPG trace
  owDiff <- c()                                                                               # Find the x-axis difference between O and W
  for(i in 1:length(wuv$wX)){
    owDiff[i] <- wuv$wX[i] - inx[o[i]]
  }
  
  oDiff <- c()                                                                                # Find the x-axis difference between successive o_points
  for(i in 1:(length(inx[o])-1)){
    oDiff[i] <- inx[o[i+1]] - inx[o[i]]
  }
  
  ibi <- c()                                                                                  # Find the x-axis difference between successive W points,  
  for(i in 1:(length(wuv$wX))-1){                                                             # and designate as inter-beat intervals
    ibi[i] <- wuv$wX[i+1] - wuv$wX[i]
  }
  
  o2 <- o2[-length(o2)]                                                                       # Remove the last W (in case final beat is incomplete)
  wuv <- wuv[-nrow(wuv), ]
  diffVU <- diffVU[-length(diffVU)]  
  
  # plot(ibi)                                                                                 # Identify abnormal inter-beat intervals indicative of 
  # points(which(ibi > 1.25*median(ibi)),                                                     # an extended time between beats (e.g due to artefact)
  # ibi[which(ibi > 1.25*median(ibi))], pch = 19, col = "red")                                # These can be plotted (see left)
  
  endWaves <- c()                                                                             # Abnormal inter-beat intervals are defined as greater 
  for(i in 2:(length(ibi)-1)){                                                                # than the median interval by a factor of 1.3, and
    if(ibi[i] > 1.3*median(ibi) & ( ibi[i] > 2*ibi[i-1] | ibi[i] > 2*ibi[i+1])){              # greater than at least one of the values adjacent
      endWaves[i] <- i                                                                        # to it by a factor of 2
    }                                                                                         
  }                                                                                           # Beats preceding abnormal intervals are removed as they
  if(length(endWaves) > 0){                                                                   # will not have an O point to mark their end
    endWaves <- endWaves[!is.na(endWaves)]
    cat("\n", length(endWaves), "/", nrow(wuv), 
        "waves removed due to the end of the wave (next o point) 
         not being corrected for baseline")
    if(q == TRUE){
      plotyyy <- 0
      while(plotyyy == 0){
        plotyy <- readline(prompt = "Would you like to view? (enter yes or no)")
        if(plotyy == "yes"){
          for(i in 1:length(endWaves)){
            if(endWaves[i] == 1){
              plot(1:(samp*10), bc[1:(samp*10)], type = "l")
              points(wuv$wX[endWaves[i]], wuv$wY[endWaves[i]], pch =19, col = 'red')
            }else{
              plot((wuv$wX[endWaves[i]]-samp*5):(wuv$wX[endWaves[i]]+samp*5), 
                   bc[(wuv$wX[endWaves[i]]-samp*5):(wuv$wX[endWaves[i]]+samp*5)], type = "l")
              points(wuv$wX[endWaves[i]], wuv$wY[endWaves[i]], pch =19, col = 'red')
              points(wuv$wX[endWaves[i]+1], wuv$wY[endWaves[i]+1], pch = 19, col = 'red')
              points(wuv$wX, wuv$wY)
            }
          }
          plotyyy <- 1
        }
        if(plotyy == "no"){
          cat("\n", "ok") 
          plotyyy <- 1
        }
        if(plotyy != "yes" & plotyy != "no"){cat("\n", "please enter 'yes' or 'no'")}
      }
    }
    if(length(endWaves) > 0){
      o2 <- o2[-endWaves]
      wuv <- wuv[-endWaves, ]
      diffVU <- diffVU[-endWaves] 
    }
  }
  d <- cbind(wuv, diffVU, o2)                                                                 # Bind vectors with element number corresponding to beat number 
  dat <- list(d, ibi, oDiff)                                                                  # into a single dataframe before outputting 
  return(dat)
}


sep_beats <- function(odiff, bc, samp, wuv, wvlen, inx, o, ibi, scale = TRUE, q = FALSE, subset = FALSE, boundaries = NULL){   
  
  # Optional subsetting of Iso Waves:
  if(subset == T){
    #Find the Iso specific part of the data
    ab_ibi <- which(ibi > mean(ibi) + (4*sd(ibi)) | ibi < mean(ibi) - (4*sd(ibi))) #excluding ibi points 4 sds above mean
    ibi[ab_ibi] <- NA 
    ibi <- ibi[!is.na(ibi)]
    meds <- rollmedian(ibi, k = 19) #rolling median
    basemed <- mean(meds[1:50]) #baseline defined as the average of the first 50 rolling median points... I'm not sure if this is sensible
    halfs <- (basemed - min(meds))/2 #half diff between minimum and baseline
    post <- min(meds) + halfs   # the half way point (the y axis point that represents the half height)
    # Sometimes the minimum is incorrect (due to outlier ibi values):
    # Suggest first finding the minimum of meds and then finding a minimum ibi local to that (+/- 10 beats):
    meds_roi <- (which(meds == min(meds))[1] - 10) : (which(meds == min(meds))[1] + 10)
    # If no discernable minimum in the time series, it is possible a value far from the middle would be included, which would introduce negative values to meds_roi
    if( sum(which(meds_roi < 1)) > 0  | sum(which(meds_roi > length(ibi))) > 0 ){
      message("Warning: No discernable minimum in 2mg trace")
      Sys.sleep(10)
      meds_roi <- (round(length(ibi)/2) - 10) : (round(length(ibi)/2) + 10)
    }
    indx <- meds_roi[which(ibi[meds_roi] == min(ibi[meds_roi]))]    # indx is taken as the minimum
    # Check again that minimum is not too near to beginning / end of trace:
    if(indx < 15 | indx > (length(ibi) - 15)){
      message("Warning: No discernable minimum in 2mg trace")
      Sys.sleep(10)
      meds_roi <- (round(length(ibi)/2) - 10) : (round(length(ibi)/2) + 10)
    }
    pre_indx <- which(abs(ibi[1:indx] - post) < 0.5) + 0 # finds the values that are closest to the half point  # To ensure pre_indx is not < 75, have change ibi[1:indx] to ibi[75:indx], and 0 to 74...
    # If no values within 0.5 of the half-way point, extend the range to 1.5 and so on:
    a <- 0.5
    while(length(pre_indx) < 1){
      a <- a + 1
      pre_indx <- which(abs(ibi[1:indx] - post) < a) + 0 # Have to adjust for the fact that 1 = 75 now
    }
    in2 <- which(abs(pre_indx - indx) == min(abs(pre_indx - indx))) # take the one closest to the minimum
    pre_indx <- pre_indx[in2]
    
    # Find the closest to the half point from the points after the minimum
    post_indx <- indx + which(abs(ibi[indx:length(ibi)] - post) < 0.5)   # changed length(ibi) to 210
    a <- 0.5
    while(length(post_indx) < 1){
      a <- a + 0.25
      post_indx <- indx + which(abs(ibi[indx:length(ibi)] - post) < a)  # changed length(ibi) to 210
    }
    while(post_indx[1] < indx + 15){
      post_indx <- post_indx[-1]
    }
    post_indx <- post_indx[1]
    rm(meds, basemed, halfs, post, indx, in2)
    ppg_pre <- (which((ppg[,1]) == beat[,1][pre_indx - 1]))    # Before infusion
    ppg_post <- which((ppg[,1]) == beat[,1][post_indx + 1])    # After infusion wears off
    # Sometimes beat[, 1] will not contain the values specified to subset, in which case... 
    if(length(ppg_pre) < 1 | length(ppg_post) < 1){
      message("Subsetting failed: review of time series suggested")
      Sys.sleep(10)
      ppg_pre <- which((ppg[,1]) == beat[1,1])
      ppg_post <- which((ppg[,1]) == beat[nrow(beat),1])
    }
    
    # You might have to re-derive pre and post indx from ppg pre and post
    subs <- which(wuv$wX > ppg_pre & wuv$wX < ppg_post)
    cat(length(subs), "beats in subset (pre-cleaning)")
    #Sys.sleep(2)
    
    # Subset:
    wuv <- wuv[subs, ]
    # Subset ibi:
    ibi <- ibi[subs]
  }
  
  if(subset == "rep"){
    
    # Subset according to previous time series:
    subs <- which(wuv$wX > boundaries[1])
    subs_pre <- subs[1:boundaries[2]]
    # There's a chance there aren't enough beats remaining in the placebo time series, in which case take as many as possible:
    tmp <- sum(is.na(subs_pre))
    if(tmp > 0){
      subs <- subs[1:(length(subs_pre)-tmp)]
      boundaries[2] <- length(subs)
    }else{
      subs <- subs[1:boundaries[2]]
    }
    
    # Subset:
    wuv <- wuv[subs, ]
    # Subset:
    ibi <- ibi[subs]
  }
  
  
  # Redefine baseline corrected data:
  sourcedata <- bc[1:length(undetrended)]
  
  # Define a dataframe to contain individual waves (first column is the x-axis (in seconds) - currently set for bioradio data):
  pulse <- data.frame(seq((-141/(samp*10)), ((wvlen*15 -9)-142)/(samp*10), by = 1/(samp*10)))   
  
  afterO <- list()
  beforeO <- list()
  extra_long_wave <- c()
  for(i in 1:length(wuv$wX)){  
    
    splPolySub <- CubicInterpSplineAsPiecePoly((round(wuv$uX[i])-15):(round(wuv$uX[i]) + (wvlen+10)), sourcedata[(round(wuv$uX[i])-15):(round(wuv$uX[i]) + (wvlen+10))], "natural")
    # Turn into discrete form
    splSub <- predict(splPolySub, c(seq((wuv$uX[i]-14), (wuv$uX[i]+(wvlen+5)), 0.1)))
    # Make into dataframe:
    splSub <-  as.data.frame(splSub)
    splSub <- cbind(splSub, c(seq((wuv$uX[i]-14), (wuv$uX[i]+(wvlen+5)), 0.1)))
    colnames(splSub) <- c('y', 'x')  
    # Scale so that v-u = 1
    if(scale == TRUE){
      splSub$y <- splSub$y/(wuv$diffVU[i])      
    }
    # Adjust such that u = 0, v = 1 on y-axis
    yDiff <- splSub$y[141] 
    splSub$y <- splSub$y - yDiff
    
    if(scale == TRUE){
      # Find the x-value for each wave that corresponds to when it = 0.5 in height (this requires making a spline):
      splPolySub2 <- CubicInterpSplineAsPiecePoly(splSub$x, splSub$y, "natural")
      halfCross <- solve(splPolySub2, b = 0.5, deriv = 0)
      halfCross <- halfCross[which(abs(halfCross - wuv$wX[i]) == min(abs(halfCross - wuv$wX[i])))]    
      
      # Convert to discrete form again: (need to redefine splSub)
      splSub2 <- predict(splPolySub, c(seq((halfCross-14), (halfCross+(wvlen+5)), 0.1)))  
      splSub2 <-  as.data.frame(splSub2)
      splSub2 <- cbind(splSub2, c(seq((halfCross-14), (halfCross+(wvlen+5)), 0.1)))  
      colnames(splSub2) <- c('y', 'x') 
      
      # Scale again
      splSub2$y <- splSub2$y/(wuv$diffVU[i]) 
      # Adjust y-axis such that u = 0, v = 1
      yDiff <- wuv$uY[i] / wuv$diffVU[i]
      splSub2$y <- splSub2$y - yDiff
    }else{
      splSub2 <- splSub
    }
    
    # Find next_o
    afterO[[i]] <- which(splSub2$x > inx[o][min(which(inx[o] > wuv$wX[i]))])
    
    # Find values before the o of the wave itself 
    beforeO[[i]] <- which(splSub2$x < inx[o][max(which(inx[o] < wuv$wX[i]))])
    
    # Correct such that x column and wave column are correctly aligned
    splSub3 <- c()
    for(j in 1:nrow(splSub2)){
      splSub3[j+1] <- splSub2$y[j]
    }
    
    # If splSub3 and nrow(pulse) are the same length, you need only adjust afterO
    if(length(splSub3) == nrow(pulse)){
      if(length(afterO[[i]]) > 0){
        diff2 <- abs(length(splSub3) - max(afterO[[i]]))
        for(j in 1:diff2){
          afterO[[i]] <- c(afterO[[i]], (max(afterO[[i]]) + 1))
        }
      }
    }
    
    # Or
    # Adjust such that splSub3 is the same length as pulse
    if(length(splSub3) > nrow(pulse)){
      diff <- length(splSub3) - nrow(pulse)
      len <- length(splSub3)
      splSub3 <- splSub3[-((len - (diff-1)):len)]
      if(diff > 1){     # must correct the afterO values so that they also do not contain values beyond the length of splSub3 (include case where length of afterO[[i]] is one so the code works...)
        if(length(afterO[[i]]) > 1 ){
          afterO[[i]] <- afterO[[i]][-(which(afterO[[i]] > length(splSub3)))]  #afterO[[i]][1:(which(afterO[[i]] == (len - (diff-1))) - 1) 
        }else{
          afterO[[i]] <- afterO[[i]][-(which(afterO[[i]] > length(splSub3)))]
        }
      }
    }
    
    
    if(length(splSub3) < nrow(pulse)){
      diff <- nrow(pulse) - length(splSub3)  # diff = how much longer pulse is than splSub3...
      splSub3 <- c(splSub3, rep(NA, diff))   # The extra rows in pulse are thus made up by NAs
      if(length(afterO[[i]]) > 0){
        diff2 <- length(splSub3) - max(afterO[[i]])   # This = how many elements at the end of splSub3 are after O.
        # This for loop adds values to after_o to make up the difference
        for(j in 1:diff2){
          afterO[[i]] <- c(afterO[[i]], (max(afterO[[i]]) + 1))
        }
      }
    }
    
    # Add column to dataframe
    pulse <- cbind(pulse, splSub3)
  }
  
  for(i in 1:(ncol(pulse) -1)){ 
    colnames(pulse)[i+1] <- paste("wave", i, sep = "_")       
  }
  colnames(pulse)[1] <- "x"
  
  # Remove any values after O for each wave:
  for(i in 2:(ncol(pulse))){
    pulse[, i][afterO[[(i-1)]][-1]] <- NA  
  }
  
  # Remove values before O before each wave:
  for(i in 2:(ncol(pulse))){
    pulse[, i][beforeO[[(i-1)]][-1]] <- NA  
  }
  
  # Remove any extra long waves (i.e where a distance of 2 Os has been counted as one wave):
  # find average length of waves (without NAs) in pulse
  wavelengths <- c()
  for(i in 2:ncol(pulse)){
    wavelengths[i] <- length(pulse[, i][!is.na(pulse[, i])])  
  }
  wavelengths <- wavelengths[!is.na(wavelengths)]
  
  
  ########## Cleaning: #########
  
  extra_long_wave <- c()
  for(i in 2:length(wavelengths)){
    if(wavelengths[i] > (mean(wavelengths) + sd(wavelengths)) & wavelengths[i] > 1.8*(wavelengths[i-1]) ){
      extra_long_wave[i] <- i
    }
  }
  extra_long_wave <- extra_long_wave[!is.na(extra_long_wave)]
  if(length(extra_long_wave) > 0){
    #rejected_waves <- rejected_waves
    cat("\n", length(extra_long_wave), "/", (ncol(pulse)-1), "waves removed for being abnormally long")
    if(q == TRUE){
      plotyyy <- 0
      while(plotyyy == 0){
        plotyy <- readline(prompt = "Would you like to view? (enter yes or no)")
        if(plotyy == "yes"){
          for(i in 1:length(extra_long_wave)){
            plot((wuv$wX[extra_long_wave[i]]-samp*2):(wuv$wX[extra_long_wave[i]]+samp*2), bc[(wuv$wX[extra_long_wave[i]]-samp*2):(wuv$wX[extra_long_wave[i]]+samp*2)], type = "l")
            points(wuv$wX[extra_long_wave[i]], wuv$wY[extra_long_wave[i]], pch =19, col = 'red')
          }
          plotyyy <- 1
        }
        if(plotyy == "no"){
          cat("\n", "ok") 
          plotyyy <- 1
        }
        if(plotyy != "yes" & plotyy != "no"){cat("\n", "please enter 'yes' or 'no'")}
      }
    }
    wuv <- wuv[-(extra_long_wave[!is.na(extra_long_wave)]-1), ]
    pulse <- pulse[, -(extra_long_wave[!is.na(extra_long_wave)])]  
    ibi <- ibi[-(extra_long_wave[!is.na(extra_long_wave)]-1)]
  }
  
  
  
  extra_short_waves <- c()
  for(i in 2:length(wavelengths)){
    if(wavelengths[i] < mean(wavelengths[-1]) - 3*sd(wavelengths[-1])){
      extra_short_waves[i] <- i
    }
  }
  extra_short_waves <- extra_short_waves[!is.na(extra_short_waves)]
  if(length(extra_short_waves) > 0){
    cat("\n", length(extra_short_waves), "/", (ncol(pulse)-1), "waves removed for being abnormally short in duration")
    if(q == TRUE){
      plotyyy <- 0
      while(plotyyy == 0){
        plotyy <- readline(prompt = "Would you like to view? (enter yes or no)")
        if(plotyy == "yes"){
          for(i in 1:length(extra_short_waves)){
            plot((wuv$wX[extra_short_waves[i]]-samp*2):(wuv$wX[extra_short_waves[i]]+samp*2), bc[(wuv$wX[extra_short_waves[i]]-samp*2):(wuv$wX[extra_short_waves[i]]+samp*2)], type = "l")
            points(wuv$wX[extra_short_waves[i]], wuv$wY[extra_short_waves[i]], pch =19, col = 'red')
            points()
          }
          plotyyy <- 1
        }
        if(plotyy == "no"){
          cat("\n", "ok") 
          plotyyy <- 1
        }
        if(plotyy != "yes" & plotyy != "no"){cat("\n", "please enter 'yes' or 'no'")}
      }
    }
    wuv <- wuv[-(extra_short_waves-1), ]
    pulse <- pulse[, -(extra_short_waves)]  
    ibi <- ibi[-(extra_short_waves-1)]
  }
  
  
  double_segments <- c()
  double_peaks_wave <- list()
  for(i in 2:ncol(pulse)){
    # Make a spline to identify inflection points:
    wave <- pulse[, i][!is.na(pulse[, i])]
    sfunction <- splinefun(1:length(wave), wave, method = "natural")
    spline1 <- sfunction(seq(1, length(wave)), deriv = 0)
    splinePoly <- CubicInterpSplineAsPiecePoly(1:length(wave), spline1, "natural") 
    inflexX <- solve(splinePoly, b = 0, deriv = 1)
    inflexY <- predict(splinePoly, inflexX)
    peaks <- which(inflexY > mean(inflexY) + 1.7*sd(inflexY))
    
    # plot(splinePoly)
    # points(inflexX, inflexY)
    
    # If there are more than one inflection points above threshold and the x-axis distance between the first and last is greater than 100, assume a double segment
    if(length(peaks) > 0){
      if( length(peaks) > 1 & ((inflexX[peaks][length(inflexX[peaks])] - inflexX[peaks[1]]) > 100) ){
        double_segments[i] <- i 
        double_peaks_wave[[i]] <- wave
      }
    }
  }
  double_segments <- double_segments[!is.na(double_segments)]
  if(length(double_peaks_wave) > 1){double_peaks_wave <- double_peaks_wave[-(which(sapply(double_peaks_wave, is.null)))]}        
  if(length(double_segments) > 0){
    cat("\n", length(double_segments), "/", (ncol(pulse)-1), "waves removed for containing two systolic peaks")
    if(q == TRUE){
      plotyyy <- 0
      while(plotyyy == 0){
        plotyy <- readline(prompt = "Would you like to view? (enter yes or no)")
        if(plotyy == "yes"){
          for(i in 1:length(double_segments)){
            plot(1:length(double_peaks_wave[[i]]), double_peaks_wave[[i]], t = "l")
          }
          plotyyy <- 1
        }
        if(plotyy == "no"){
          cat("\n", "ok") 
          plotyyy <- 1
        }
        if(plotyy != "yes" & plotyy != "no"){cat("\n", "please enter 'yes' or 'no'")}
      }
    }
    wuv <- wuv[-(double_segments-1), ]
    pulse <- pulse[, -(double_segments)]  
    ibi <- ibi[-(double_segments-1)]
  }
  
  
  # Waves that enter a second systolic upstroke
  systolic_endings <- c()
  for(i in 2:ncol(pulse)){
    wave <- pulse[, i] 
    wave <- wave[!is.na(wave)]
    if( wave[length(wave)] > 0.25 | (length(wave) > mean(wavelengths) & (max(wave[round(0.75*length(wave)):length(wave)]) > 0.8)) ){            # Adding additional criterion here: if the wave is above average length and the last quarter has a value above 0.8, consider it also a second systolic peak
      systolic_endings[i] <- i
    }   
  }
  systolic_endings <- systolic_endings[!is.na(systolic_endings)]
  if(length(systolic_endings) > 0){
    cat("\n", length(systolic_endings), "/", (ncol(pulse)-1), "waves removed for including the following systolic upstroke")
    if(q == TRUE){
      plotyyy <- 0
      while(plotyyy == 0){
        plotyy <- readline(prompt = "Would you like to view? (enter yes or no)")
        if(plotyy == "yes"){
          for(i in 1:length(systolic_endings)){
            plot(pulse[, systolic_endings[i]])
          }
          plotyyy <- 1
        }
        if(plotyy == "no"){
          cat("\n", "ok") 
          plotyyy <- 1
        }
        if(plotyy != "yes" & plotyy != "no"){cat("\n", "please enter 'yes' or 'no'")}
      }
    }
    wuv <- wuv[-(systolic_endings-1), ]
    pulse <- pulse[, -systolic_endings]
    ibi <- ibi[-(systolic_endings-1)]
  }
  
  # Find average wave:
  average_wave <- find_average(p = pulse, ao = afterO)
  
  # Waves that fall significantly below O:
  drops_below_o <- c()
  for(i in 2:ncol(pulse)){
    wave <- pulse[, i][!is.na(pulse[, i])]
    thd <- 4   # threshold subject to change (default 2)
    if((min(wave) < wave[1]*thd | min(wave) < -0.3 | min(wave) < min(average_wave[!is.na(average_wave)])*thd) & min(wave) < 0){  
      drops_below_o[i] <- i
    }
  }
  drops_below_o <- drops_below_o[!is.na(drops_below_o)]
  if(length(drops_below_o) > 0){
    cat("\n", length(drops_below_o), "/", (ncol(pulse)-1), "waves removed for having values significantly below baseline")
    if(q == TRUE){
      plotyyy <- 0
      while(plotyyy == 0){
        plotyy <- readline(prompt = "Would you like to view? (enter yes or no)")
        if(plotyy == "yes"){
          for(i in 1:length(drops_below_o)){
            plot(pulse[, drops_below_o[i]])
          }
          plotyyy <- 1
        }
        if(plotyy == "no"){
          cat("\n", "ok") 
          plotyyy <- 1
        }
        if(plotyy != "yes" & plotyy != "no"){cat("\n", "please enter 'yes' or 'no'")}
      }
    }
    wuv <- wuv[-(drops_below_o-1), ]        # Remember, the pulse has the 1st column, meaning pulse[, drops_below_o] = wuv[drops_below_o - 1, ] (hence the -1 here) 
    pulse <- pulse[, -drops_below_o]    
    ibi <- ibi[-(drops_below_o-1)]
  }
  
  
  # Recalculate average after removing the most aberrant waves:
  average_wave <- find_average(p = pulse, ao = afterO)
  
  # Find the SD wave:
  sd_wave <- find_sd(p = pulse, ao = afterO)
  
  # Remove waves that have a high SD of residuals:
  resid_sd <- c()
  for(i in 2:ncol(pulse)){
    residuals <- average_wave[142:length(average_wave)] - pulse[, i][142:length(pulse[, i])]
    resid_sd[i] <- sd(residuals[!is.na(residuals)][-c(1:100)])  # removed 1:100 since at this point all values are close to the average on the systolic upstroke
  }
  resid_sd <- sqrt(resid_sd)  # equivalent of squaring, since the resifual values are fractions when scaled
  resid_sd <- resid_sd[-1] # remove the NA
  # thld <- mean(resid_sd) - sd(resid_sd)*0.5   # threshold could be adjusted (has been increased from 2)
  thld <-  0.35  
  hrsd_waves <- which(resid_sd > thld) + 1  # +1 to account for the removed NA value
  if(length(hrsd_waves) > 0){
    cat("\n", length(hrsd_waves), "/", (ncol(pulse)-1), "waves removed for having high residual SD")
    if(q == TRUE){
      plotyyy <- 0
      while(plotyyy == 0){
        plotyy <- readline(prompt = "Would you like to view? (enter yes or no)")
        if(plotyy == "yes"){
          for(i in 1:length(hrsd_waves)){
            plot(pulse[, hrsd_waves[i]])
          }
          plotyyy <- 1
        }
        if(plotyy == "no"){
          cat("\n", "ok") 
          plotyyy <- 1
        }
        if(plotyy != "yes" & plotyy != "no"){cat("\n", "please enter 'yes' or 'no'")}
      }
    }
    wuv <- wuv[-(hrsd_waves-1), ]
    pulse <- pulse[, -hrsd_waves]
    ibi <- ibi[-(hrsd_waves-1)]
  }
  
  
  # Waves that are beyond 5 SDs away from the average:
  outlier_waves <- c()
  for(i in 2:ncol(pulse)){
    breaches <- c()
    wave <- pulse[, i]
    for(j in 142:length(wave)){   # look only after w
      if(is.na(pulse[j, i]) == FALSE & is.na(average_wave[j] + 4*sd_wave[j]) == FALSE){
        if(pulse[j, i] > (average_wave[j] + 4*sd_wave[j]) | pulse[j, i] < (average_wave[j] - 4*sd_wave[j])){
          breaches[j] <- 1
        }
      }
    }
    if(sum(breaches[!is.na(breaches)]) > 0){
      outlier_waves[i] <- i
    }
  }
  outlier_waves <- outlier_waves[!is.na(outlier_waves)]
  if(length(outlier_waves) > 0){
    cat("\n", length(outlier_waves), "/", (ncol(pulse)-1), "waves removed for having values beyond 5 SDs from the average")
    if(q == TRUE){
      plotyyy <- 0
      while(plotyyy == 0){
        plotyy <- readline(prompt = "Would you like to view? (enter yes or no)")
        if(plotyy == "yes"){
          for(i in 1:length(outlier_waves)){
            plot(pulse[, outlier_waves[i]])
          }
          plotyyy <- 1
        }
        if(plotyy == "no"){
          cat("\n", "ok") 
          plotyyy <- 1
        }
        if(plotyy != "yes" & plotyy != "no"){cat("\n", "please enter 'yes' or 'no'")}
      }
    }
    wuv <- wuv[-(outlier_waves-1), ]
    pulse <- pulse[, -outlier_waves]
    ibi <- ibi[-(outlier_waves-1)]
  }
  
  # Final recalculation of average wave:
  average_wave <- find_average(p = pulse, ao = afterO) 
  
  # print beats carried over:
  cat(nrow(wuv), "beats carried over for 2mg subsetting (post-cleaning)")
  
  # Save rejected waves:
  rejects <- list(extra_long_wave, extra_short_waves, double_segments, systolic_endings, drops_below_o, hrsd_waves, outlier_waves)
  
  dat <- list(average_wave, pulse, wuv, rejects)
  if(subset == T){dat <- list(average_wave, pulse, wuv, rejects, c(ppg_pre, nrow(wuv)), ibi)}
  
  return(dat)
}


find_average <- function(p, ao){
  
  # Find the last 10 values of each wave:
  last10 <- list()
  for(i in 2:(ncol(p))){
    last10[[i-1]] <- p[, i][!is.na(p[, i])][(length(p[, i][!is.na(p[, i])])-10):(length(p[, i][!is.na(p[, i])]))]
  }
  
  # Find the average last 10 values
  avLast10 <- c()
  for(i in 1:10){   
    rowVec <- c()
    for(j in 1:length(last10)){      
      rowVec[j] <- last10[[c(j, i)]] 
    }
    avLast10[i] <- mean(rowVec[!is.na(rowVec)])   
  }
  
  # Find the consecutive y-axis differences (gradient essentially) of last 10 values:
  avDiff <- c()
  for(i in 1:9){                    # 9 here since 10 values will give 9 values for differences between them
    avDiff[i] <- avLast10[i+1] - avLast10[i]
  }
  
  # Simply calculating the average wave by averaging each row of the dataframe doesn't work at the end of the wave, 
  # since as waves end, they no longer factor in to the calculation of the average. The average would be erroneously 
  # drawn out until the end of the last wave. Therefore, a clone dataframe of p is used to continue the trajectories
  # of waves after they have ended (by using the average gradient above) and thus make the end of the average wave
  # a more accurate approximation of the true average. 
  
  # Replace NA values with continuing downward gradient in a clone dataframe:
  p2 <- p
  for(i in 2:(ncol(p2))){
    for(j in 1:length(p2[, i][ao[[(i-1)]][-1]])){
      p2[, i][ao[[(i-1)]][-1]][j] <- p2[, i][ao[[(i-1)]][1]] + j*mean(avDiff[1:5])
    }
  }
  
  # Calculate the average wave by row:
  avWav <- c()
  sdWav <- c()
  medWav <- c()
  for(i in 1:nrow(p2)){
    rowVec <- c()
    for(j in 2:(ncol(p2))){
      rowVec[j-1] <- p2[i, j] 
    }
    avWav[i] <- mean(rowVec[!is.na(rowVec)])   
    sdWav[i] <- sd(rowVec[!is.na(rowVec)]) 
    medWav[i] <- median(rowVec[!is.na(rowVec)])
  }
  
  # Find where the end of the average wave should end: 
  # This is done by finding the average (mode) y-value of the final value of each wave (before trajectory continuation)
  end <- c()
  for(i in 2:ncol(p)){
    end[i-1] <- p2[, i][ao[[(i-1)]][1]]
  }
  # define mode function
  mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  # Made this into a for loop because there can be waves that have no values after o, and if all waves are like that, end will be null
  if(length(end) > 1){     
    avEnd <- mode(round(end[!is.na(end)], digits = 2))
    # Remove elements of average wave after where its end should be:
    belowO <- which(avWav < avEnd)
    # Find the first point in the last quarter of the wave that goes below baseline
    if(sum(which(belowO > (length(avWav)/(4/3)))) > 0){
      b. <- belowO[min(which(belowO > (length(avWav)/(4/3))))]
    }else{
      b. <- NA
    }
    # If the wave doesn't go below baseline, don't remove any 
    if(is.na(b.) == FALSE){
      # otherwise, remove those elements from the average:
      avWav[b.:length(avWav)] <- NA
    }
  }
  
  # Find the median of the first x-values - make that the start of average:
  first <- c()
  for(i in 2:ncol(p)){
    stpqw <- p[, i][1:200]
    stq <- stpqw[is.na(stpqw)]
    first[i-1] <- length(stq) + 1
  }
  start <- median(first)-2
  avWav[1:start] <- NA
  
  # From the last value moving backwards, if it's below the value of the first avWave value, remove it... 
  avWave <- as.vector(avWav)
  
  avWave_tmp <- avWave[!is.na(avWave)]
  
  while(avWave_tmp[length(avWave_tmp)] < avWave_tmp[1]){
    avWave_tmp <- avWave_tmp[-length(avWave_tmp)]
    avWave[length(avWave_tmp) + which(!is.na(avWave))[1]] <- NA
  }
  
  return(avWave)
}



find_sd <- function(p, ao){
  
  # Find the last 10 values of each wave:
  last10 <- list()
  for(i in 2:(ncol(p))){
    last10[[i-1]] <- p[, i][!is.na(p[, i])][(length(p[, i][!is.na(p[, i])])-10):(length(p[, i][!is.na(p[, i])]))]
  }
  
  # Find the average last 10 values
  mean_last10 <- c()
  for(i in 1:10){   
    row_vector <- c()
    for(j in 1:length(last10)){      
      row_vector[j] <- last10[[c(j, i)]] 
    }
    mean_last10[i] <- mean(row_vector[!is.na(row_vector)])   
  }
  
  # Find the consecutive y-axis differences (gradient essentially) of last 10 values:
  mean_diff_last10 <- c()
  for(i in 1:9){                    # 9 here since 10 values will give 9 values for differences between them
    mean_diff_last10[i] <- mean_last10[i+1] - mean_last10[i]
  }
  
  # Replace NA values with continuing downward gradient in a clone dataframe:
  p_for_finding_average <- p
  for(i in 2:(ncol(p_for_finding_average))){
    for(j in 1:length(p_for_finding_average[, i][ao[[(i-1)]][-1]])){
      p_for_finding_average[, i][ao[[(i-1)]][-1]][j] <- p_for_finding_average[, i][ao[[(i-1)]][1]] + j*mean(mean_diff_last10[1:5])
    }
  }
  
  # Calculate the average wave by row:
  average_wave <- c()
  sd_wave <- c()
  median_wave <- c()
  for(i in 1:nrow(p_for_finding_average)){
    row_vector <- c()
    for(j in 2:(ncol(p_for_finding_average))){
      row_vector[j-1] <- p_for_finding_average[i, j] 
    }
    average_wave[i] <- mean(row_vector[!is.na(row_vector)])   
    sd_wave[i] <- sd(row_vector[!is.na(row_vector)]) 
    median_wave[i] <- median(row_vector[!is.na(row_vector)])
  }
  
  return(as.vector(sd_wave))
}



diast_pk <- function(avw, sr, scale = F, dias_param = NULL){
  # Find the diastolic peak on the average wave to inform OSND finding (also some adjusment of x-values for removal of NA values):
  avw <- avw[!is.na(avw)]
  
  # Need to find new W position (0.5) after removing NAs
  if(scale == TRUE){
    xShift <- which(abs(avw-0.5) == min(abs(avw - 0.5)))
  }else{
    xShift <- which.min(abs(avw)) 
  }
  avWavePoly <- CubicInterpSplineAsPiecePoly(1:length(avw), avw, "natural")
  avInflexX <- solve(avWavePoly, b = 0, deriv = 1)
  avInflexY <- predict(avWavePoly, avInflexX)
  
  # plot(avWavePoly)
  # points(avInflexX, avInflexY) 
  
  # Specify limitations for where the diastolic peak can first be found i.e between 120:230 on x-axis, and below 1 on y-axis:
  #peaks <- order(avInflexY[which(avInflexX < 215 & avInflexX > 120 & avInflexY < 1)], decreasing = TRUE)  
  #diastPk <- avInflexX[which(avInflexX < 215 & avInflexX > 120 & avInflexY < 1)][peaks[1]]
  
  # OR:
  # Find any peaks (above 0):
  peaks_above_0 <- which(avInflexY > 0)
  
  # If a fitted wave is being assessed, find out what the dias_param corresponds to:
  if(!is.null(dias_param)){
    # peaks_above_0[1] can either be O or S - make sure it is being assigned to S:
    if( avInflexY[peaks_above_0[1]] <  mean(avw) ){
      peaks_above_0 <- peaks_above_0[-1]
    }
    dias_param <- dias_param + avInflexX[peaks_above_0[1]]
    diastPk <-  avInflexX[peaks_above_0[which( abs(avInflexX[peaks_above_0] - dias_param) == min(abs(avInflexX[peaks_above_0] - dias_param)))]] 
    # If the value assigned is just the max of the wave (i.e systolic peak, default to sr*500):
    if(avInflexY[peaks_above_0[which( abs(avInflexX[peaks_above_0] - dias_param) == min(abs(avInflexX[peaks_above_0] - dias_param)))]] - max(avw) < 2){
      diastPk <- 5*sr
    }
    peaks <- peaks_above_0
    fitted <- 1
  }else{
    fitted <- 0
  }
  
  if(fitted != 1){
    # Order peaks above 0 in descending order
    peaks <- order(avInflexY[peaks_above_0], decreasing = TRUE)
    # This ensures peaks refers to the correct avInflexX points
    peaks <- peaks_above_0[peaks]
    # But reducing them all as required so that there is no unnamed amount of non-peaks at the beginning of the wave, might be neater
    # So make the first peak 1:
    dif <- peaks[1] - 1
    peaks <- peaks - dif
    # Correct avInflexX / avInflexY:
    if(dif > 0){
      avInflexX <- avInflexX[-c(1:dif)]
      avInflexY <- avInflexY[-c(1:dif)] 
    }
    # Remove 0 or negative values if they have resulted from the above subtraction:
    if(sum(peaks < 1) > 0){
      peaks <- peaks[-which(peaks < 1)]   
    }
    
    # Take the second highest inflection point as the diastolic peak (unless diastolic is taller than systolic, then then take the highest)
    if(peaks[1] != 1 & length(peaks) > 1){diastPk <- avInflexX[peaks[1]]}else{diastPk <- avInflexX[peaks[2]]}
    # unless it is unreasonably high i.e a diastolic to systolic peak ratio of > 0.75 (the threshold for unreasonably high must be dependent on how close the peak is to S on the x-axis (sd_time) - lower threshold for closer)   
    if(!is.na(diastPk) & diastPk == avInflexX[peaks[2]]){
      ai. <- avInflexY[peaks[2]] / avInflexY[peaks[1]]
      sd_time <- avInflexX[peaks[2]] - avInflexX[peaks[1]]
      # We move forward through the peaks until finding one that is not unreasonably high (or unreasonably low):
      while( (ai. > 0.45 & sd_time < 75) | ai. > 0.7 | ai. < 0.05 ){   # these tresholds are difficult to get right
        # Assign the previous Dpeak as old, and remove it from the peaks list
        old <- peaks[2]
        peaks <- peaks[-2]
        # Replace the presumed D peak with one that later in time
        # Is the new 2nd highest peak earlier in time? If so we skip it and assign the 3rd value as the dpeak, if not, we assign it as the dpeak. 
        # But if there is no second highest peak now (i.e only systolic remains), se default value:
        if(is.na(peaks[2])){
          diastPk <- 5*sr
          ai. <- 0.5
          sd_time <- 100
        }else{
          if(peaks[2] < old){
            peaks <- peaks[-2]
            if(is.na(peaks[2])){
              diastPk <- 5*sr
              ai. <- 0.5
              sd_time <- 100
            }else{
              diastPk <- avInflexX[peaks[2]] 
              ai. <-  avInflexY[peaks[2]] / avInflexY[peaks[1]]
              sd_time <- avInflexX[peaks[2]] - avInflexX[peaks[1]]
            }
          }else{
            diastPk <- avInflexX[peaks[2]] 
            ai. <- avInflexY[peaks[2]] / avInflexY[peaks[1]]
            sd_time <- avInflexX[peaks[2]] - avInflexX[peaks[1]]
          }
        }
        # If no valid peaks remain, choose the default estimation:
        if(is.na(ai.)){
          diastPk <- 5*sr
          ai. <- 0.5
          sd_time <- 100
        }
      }
    }
  }
  
  # diastPk will be NA for class 3 waveforms, in which case set a default value
  if(is.na(diastPk) | diastPk < avInflexX[peaks[1]]){
    diastPk <- 5*sr
  }
  return(c(diastPk, xShift))
}



osnd_of_average <- function(aw, dp, diff, sr, plot = TRUE){
  
  switch <- 0
  aw <- aw[!is.na(aw)]
  
  avWavPoly <- CubicInterpSplineAsPiecePoly(1:length(aw), aw, "natural")
  sfunction <- splinefun(1:length(aw), aw, method = "natural")
  d1Wav <- sfunction(1:length(aw), deriv = 1)
  d1WavPoly <- CubicInterpSplineAsPiecePoly(1:length(aw), d1Wav, "natural") 
  
  # plot(avWavPoly)
  
  # Find inflexion points on d1WavPoly
  d1InflxX <- solve(d1WavPoly, b = 0, deriv = 1)
  d1InflxY <- predict(d1WavPoly, d1InflxX)
  
  # plot(d1WavPoly)
  # points(d1InflxX, d1InflxY)
  
  # Find OSND
  wavInflxX <- solve(avWavPoly, b = 0, deriv = 1)
  wavInflxY <- predict(avWavPoly, wavInflxX)
  
  # Finding notch based on x-axis:
  # Find inflexion point closest to where the notch usually is (aka 75-80)
  notchRange <- which(d1InflxX > (3.104572*sr - diff) & d1InflxX < dp) # 3.104572 used to be 3.5! #  dp used to be 5*sampling rate!
  # If there is no inflexion point detected within the notch range, this could be because there is a plateu rather than a peak
  # In this case, taking the mean value of the notch range boundaries gives a reasonable approximation
  if(length(notchRange) < 1 | (length(notchRange) == 1 & d1InflxY[notchRange][1] < -0.02)){
    new.n <- ((3.104572*sr) + dp)/2
    # In case new.n is greater than the number of datapoints (or very close to the end), set D to last inflection point on derivative (assuming the wave is very short)
    if(new.n > length(aw) | length(aw) - new.n < (5*(length(aw)/100)) ){
      new.n <- d1InflxX[length(d1InflxX)]
    }
  }else{
    if(length(notchRange) != 1){
      a. <- which(d1InflxY[notchRange] == max(d1InflxY[notchRange]))
      # In cases where the renal peak is higher on 1st deriv than the notch peak, make sure the notch peak is limited by x-axis
      while(d1InflxX[notchRange[a.]] < sr*3){   # was 115 instead of sr*3
        b. <- 2
        a. <- order(d1InflxY[notchRange], decreasing = TRUE)[b.]
        b. <- 3
      }
      # Make sure the 1st peak is not the notch:
      if(which(d1InflxY == max(d1InflxY)) == notchRange[a.]){
        if(length(notchRange[which(notchRange > notchRange[a.])]) > 0){
          notchRange <- notchRange[which(notchRange > notchRange[a.])] 
        }
      }
      if(is.na(d1InflxX[notchRange[a.]])){  # Added this due to bug with participant 7
        a. <- which(d1InflxY[notchRange] == max(d1InflxY[notchRange]))
        new.n <- d1InflxX[notchRange[a.]-1]
      }else{
        new.n <- d1InflxX[notchRange[a.]]
      }
    }else{
      new.n <- d1InflxX[notchRange[1]]
    }
  }
  new.ny <- predict(avWavPoly, new.n)
  # plot(d1WavPoly)
  # points(new.n, predict(d1WavPoly, new.n))
  
  
  # If there are no inflection points before new.n, assume it is incorrect and move to the next inflection point that has a lower y-value than it:
  if(sum(which(wavInflxX < new.n)) < 1){ 
    new.n <- wavInflxX[min(which(wavInflxX > new.n & wavInflxY < new.ny))]
    new.ny <- predict(avWavPoly, new.n)
  } 
  
  
  # After having found the notch on all waves, you can see if there are inflexion points either side:
  # If inflexion point before is lower and inflexion point after is higher (on y axis), this must mean a second peak aka canonical wave
  # Thus if this criterion is fulfilled you can create N and D separately 
  
  
  if(length(wavInflxX) > 1 & wavInflxY[max(which(wavInflxX < new.n))] < new.ny){
    new.n <- wavInflxX[max(which(wavInflxX < new.n))]
    new.ny <-  predict(avWavPoly, new.n)
    
    # Take the inflection point after the notch as the d peak. If there is no inflection point after the notch, take instead the next inflection point on deriv1
    if(sum(wavInflxX > new.n) > 0){
      d. <- wavInflxX[min(which(wavInflxX > new.n))]
      d.y <- predict(avWavPoly, d.)
    }else{
      d. <- d1InflxX[min(which(d1InflxX > new.n))]
      d.y <- predict(avWavPoly, d.)
    }
    switch <- 1
  }
  
  
  
  if(plot == TRUE){
    plot(avWavPoly)
    points(wavInflxX, wavInflxY)
    points(new.n, new.ny, col = "red")
    if(switch == 1){
      points(d., d.y, col = "blue")
    }
  }
  
  
  # Find W: the max inflection point on first deriv:
  w. <- d1InflxX[which( d1InflxY ==  max(d1InflxY) & d1InflxX[which(d1InflxY ==  max(d1InflxY))] < new.n)]
  if(length(w.) < 1){
    w. <- 1
  }
  w.y <- predict(avWavPoly, w.)
  if(plot ==TRUE){points(w., w.y)}
  
  # Find U and V:
  # Find half the height of w (on derivative y-axis)
  hhaw <- max(d1InflxY)/2
  # Find u and v for derivative:
  halfHeights <- solve(d1WavPoly, b = hhaw)
  #halfHeightsY <- predict(d1WavPoly, halfHeights)
  
  # If only one half height detected, assume the segment had quite a high O value such that O was higher than U:
  if(length(halfHeights) < 2){
    halfHeights[2] <- halfHeights[1]
    halfHeights[1] <- solve(d1WavPoly, b = d1Wav[1])[1]
  }
  
  # If more than one half height detected:
  if(length(halfHeights) > 2){
    # Find the distance between each detected half height, compare their xvals to the peak, and keep only the two that have the most similar distance to the peak... 
    # bear in mind, one has to be either side of w...
    a <- halfHeights - w.
    postW <- which(a > 0)
    preW <- which(a < 0)
    # If no points are identified before W (e.g if the wave is too steep to begin with, take the first value as you)
    if(length(preW) < 1){
      print("No valid U found, first value taken as proxy")
      halfHeights <- c(1, min(halfHeights[postW]))
    }else{
      halfHeights <- c(max(halfHeights[preW]), min(halfHeights[postW])) 
    }
  }
  
  u <- halfHeights[1]
  v <- halfHeights[2] 
  # Find u and v y-values for original wave:
  uvY <- predict(avWavPoly, halfHeights)
  uY <- uvY[1]
  vY <- uvY[2]
  if(plot == TRUE){
    points(u, uY)
    points(v, vY) 
  }
  
  # If there are inflection points before w, use the maximum one as 0:
  if(length(which(wavInflxX < w.)) > 0){
    o. <- wavInflxX[max(which(wavInflxX < w.))]
  }else{
    # If there are no inflection points before w, check if there are any in the first deriv before w:
    if(length(which(d1InflxX < w.)) < 1){
      # If not, assign O to the first value
      o. <- 1
    }else{
      # If there are, make sure the max inflection point is not greater than 0 on the original y-axis (which would be higher than U):
      if((predict(avWavPoly, d1InflxX[max(which(d1InflxX < w.))]) > 0)){
        # If the max inflection point is above 0, assign O to the first value
        o. <- 1
      }else{
        # If the max inflection point is not above 0, assign O to it:
        o. <- d1InflxX[max(which(d1InflxX < w.))]
      }
    }
  }
  o._yval <- predict(avWavPoly, o.)
  if(plot == TRUE){points(o., o._yval, pch = 19)}
  
  
  ## Find S:
  # Find new S:
  if( (v - w.) > 100 ){    # if v is erroneous, then this part will not generate a viable alternative to s1
    s2Y <- max(wavInflxY)
    s2 <- wavInflxX[which(wavInflxY == max(wavInflxY))]
  }else{
    s2 <- w. + 2*(abs(v - w.))
    s2Y <- predict(avWavPoly, s2)
  }        
  # points(s2, s2Y)
  # Define old s:
  s1Y <- max(wavInflxY)
  s1 <- wavInflxX[which(wavInflxY == max(wavInflxY))]
  # Decide which S to use...
  if((s1 - w.) < (s2 - w.)){
    s. <- s1
    s.y <- s1Y
  }else{
    s. <- s2
    s.y <- s2Y
  }
  if(plot == TRUE){points(s., s.y, pch = 19)} 
  
  if(switch == 0){
    x <- c(o., s., new.n, new.n)
  }else{
    x <- c(o., s., new.n, d.)
  }
  
  if(switch == 0){
    y <- c(o._yval, s.y, new.ny, new.ny)
  }else{
    y <- c(o._yval, s.y, new.ny, d.y)
  }
  
  osnd <- data.frame(x, y)
  return(osnd)
}


feature_extract <- function(oa, p, pw){
  # S_values:
  s_vals <- c()
  for(i in 1:length(oa)){
    s_vals[i] <- oa[[i]]$y[2]
  }
  
  # N_values:
  n_vals <- c()
  for(i in 1:length(oa)){
    n_vals[i] <- oa[[i]]$y[3]
  }
  
  # D_values:
  d_vals <- c()
  for(i in 1:length(oa)){
    d_vals[i] <- oa[[i]]$y[4]
  }
  
  # NP_ratio:
  np_ratio <- c()
  for(i in 1:length(oa)){
    np_ratio[i] <- oa[[i]]$y[3] / oa[[i]]$y[2]
  }
  
  # PPT:
  ppt <- c()
  for(i in 1:length(oa)){
    ppt[i] <- oa[[i]]$x[4] - oa[[i]]$x[2]
  }
  
  # Maximum amplitude:
  max_amp <- c()
  for(i in 1:length(oa)){
    inx <- solve(pw[[i]], b = 0, deriv = 1)
    iny <- predict(pw[[i]], inx)
    max_amp[i] <- max(iny)
  }
  
  # Total AUC:
  auc <- c()
  for(i in 1:length(oa)){
    wave <- p[, (i+1)]
    v <- which(!is.na(wave))
    auc[i] <- AUC(x = v, y = wave[v], method = "spline")
  }
  
  # AUC after peak of systole:
  auc_s <- c()
  for(i in 1:length(oa)){
    wave <- p[, (i+1)]
    wave <- wave[!is.na(wave)]
    s <- oa[[i]]$x[2]
    v <- 1:length(wave)
    v <- which(v > s)
    auc_s[i] <- AUC(x = v, y = wave[v], method = "spline")
  }
  
  # Length:
  l <- c()
  for(i in 1:length(oa)){
    l[i] <- length(p[, (i+1)][!is.na(p[, (i+1)])])
  }
  
  # Inflexion point area ratio (For canonical waveform use x[3], if using inflection point then use x[4]):
  ipa_ratio <- c()
  for(i in 1:length(oa)){
    wave <- p[, (i+1)][!is.na(p[, (i+1)])]
    ip <- oa[[i]]$x[3]
    v <- 1:length(wave)
    v <- which(v < ip)
    s_auc <- AUC(x = v, y = wave[v], method = "spline")
    v <- 1:length(wave)
    v <- which(v > ip)
    d_auc <- AUC(x = v, y = wave[v], method = "spline")
    ipa_ratio[i] <- s_auc/d_auc
  }
  
  # Peak to Notch time (relative to length of wave):
  pn_time <- c()
  for(i in 1:length(oa)){
    pn_time[i] <- (oa[[i]]$x[3] - oa[[i]]$x[2])/l[i]
  }
  
  # Notch-time ratio = time interval from notch to end of p / time interval from notch to beginning of p:
  nt_ratio <- c()
  for(i in 1:length(oa)){
    wave <- p[, (i+1)][!is.na(p[, (i+1)])]
    # Find next_o i.e the last value of the wave:
    next_o <- max(which(is.na(wave) ==  0)) 
    nt_ratio[i] <- (next_o - oa[[i]]$x[3]) / (oa[[i]]$x[3] -  oa[[i]]$x[1])
  }
  
  # Reflectance peak to forward peak ratio (augmentation index):
  ai <- c()
  for(i in 1:length(oa)){
    ai[i] <- oa[[i]]$y[4] / oa[[i]]$y[2]
  }
  
  # Alternative augmentation index:
  aai <- c()
  for(i in 1:length(oa)){
    aai[i] <- (oa[[i]]$y[2] - oa[[i]]$y[4]) / oa[[i]]$y[2]
  }
  
  # Crest time (O to S time):
  ct <- c()
  for(i in 1:length(oa)){
    ct[i] <- oa[[i]]$x[2] - oa[[i]]$x[1]
  }
  
  features <- data.frame(s_vals, n_vals, d_vals, np_ratio, ppt, max_amp, auc, auc_s, l, ipa_ratio, pn_time, nt_ratio, ai, aai, ct)
  return(features)
}



# References:
# Elgendi et al, 2018: Toward generating more diagnostic features from photoplethysmogram waveforms



# HED functions:

# 1. GetParticipants (ISO dataset specific)
# 2. GetDirec (ISO dataset specific)
# 3. GetPairs (ISO dataset specific)
# 4. Undetrend (ISO dataset specific)
# 5. FactorAdjust (ISO dataset specific)
# 6. OffsetAdjust (ISO dataset specific)
# 7. FindUndetrendingParams (ISO dataset specific)
# 8. AddOutput
# 9. FindStartParams
# 10. FindWithinParams
# 11. make_matrix
# 12. extractOutput
# 13. FixOutput
# 14. UpdateBeat
# 15. FixBaseline
# 16. PlotFits
# 17. osnd_fit
# 18. ArrangeOutputs
# 19. model2.GetSegment
# 20. model2.Excess
# 21. model2.Peak
# 22. model2.SubtractExcessPeak
# 23. model2.ChiSq3
# 24. model2.ChiSq4
# 25. model2.Rebuild2
# 26. model2.Excess.Inv2
# 27. model2.FIX_PAR3
# 28. model2.FixParams3
# 29. simplex.MakeSimplex2
# 30. simplex.MakeSimplex3
# 31. simplex.Run2
# 32. simplex.HypoCentre
# 33. simplex.SortHighLow
# 34. PlotRejects
# 35. PlotWavesCarriedForward



GetParticipants <- function(direc){
  string_list <- list.files(path = direc)
  string_list <- string_list[-119]
  string_list <- substr(string_list, 0, 5)
  rejected_ts <- c("AL826", "AO139", "AP493", "AU602", "AZ985", "AZ883")
  to_remove <- c()
  for(i in 1:length(rejected_ts)){
    tmp <- which(string_list == rejected_ts[i])
    if(length(tmp) > 0){to_remove <- c(to_remove, tmp)}
  }
  string_list <- string_list[-to_remove]
  return(string_list)
}


GetDirec <- function(run, Participants, dir){
  subjectID <- Participants[run]
  participant_number <- run
  direc <- paste(dir, subjectID, sep = "")
  scan_no <- list.files(direc)[grep(list.files(direc), pattern = "scan")]
  # If it's an 'ISO_ONLY' file another step is needed:
  if(length(scan_no) < 1){
    a <- list.files(substr(direc, 1, nchar(direc)-5))
    b <- grep(a, pattern = subjectID)
    direc <- paste(dir, a[b], sep = "")
    scan_no <- list.files(direc)[grep(list.files(direc), pattern = "scan")]
  }
  direc <- paste(direc, "/", scan_no, sep = "")
  physio_file <- list.files(direc)[grep(list.files(direc), pattern = "physio")]
  direc <- paste(direc, "/", physio_file, sep = "")
  temp <- c(direc, subjectID)
  return(temp)
}



GetPairs <- function(direc, run_order, participant_number, subjectID){
  
  str <- direc
  direc <- list.files(path = direc)
  
  # Refine direc:
  direc <- direc[-grep("REST", direc)]                                 # removing irrelevant files... 
  direc <- direc[-grep(".fig", direc)]  
  direc <- direc[-grep(".ecgpk", direc)]  
  direc <- direc[-grep(".resppk", direc)]  
  direc <- direc[-grep(".hrrv", direc)]  
  direc <- direc[grep("ECG", direc)]                                  
  
  dose_numbers <- parse_number(direc)                                  # extracting numbers from files
  direc <- direc[rev(order(dose_numbers, decreasing = T))]             # Arrange files in numerical order
  if(length(direc) < 6){message("Warning: Some time series are missing from this folder")}
  # Extract dose order:
  dose_order <- run_order[which(run_order$subj_ID == subjectID), 3]    
  dose_order <- as.numeric(unlist(strsplit(dose_order,",")))
  
  # Arrange files in order of escalating dose:
  direc <- direc[rev(order(dose_order, decreasing = T))]
  
  # Randomized pairing (probability 0.5):                              # Setting the same seed and having individual participant numbers will mean randomization is reproducible and varied across participants
  set.seed(32)
  
  if(rnorm(participant_number)[participant_number] > 0){
    # Pair 1:5 and 2:6
    pair1 <- c(direc[5], direc[1])
    pair2 <- c(direc[6], direc[2])
  }else{
    # Pair 1:6 and 2:5
    pair1 <- c(direc[6], direc[1])
    pair2 <- c(direc[5], direc[2])
  }
  
  pairs <- list(pair1, pair2)
  
  # Return pairs:
  return(pairs)
}



UnDetrend <- function(ppg,factor=0,offset=1)    
{
  k <- offset * (1-factor)
  n <- nrow(ppg)
  result <- (1:n)*0
  result[1] = ppg[1,2]
  
  for (i in 2:n)
  {
    result[i] = ppg[i,2] - ppg[i-1,2] * factor - k + result[i-1]
  }
  
  return(result)
}


FactorAdjust <- function(data, factorCutoff, ppg, u, beat, a., test, gs=gs, beatTime, nextTime, plot = T){
  
  # Calculate the Gradient of the tail of the beat:
  tail <- c(data[nrow(data), 2], data[nrow(data)-1, 2], data[nrow(data)-2, 2], data[nrow(data)-3, 2], 
            data[nrow(data)-4, 2])  
  
  # This will determine whether the gradient is initially found to be positive or not, which will influence the events that follow. 
  
  # So make initial assessment based on last 5 values. 
  # If the gradient of last 5 values is positive, use values from before the last 5. 
  # If the last 5 were positive due to noise, this will make sure it is counted as negative from that point on. 
  # If the last 5 were genuinely positive, then the positivity should extend beyond the last 5
  # If the gradient of last 5 values is negative, use the last values. 
  # Don't tend to see falsely negative values, only falsely positive. 
  # Choosing values earlier than the last 5 
  
  # There are actually two criteria I'm now using to iterate the while loop. The first is the gradient, the second 
  # is a check to make sure the minimum of the data segment is not in the middle 50% i.e due to the notch. 
  
  # This means that the second gradient check and subsequent ifelse statements are necessary, because positive and negative gradients
  # could enter the while loop if only the low notch criteria is met. 
  
  # As far as I can tell, there needs to be an initial check of the last 5 data segments to determine which section of the tail will then
  # be iterated on in the while loop. 
  
  # The reasoning for this is as follows:
  
  #  Basically, at the end of waves you can have 4 shapes:
  # 1. prolonged negative slope (morphologically normal)
  # 2. prolonged positive slope (due to heavy detrending)
  # 3. ~5 positive slope coming out of the notch, followed by ~5 negative slope as a rather short tail i.e an n shape (these tend to be ISO waves as they are short due to high HR) 
  # 4. ~5 positive slope at the end of the tail, preceded by a prolonged negative slope i.e a v shape (due to tail noise)
  
  # With an initial assessment determining the iterations that follow:
  # 1. and 2. will be processed the same regardless of which of the last 10 datapoints the gradient is taken from
  # 3. In the initial 'last 5' datapoint assessment, these will be found to be negative and thus will continue being iterated on based on the last 5 datapoints (which is what we want)
  # 4. In the initial 'last 5' datapoint assessement, these will be found to be positive, and so will be iterated on based on NOT the last 5 datapoints (also what we want)
  
  tail <- rev(tail)
  xx. <- 1:length(tail)
  y. <- lm(tail~xx.)
  # Adjust the factor value until the gradient of the tail reaches an appropriate threshold, and the min value is not the notch:
  factor_value <- 1  
  while(y.[[1]][2] > factorCutoff | (which.min(data[, 2]) > quantile(1:nrow(data))[[2]] & which.min(data[, 2]) < quantile(1:nrow(data))[[4]]) ){    
    if(factor_value < 0.7){break}
    factor_value <- factor_value - 0.01
    ppg2 <- ppg
    ppg2[, 2] <- u(ppg,factor=factor_value,offset=1)
    #beatTime <- beat[test + a., 1]   # same adjustments made here as first two lines
    #nextTime <- beat[test + (a. + 1), 1] 
    seg <- c(which(ppg$`time (s)` ==  beatTime), 0, which(ppg$`time (s)` == nextTime))
    data <- gs(ppg2,seg)
    if(plot == TRUE){plot(data, pch = 19)}
    # If the gradient is positive, ignore the last 5 values... 
    if(y.[[1]][2] > 0){
      tail <- c(data[nrow(data)-5, 2], data[nrow(data)-6, 2], data[nrow(data)-7, 2], data[nrow(data)-8, 2], 
                data[nrow(data)-9, 2], data[nrow(data)-10, 2], data[nrow(data)-11, 2], data[nrow(data)-12, 2],
                data[nrow(data)-13, 2], data[nrow(data)-14, 2])
    }else{
      tail <- c(data[nrow(data), 2], data[nrow(data)-1, 2], data[nrow(data)-2, 2], data[nrow(data)-3, 2], 
                data[nrow(data)-4, 2], data[nrow(data)-5, 2], data[nrow(data)-6, 2], data[nrow(data)-7, 2], 
                data[nrow(data)-8, 2], data[nrow(data)-9, 2])
    }
    tail <- rev(tail)
    xx. <- 1:length(tail)
    y. <- lm(tail~xx.)
  }
  
  return(factor_value)
}


OffsetAdjust <- function(ppg3, ppg, u = UnDetrend, factor_value, plot = F){
  if(factor_value == 1){   # if no changes in factor value were needed, no need to correct offset
    print("no adjustment needed")
    return(1)
  } 
  vv. <- ppg3[, 1]      
  yv. <- lm(ppg3[, 2]~vv.)
  offset_value <- 1
  while(yv.[[1]][2] > 0){
    if(yv.[[1]][2] > 5){
      offset_value <- offset_value + 1 
    }else{
      if(yv.[[1]][2] > 1){
        offset_value <- offset_value + 0.5
      }else{
        offset_value <- offset_value + 0.05
      }
    }
    ppg3 <- data.frame(ppg[,1], u(ppg,factor=factor_value,offset=offset_value))
    vv. <- ppg3[, 1]      
    yv. <- lm(ppg3[, 2]~vv.)
    if(plot == TRUE){
      if(yv.[[1]][2]>0){plot(ppg[,1],u(ppg,factor=factor_value,offset=offset_value), type = "l")}
      if(yv.[[1]][2]>0){abline(a = yv.[[1]][1], b = yv.[[1]][2], col = "red")}
      if(yv.[[1]][2]>0){print(yv.[[1]][2])} 
    }
  }
  return(offset_value)
}


FindUndetrendingParams <- function(direc, gs = model2.GetSegment, u = UnDetrend, oa = OffsetAdjust, fa = FactorAdjust, factorCutoff = 0, sr = samplingRate, pk_thrshd, pairs, plot = T){
  
  if(plot == TRUE){p <- TRUE}else{p <- FALSE}
  
  # Factor Value Adjustment (4 waves from each time series):
  
  factor_value_vec <- c()
  for(ps in 1:2){
    if(ps == 1){pair <- pairs[[1]]}else{pair <- pairs[[2]]}
    for(pr in 1:2){
      new_direc <- paste(direc, "/", pair[pr], sep = "")
      ppg <- read.csv(new_direc, sep = "")   
      ppg <- data.frame(
        time = (0:(nrow(ppg)-1)) / samplingRate,
        ppg = ppg[,1]
      )
      names(ppg)[1] <- "time (s)"
      names(ppg)[2] <- "Detrended"
      
      # Find beats:    
      n <- dim(ppg)[1]
      vpg <- ppg[2:n,2] - ppg[1:(n-1),2]
      beat <- data.frame(ppg[which(vpg[1:(n-1)] < pk_thrshd & vpg[2:n] >= pk_thrshd),1])  
      if(nrow(beat) < length(vpg)/(sr*2)){    # if number of beats found suggests a HR of < 30bpm, trigger warning
        message("Warning: Minimal peaks found - consider resetting vpg peak threshold")
        Sys.sleep(10)
      }
      
      t_value <- c()
      for(i in 1:4){
        print(i)
        # We don't know the exact time of onset of iso at this point, but we can infer the IBI from beat and identify
        # waves around the minimum point:
        
        # Find rolling median:
        pre_ibi <- abs(beat[1:(nrow(beat)-1), 1] - beat[2:nrow(beat), 1])
        meds <- rollmedian(pre_ibi, k = 19)
        # plot(meds)
        
        # Failsafes in case the minimum is close to the end / beginning
        min <- which.min(meds)
        if((min - 10) < 1){min = 11}
        if((min + 10) > length(meds)){min = length(meds) - 11}
        test <- round(quantile((min-10):(min+30)))[[i]]
        
        # Extract the relevant beat:
        beatTime <- beat[test, 1]   
        nextTime <- beat[(test + 1), 1]
        seg <- c(which(ppg$`time (s)` ==  beatTime), 0, which(ppg$`time (s)` == nextTime))
        data <- gs(ppg,seg)
        
        # Since we are using a less robust method to find beats, there is a chance of multi-beat segments:
        # Detect them and choose the next segment if so...:
        a. <- 0
        while (nrow(data) > sr*1.5 | nrow(data) < (0.375*sr) | 
               sum(c(order(data[, 2], decreasing = T)[1:5][-1], order(data[, 2], decreasing = T)[1:5][1]) > 20) > 0 ){    # this line insures that the max (5) points of the data are not far apart in time (as would be the case if two peaks were present)
          beatTime <- beat[test + a., 1]  
          nextTime <- beat[test + (a. + 1), 1]
          seg <- c(which(ppg$`time (s)` ==  beatTime), 0, which(ppg$`time (s)` == nextTime))
          data <- gs(ppg,seg)
          a. <- a. + 1
        }
        if(a. != 0){
          a. <- a. - 1
        }
        
        if(plot == TRUE){plot(data, pch = 19)}
        t_value[i] <- fa(data, factorCutoff, ppg, u, beat, a., test, gs, beatTime, nextTime, plot = p)
      }
      
      factor_value_vec <- c(factor_value_vec, t_value) 
    }
  }
  
  # Find the 2nd from minimum factor value (and check which dose level time series it comes from):
  min2 <- order(factor_value_vec)[2]  # ascending order
  if(sum(min2 == 5:8) > 0 | sum(min2 == 13:16) > 0){
    message("Warning: factor value taken from 0mg time series")
    Sys.sleep(8)
    min1 <- order(factor_value_vec)[1]
    if(sum(min1 == 5:8) > 0 | sum(min1 == 13:16) > 0){
      message("though minimum factor value found in 2mg time series")
      Sys.sleep(8)
    }
  }
  factor_value <- factor_value_vec[min2]
  
  # Offset Value Adjustment:
  
  offset_value <- c()
  for(ps in 1:2){
    if(ps == 1){pair <- pairs[[1]]}else{pair <- pairs[[2]]}
    for(pr in 1:2){
      new_direc <- paste(direc, "/", pair[pr], sep = "")
      # Load time series:
      ppg <- read.csv(new_direc, sep = "")   
      ppg <- data.frame(
        time = (0:(nrow(ppg)-1)) / samplingRate,
        ppg = ppg[,1]
      )
      names(ppg)[1] <- "time (s)"
      names(ppg)[2] <- "Detrended"
      # Adjust factor:  
      ppg3 <- data.frame(ppg[,1],UnDetrend(ppg,factor=factor_value,offset=1))
      # Adjust offset:
      offset_value <- c(offset_value, oa(ppg3, ppg, u = UnDetrend, factor_value, plot = p))
    }
  }
  
  offset_value <- median(offset_value)
  values <- c(factor_value, offset_value)
  return(values)
}



AddOutput <- function(beat){
  beat$First      = 1:nrow(beat) * 0
  beat$Last       = 1:nrow(beat) * 0
  beat$Baseline   = 1:nrow(beat) * 0
  beat$Baseline2  = 1:nrow(beat) * 0
  beat$STime      = 1:nrow(beat) * 0
  beat$SAmplitude = 1:nrow(beat) * 0
  beat$SWidth     = 1:nrow(beat) * 0
  beat$DTime      = 1:nrow(beat) * 0
  beat$DAmplitude = 1:nrow(beat) * 0
  beat$DWidth     = 1:nrow(beat) * 0
  beat$NTime      = 1:nrow(beat) * 0
  beat$NAmplitude = 1:nrow(beat) * 0
  beat$NWidth     = 1:nrow(beat) * 0
  beat$config.rate= rep(0.75, nrow(beat))  
  return(beat)
}


FindStartParams <- function(batch_number, beats_in, beat, ppg, gs = model2.GetSegment, e = model2.Excess, sep = model2.SubtractExcessPeak, o_points = inflexX[o_orig], wuv = wuv, inflexX = inflexX, all_beats = FALSE, plot = FALSE){
  nBeats <- nrow(beat)
  seg <- c(0,0,0)
  if((batch_number*beats_in) > nBeats){
    print("Batch and beat values request more beats than are in time series, defaulting to max number of beats")
    maxn <- nBeats 
  }else{
    maxn <-(batch_number*beats_in)
    if(all_beats == TRUE){maxn <- maxn + (nrow(beat) - maxn)}
  }
  
  # You can check 02 and 0 points imported with the following plot (useful for debugging):
  # plot(5900:6500, ppg$Detrended[5900:6500], type = "l")
  # points(o_points, rep(0, length(o_points)))
  # points(inflexX[wuv$o2], rep(0, length(wuv$o2)), col = "red")
  
  for (i in 1:maxn){  
    
    # Find Beat:
    beatTime <- beat[i,1]
    current_o <- which(ppg[, 1] == beatTime)
    # Find the minimum o_point that is after the current o point
    next_o <- min(which(o_points > current_o))
    next_o <- o_points[next_o]
    # Make sure that through rounding you haven't chosen the same o twice
    if ((next_o - current_o) < 5){
      next_o <- min(which(o_points > current_o)) + 1
      next_o <- o_points[next_o]
    }
    # Find the closest ppg time value to the next_o (perhaps a rounding issue was causing problems before...?)
    nextTime <- ppg[round(next_o), 1]
    seg <- c(which(ppg$`time (s)` ==  beatTime), 0, which(abs(ppg[, 1] - nextTime) == min(abs(ppg[, 1] - nextTime))))
    data <- gs(ppg,seg)
    if(plot == TRUE){plot(data)}
    
    tStart <- ppg[seg[1],1]
    yPrev <- ppg[max(seg[1]-1,1),2]
    
    amp <- max(data[, 2]) - min(data[, 2])
    constant <- 0.1092254*amp
    baseline <- min(data[,2]) - constant    
    residue <- e(data[,2], ppg[seg[1]-1,2], -0) ## - 0?     
    
    count <- nrow(data)
    excess <- 1:count * 0.0
    excess[1] = data[1,2] - (baseline + config.rate*(yPrev-baseline))
    for (j in 2:count){
      excess[j] = data[j,2] - (baseline + config.rate*(data[j-1,2]-baseline))
    }
    rm(count)
    rm(j)
    # plot(data[,1],excess, type = "l")
    par <- 1:10 * 0.
    par[1] = baseline
    residue <- excess  
    
    # S peak
    peak.w <- which(data[,1] > beat[i,1]-0.2 & data[,1] < beat[i,1]+0.2)  
    peak.t <- data[peak.w,1]
    peak.y <- residue[peak.w]
    # plot(data[,1],excess)   
    # lines(peak.t,peak.y)
    par[3] <- max(peak.y)      # height of the S peak
    par[2] <- peak.t[which(peak.y==par[3])]   # timing of the S peak
    par[4] <- 0.25    # width of the S-peak, a priori
    rm(peak.w,peak.t,peak.y)
    residue <- sep(data[,1],residue,par[2:4])
    # plot(data[,1],excess)
    # lines(data[,1],residue)  
    
    # D peak
    peak.w <- which(data[,1] > beat[i,1]+0.3 & data[,1] < beat[i,1]+0.6)  # this finds a range in where to find the peak of d
    peak.t <- data[peak.w,1]     # this finds the time corresponding to the peak
    peak.y <- residue[peak.w]
    # plot(data[,1],excess)   
    # lines(peak.t,peak.y)
    par[6] <- max(peak.y)     
    par[5] <- peak.t[which(peak.y==par[6])]   
    par[7] <- 0.25  
    rm(peak.w,peak.t,peak.y)
    residue <- sep(data[,1],residue,par[5:7])
    # plot(data[,1],excess)
    # lines(data[,1],residue)
    
    # N peak
    t <- par[2] + c(0.25,0.75) * (par[5]-par[2])
    peak.w <- which(data[,1] > t[1] & data[,1] < t[2])
    peak.t <- data[peak.w,1]
    peak.y <- residue[peak.w]
    # plot(data[,1],excess)   
    # lines(peak.t,peak.y)
    par[9] <- max(peak.y)
    par[8] <- peak.t[which(peak.y==par[9])]
    par[10] <- 0.25
    rm(peak.w,peak.t,peak.y,t)
    residue <- sep(data[,1],residue,par[8:10])
    # plot(data[,1],excess)
    # lines(data[,1],residue)
    
    # Store parameters
    w <- seg[1]:seg[3]
    ppg$Baseline[w] <- baseline
    ppg$Excess[w] <- excess
    ppg$Residue[w] <- residue
    rm(w,excess,residue,data,baseline,yPrev,nextTime,tStart)
    beat[i,3:4]  = c(seg[1],seg[3])
    beat[i, 5] <- par[1]
    beat[i,6:15] = par
    beat[i,10] = beat[i,10]-beat[i,7]
    beat[i,13] = beat[i,13]-beat[i,7]
    rm(par)
  }
  rm(seg)
  temp <- list(beat, ppg)
  return(temp)
}


FindWithinParams <- function(beats_in, ppg, beat, gs = model2.GetSegment, fp = model2.FixParams3, ms = simplex.MakeSimplex3, m2 = model2.ChiSq3, beat_vector = beat_vector, renal_param = renal_param, dias_param = dias_param, sys_time, w){
  a <- list()
  for(i in 1:beats_in){         
    
    par <- as.numeric(beat[i,5:16])
    #par <- fp(data[, 1:2], par, rp = renal_param, sys_t = sys_t)      # Do we need to fix the parameters here...?
    beat_indi <- list(1, beat_vector[[2]][i], beat_vector[[3]][i])
    
    a[[i]] <- ms(ppg = ppg, param = par, f = m2, inScale = 0.1, inTol=-1, beat_vector = beat_indi, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time[i], w = w[i]) 
  }
  return(a)
}

make_matrix <- function(sim, a){
  # Save the top row of sim for replication:
  top_row_sim <- sim[1, c(5:6, 8:9, 11:12)]              
  #Remove redundant rows and columns from sim:
  sim <- sim[c(6:7, 9:10, 12:13), c(5:6, 8:9, 11:12)]                           
  # You need the top_row of sim to replicate:
  top_row_sim <- matrix(data = top_row_sim, nrow = 6, ncol = 6, byrow = TRUE)    
  
  # Make the a values just the rows where within beat parameters are changed  
  # Save the top row of each matrix of a for replication...
  top_row <- list()
  for(i in 1:beats_in){                                                          
    top_row[[i]] <- a[[i]][1, -c(5:6, 8:9, 11:12)]
    a[[i]] <- a[[i]][-c(1, 6:7, 9:10, 12:13), -c(5:6, 8:9, 11:12)]             
  }
  # Top row needs to be replicated for each within beat row when the across-beat params are being changed:
  for(i in 1:beats_in){
    top_row[[i]] <- matrix(data = top_row[[i]], ncol = 6, nrow = 6, byrow = TRUE)  
  }
  
  
  # Assemble Matrix:
  
  # Bind replicate rows of top_row for each beat to sim:              
  for(i in 1:beats_in){
    sim <- cbind(sim, top_row[[i]])
  }
  
  
  # Create rows for each beat
  beat_rows <- list()                                                 
  for(i in 1:beats_in){                                             
    
    # Add values to the left
    beat_rows[[i]] <- top_row_sim   
    if(i != 1){
      for(j in 1:(i-1)){
        beat_rows[[i]] <- cbind(beat_rows[[i]], top_row[[j]])  
      }
    }
    
    # Add A
    beat_rows[[i]] <- cbind( beat_rows[[i]], a[[i]])
    
    # Add values to the right
    if(i == beats_in){
      break
    }else{
      for(j in (i+1):beats_in){
        beat_rows[[i]] <- cbind(beat_rows[[i]], top_row[[j]])
      }
    }
  }
  
  # Bind rows together:
  for(i in 1:beats_in){
    sim <- rbind(sim, beat_rows[[i]])
  }
  
  # Add the top row:
  final_top_row <- top_row_sim[1, ]
  for(i in 1:beats_in){
    final_top_row <- c(final_top_row, top_row[[i]][1, ])
  }
  sim <- rbind(final_top_row, sim)
  
  return(sim)
}


extractOutput <- function(beats_in, sim){
  across <- sim[1, ][1:6]
  within <- list()
  for(i in 1:beats_in){
    temp <- rep(0, 12)
    temp[c(1:4, 7, 10)] <-  sim[1, ][((i*6)+1):((i*6)+6)]  
    within[[i]] <- temp
  }
  temp <- list(across, within)
  return(temp)
}


FixOutput <- function(beats_in, beat, ppg, gs = model2.GetSegment, fp = model2.FixParams3, across = output[1], within = output[2], sys_time = sys_time){
  fixed <- list()
  for(i in 1:beats_in){
    seg <- c(beat[i,3],0,beat[i,4])
    data <- model2.GetSegment(ppg,seg)
    rm(seg)
    fixed[[i]] <- model2.FixParams3(data, params = as.numeric(within[[i]]), across_beat_params = across, sys_t = sys_time[i])
  } 
  return(fixed)
}


UpdateBeat <- function(beats_in, beat, fixed){
  new_beat <- data.frame(matrix(0, ncol = 12, nrow = beats_in))
  for(i in 1:beats_in){
    new_beat[i, ] <- fixed[[i]]
  }
  new_beat <- cbind(beat[, 1:4], new_beat)
  return(new_beat)
}


FixBaseline <- function(new_beat, f = model2.ChiSq3, renal_param, dias_param, sys_time, w){
  for(j in 1:nrow(new_beat)){
    if(abs(new_beat[j, 6] - new_beat[j, 5]) < 5){
      # Assess fit:
      wave_check <- model2.ChiSq3(data = ppg, params = as.numeric(new_beat[j, 5:16]), beats = list(1, new_beat[j, 3], new_beat[j, 4]), beat = NULL, a = NULL, plot = FALSE, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time[j], w = w[j])
      # Assess fit with baselines equal:
      wave_check2 <- model2.ChiSq3(data = ppg, params = c(rep(new_beat[j, 5], 2), as.numeric(new_beat[j, 7:16])), beats = list(1, new_beat[j, 3], new_beat[j, 4]), beat = NULL, a = NULL, plot = FALSE, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time[j], w = w[j])
      # If baselines equal gives a better value of ChiSq, fix them to be so:
      if(wave_check2 < wave_check){
        new_beat[j, 6] <- new_beat[j, 5]
      }
    }
  }
  return(new_beat)
}


PlotFits <- function(beats_in, ppg, beat2, gs = model2.GetSegment, rb = model2.Rebuild2){
  for(i in 1:beats_in){
    seg <- c(beat[i,3],0,beat[i,4])  
    data <- model2.GetSegment(ppg,seg)
    yPrev <- ppg[seg[1]-1,2]
    xPrev <- ppg[seg[1]-1, 1]
    xNext <- ppg[seg[3], 1]
    rm(seg)
    temp<-model2.Rebuild2(data, yPrev, as.double(beat2[i,]),TRUE)    
    plot(data[, 1], data[, 2], ylim = c(beat2$Baseline[1]*1.5, max(data[, 2]*1.2)), main = paste(c("batch", k, "wave", i), collapse = " "))   # ylim = c(76, 86)
    lines(data[,1],temp)
    # Plot baselines:
    lines(c(xPrev, (beat2[i, 3]  + (1*beat2[i, 6]))), rep(beat2[i, 1], 2))   
    lines(c((beat2[i, 3]  + (1*beat2[i, 6])), xNext), rep(beat2[i, 2], 2))
    # Plot systolic:
    # Always need 12 elements, but set amplitude and width to 0 for the peaks that aren't being used. 
    par <- as.double(beat2[i,])
    par[c(7:8, 10:11)] <- 0
    temp<-model2.Rebuild2(data,yPrev,par,TRUE)
    lines(data[,1],temp, col = "red")
    # Plot diastolic:
    par <- as.double(beat2[i,])
    par[c(4:5, 10:11)] <- 0
    temp<-model2.Rebuild2(data,yPrev,par,TRUE)      
    lines(data[,1],temp, col = "blue")
    # Plot renal:
    par <- as.double(beat2[i,])
    par[c(4:5, 7:8)] <- 0
    temp<-model2.Rebuild2(data,yPrev,par,TRUE)
    lines(data[,1],temp, col = "green")
  }
}

GGplotFits <- function(beats_in, ppg, beat2, gs = model2.GetSegment, rb = model2.Rebuild2, run, pr, p = F){
  
  if(pr == 1){pr = 0}
  
  for(i in 1:beats_in){
    seg <- c(beat[i,3],0,beat[i,4])  
    data <- model2.GetSegment(ppg,seg)
    yPrev <- ppg[seg[1]-1,2]
    xPrev <- ppg[seg[1]-1, 1]
    xNext <- ppg[seg[3], 1]
    rm(seg)
    temp<-model2.Rebuild2(data, yPrev, as.double(beat2[i,]),TRUE)    
    fit <- model2.Rebuild2(data, yPrev, as.double(beat2[i,]),TRUE) 
    
    # Create waves dataframe:
    waves <- data.frame(data)
    waves <- cbind(data, fit)
    par <- as.double(beat2[i,])
    par[c(7:8, 10:11)] <- 0
    temp<-model2.Rebuild2(data,yPrev,par,TRUE)
    waves <- cbind(waves, temp)
    par <- as.double(beat2[i,])
    par[c(4:5, 7:8)] <- 0
    temp<-model2.Rebuild2(data,yPrev,par,TRUE)
    waves <- cbind(waves, temp)
    par <- as.double(beat2[i,])
    par[c(4:5, 10:11)] <- 0
    temp<-model2.Rebuild2(data,yPrev,par,TRUE)      
    waves <- cbind(waves, temp)
    
    # Baselines can be added as line segments with geom_segment:
    b1y <- beat2[i, 1]
    b1x <-  c(xPrev, (beat2[i, 3]  + (1*beat2[i, 6])))  
    b2x <- c((beat2[i, 3]  + (1*beat2[i, 6])), xNext)
    b2y <-  beat2[i, 2]
    
    # Make component waves + fit into splines:
    sfunction <- splinefun(1:nrow(waves), fit, method = "natural")
    fit2 <- sfunction(seq(1, nrow(waves), 0.1), deriv = 0) 
    par <- as.double(beat2[i,]) 
    par[c(7:8, 10:11)] <- 0
    sys<-model2.Rebuild2(data,yPrev,par,TRUE)
    sfunction <- splinefun(1:nrow(waves), sys, method = "natural")
    sys2 <- sfunction(seq(1, nrow(waves), 0.1), deriv = 0) 
    par <- as.double(beat2[i,])
    par[c(4:5, 10:11)] <- 0
    dias<-model2.Rebuild2(data,yPrev,par,TRUE)      
    sfunction <- splinefun(1:nrow(waves), dias, method = "natural")
    dias2 <- sfunction(seq(1, nrow(waves), 0.1), deriv = 0) 
    par <- as.double(beat2[i,])
    par[c(4:5, 7:8)] <- 0
    R1<-model2.Rebuild2(data,yPrev,par,TRUE)
    sfunction <- splinefun(1:nrow(waves), R1, method = "natural")
    R12 <- sfunction(seq(1, nrow(waves), 0.1), deriv = 0) 
    
    # Build waves_stacked 
    time <- waves[, 1]
    time <- seq(from = time[1], to = time[length(time)], length.out = length(fit2))
    fit2 <- data.frame(time, fit2)
    fit2 <- cbind(fit2, rep("fit", nrow(fit2)))
    sys2 <- data.frame(time, sys2)
    sys2 <- cbind(sys2, rep("systolic wave", nrow(sys2)))
    dias2 <- data.frame(time, dias2)
    dias2 <- cbind(dias2, rep("2nd reflectance wave", nrow(dias2)))
    R12 <- data.frame(time, R12)
    R12 <- cbind(R12, rep("1st reflectance wave", nrow(R12)))
    colnames(fit2) <- c("x", "values", "Wave")
    colnames(sys2) <-  c("x", "values", "Wave")
    colnames(dias2) <-  c("x", "values", "Wave")
    colnames(R12) <-  c("x", "values", "Wave")
    waves_stacked_final <- rbind(fit2, sys2, dias2, R12)
    
    # stack data:
    data_stacked <- data.frame(data)
    data_stacked <- cbind(data_stacked, rep("data"))
    colnames(data_stacked) <- c("x", "values", "Wave")
    
    # Plot! 
    if(p == T){
      c <- ggplot(data = waves_stacked_final, aes(x = x, y = values, col = Wave)) + geom_line(aes(size = Wave, alpha = Wave)) +
        scale_color_manual(values = c("#03fc7b", "#03b5fc", "black", "black", "black", "#ff4242", "black")) + scale_size_manual(values = c(0.7, 0.7, 1.5, 0.7, 0.7)) + 
        scale_alpha_manual(values = c(1, 1, 1, 1, 1)) + ylab("PPG Signal") + xlab("Time") + geom_point(data = data_stacked) + 
        geom_segment(aes(x = data[1, 1], y = b1y, xend = b1x[2], yend = b1y, colour = "black")) + geom_segment(aes(x = b2x[1], y = b2y, xend = b2x[2], yend = b2y, colour = "black")) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
        theme(legend.position = "none") + ggtitle(paste0("Participant", " ", run, "\n", pr, " ", "mcg", sep = "")) # + ylim(-350, 1600) 
    }else{
      ggplot(data = waves_stacked_final, aes(x = x, y = values, col = Wave)) + geom_line(aes(size = Wave, alpha = Wave)) +
        scale_color_manual(values = c("green", "blue", "black", "black", "black", "red", "black")) + scale_size_manual(values = c(0.7, 0.7, 1.5, 0.7, 0.7)) + 
        scale_alpha_manual(values = c(1, 1, 1, 1, 1)) + ylab("PPG Signal") + xlab("Time") + geom_point(data = data_stacked) + 
        geom_segment(aes(x = b1x[1], y = b1y, xend = b1x[2], yend = b1y, colour = "black")) + geom_segment(aes(x = b2x[1], y = b2y, xend = b2x[2], yend = b2y, colour = "black")) + 
        theme(legend.position = "none") # + ylim(-1000, 2000)
    }
  }  
}



osnd_fit <- function(bf = beat_final, ppg, gs = model2.GetSegment, r = model2.Rebuild2, sf = splinefun, dp = diast_pk, oa = osnd_of_average, sr = samplingRate, plot = FALSE){
  
  osnd_diff <- list()
  for(i in 1:nrow(bf)   ){  #nrow(bf)  
    # Find the correct data segments and corresponding model fit:
    seg <- c(bf[i,3],0,bf[i,4])  
    data <- gs(ppg,seg)
    yPrev <- ppg[seg[1]-1,2]
    xPrev <- ppg[seg[1]-1, 1]
    xNext <- ppg[seg[3], 1]
    rm(seg)
    temp <- r(data, yPrev, as.double(bf[i,-c(1:4)]),TRUE)   
    
    # Upsample fit - so that fit OSND can be calculated as the data would have been in the main script
    sfunction <- sf(1:length(temp), temp, method = "natural")
    fit <-  sfunction(seq(1, length(temp), 0.1), deriv = 0)
    # Upsample data segment:
    sfunction <- sf(1:length(data[, 2]), data[, 2], method = "natural")
    dat <-  sfunction(seq(1, length(data[, 2]), 0.1), deriv = 0)
    
    # plot(dat)
    # plot(fit)
    
    # Find OSND of fit:
    tmp <- dp(avw = fit, sr = sr, scale = T, dias_param = bf[i, 10]*sr*10)      # pass in the Dparam to help find D on fit waves
    dPeak <- tmp[1]
    xShift <- tmp[2]
    rm(tmp)
    osnd_fit <- oa(fit, dp = dPeak, diff = 0, sr = sr, plot = F)
    # Find OSND of data:
    tmp <- dp(avw = dat, sr = sr, scale = T, dias_param = bf[i, 10]*sr*10)  # why not make dias_param also informed by the model here? instead of null? try it!
    dPeak <- tmp[1]
    xShift <- tmp[2]
    rm(tmp)
    osnd_dat <- oa(dat, dp = dPeak, diff = 0, sr = sr, plot = F)
    
    # Adjust x-axis for upsampling and sampling rate, then find the difference between OSND points:
    osnd_dat$x <- osnd_dat$x/(10*sr)
    osnd_fit$x <- osnd_fit$x/(10*sr)
    osnd_diff[[i]] <- osnd_dat - osnd_fit
    
    if(plot == TRUE){
      # Plot them together:
      plot((1:length(dat))/(10*samplingRate), dat, type = "l", xlab = "time", ylab = "")
      lines((1:length(dat))/(10*samplingRate), fit, col = "red")
      points(osnd_dat, pch = 19)
      points(osnd_fit, col = "red", pch = 19)
    }
  }
  return(osnd_diff)
}


ArrangeOutputs <- function(beat_final, beat_orig, features, pulse, fit_check, ps, pr){
  
  # Rename rows
  rownames(beat_final) <- colnames(pulse)[-1][1:nrow(beat_final)]
  rownames(features) <- colnames(pulse)[-1]
  
  # Reorganize fit_check:
  # Individual wave fits:
  wave_fits <- c()
  for(i in 1:length(fit_check)){
    wave_fits <- c(wave_fits, as.numeric(fit_check[[i]][[2]]))
  }
  # Maximum error:
  max_err <- c()
  for(i in 1:length(fit_check)){
    max_err <- c(max_err, as.numeric(fit_check[[i]][[3]]))
  }
  # NRMSE:
  NRMSE <- c()
  for(i in 1:length(fit_check)){
    NRMSE <- c(NRMSE, as.numeric(fit_check[[i]][[4]]))
  }
  # aNRMSE:
  aNRMSE <- c()
  for(i in 1:length(fit_check)){
    aNRMSE <- c(aNRMSE, as.numeric(fit_check[[i]][[5]]))
  }
  fit_check <- list(wave_fits, max_err, NRMSE, aNRMSE)
  
  tmp <- list(beat_final, features, fit_check)
  return(tmp)
}


model2.GetSegment <- function(ppg,limits){
  
  w <- c(limits[1]:limits[3])
  result <- matrix(nrow=length(w),ncol=2)
  result[,1] <- ppg[[lab.time]][w]
  result[,2] <- ppg[[lab.ppg]][w]
  
  return(result)
}

model2.Excess <- function(y,offset,baseline){
  count <- length(y)
  if (length(offset) == 0){
    print("Help")
  }
  
  result <- 1:count * 0.0
  result[1] = y[1] - (baseline + config.rate*(offset-baseline))
  for (j in 2:count){
    result[j] = y[j] - (baseline + config.rate*(y[j-1]-baseline))
  }
  
  return(result)
}

model2.Peak <- function(time,peakParams){
  temp <- 2*const.pi*(time - as.double(peakParams[1]))/as.double(peakParams[3])
  temp[which(temp < -const.pi)] = -const.pi
  temp[which(temp >  const.pi)] =  const.pi
  result <- as.double(peakParams[2]) * (0.5 * (1+cos(temp)))^2
  
  return(result)
}

model2.SubtractExcessPeak <- function(time,residue,peakParams){
  result <- residue - model2.Peak(time,peakParams)
  
  return(result)
}


model2.ChiSq3 <- function(data, params, debug=FALSE, beats, optional = NULL, beat = NULL, a = NULL, plot = FALSE, renal_param, dias_param, sys_time, w){  
  
  # makesimplex3 inputs a single set of parameters (tParam) to test, whereas run.simplex2 inputs a matrix
  # This determines where within beat parameters should be extracted from, hence beat and a are NULL unless otherwise specified. 
  
  # Across-beat parameter extraction:
  if(!is.null(a)){                                            # If a 66 parameter vector has been supplied, extract the first 6 
    across_beat_params <- a[1:6]
  }else{                                                      # If not, take them from the params input
    par <- params
    across_beat_params <- par[c(5, 6, 8, 9, 11, 12)]
  }
  
  # Calculation of ChiSq for all beats:
  beat_fit <- list()
  for(i in 1:beats[[1]]){                                          # The number of beats is determined by the first object of beats
    
    # Within-beat parameter extraction:
    if(!is.null(a)){                                               # If a 66 parameter vector has been supplied, take those values
      par2 <- a[((i*6)+1):((i*6)+6)]    
    }else{                                   
      if(!is.null(beat)){                                           # If not, take the values of beat. 
        par2 <- as.numeric(beat[i, c(5:8, 11, 14)])
      }else{                                                        # If beat is also not provided, take them from the params input
        par2 <- par[c(1:4, 7, 10)]
      }                                                            
    }
    
    # Extract individual beat data:  
    seg <- c(beats[[2]][i],0,beats[[3]][i])
    dat <- model2.GetSegment(data,seg)
    rm(seg)
    
    # Extract systolic and diastolic parameters:
    sys <- par2[3]
    #dias <- across_beat_params[2] + par2[3]
    dias <- par2[3] + dias_param
    end <- which(abs(dat[, 1]-dias) == min(abs(dat[, 1] - dias)))
    
    # Find W:
    w. <- w[i]
    w. <- which(abs(dat[, 1] - w.) == min(abs(dat[, 1] - w.))) 
    
    sys_t <- sys_time[i]
    start <- which(abs(dat[, 1]-sys_t) == min(abs(dat[, 1] - sys_t))) 
    
    # Fix parameters and calculate penalty:
    temp <- model2.FIX_PAR3(time = dat[, 1], within_beat_params = par2, across_beat_params = across_beat_params, renal_param = renal_param, sys_t = sys_t)  
    penalty <- temp[1]
    fixedPar <- temp[2:length(temp)]     
    rm(temp)
    
    # Calculate fit and residue:
    fit <- model2.Rebuild2(dat,dat[1,2],params = fixedPar)    
    residue <- dat[ ,2] - fit
    
    # Add to the penalty if the residual at sys_t is high:
    if(residue[start]*residue[start] > 0){
      penalty <- penalty + residue[start]*residue[start]
    }
    
    # Weighted region is W -> D (with slope)
    residue[w.:end[1]] <-  residue[w.:end[1]]*3
    if(length(residue) > end[1]){
      tail <- (end[1]+1):length(residue)
      for(j in 1:length(tail)){
        wgt <- 3 - (0.1*j)
        if(wgt < 1){wgt <- 1}
        residue[tail[j]] <- residue[tail[j]]*wgt
      }
    }
    
    # Calculate Reduced Chi-Square for the beat:
    nData <- nrow(dat)    
    nPar <- length(par2) 
    beat_fit[[i]] <- (sum(residue*residue) / (nData-nPar)) + as.numeric(penalty)
    
    if(plot == TRUE){
      plot(dat,  ylim = c(-150, 2000))      #ylim = c(76, 86) for bioradio data, ylim = c(-150, 1600) for ISO
      lines(dat[, 1], fit)
      #lines(dat[, 1], residue + dat[1, 2])
    }
  }
  
  # Summate individual beat ChiSq values:
  temp <- c()
  for(i in 1:length(beat_fit)){
    temp[i] <- beat_fit[[i]][1]
  }
  ts_fit <- sum(temp)
  return(ts_fit)
}


model2.ChiSq4 <- function(data, params,debug=FALSE, beats, beat, a = NULL, plot = FALSE, renal_param, dias_param, sys_time, w){  
  
  # Across-beat parameter extraction:
  if(!is.null(a)){                                            # If a 66 parameter vector has been supplied, extract the first 6 
    across_beat_params <- a[1:6]
  }else{                                                      # If not, take them from the params input
    par <- params
    across_beat_params <- par[c(5, 6, 8, 9, 11, 12)]
  }
  
  # Calculation of ChiSq for all beats:
  beat_fit <- c()   
  max_error <- c()   
  NRMSE <- c()
  aNRMSE <- c()
  for(i in 1:beats[[1]]){                                          # The number of beats is determined by the first object of beats
    
    # Within-beat parameter extraction:
    if(!is.null(a)){                                               # If a 66 parameter vector has been supplied, take those values
      par2 <- a[((i*6)+1):((i*6)+6)]    
    }else{                                   
      if(!is.null(beat)){                                           # If not, take the values of beat. 
        par2 <- as.numeric(beat[i, c(5:8, 11, 14)])
      }else{                                                        # If beat is also not provided, take them from the params input
        par2 <- par[c(1:4, 7, 10)]
      }                                                            
    }
    
    # Extract individual beat data:  
    seg <- c(beats[[2]][i],0,beats[[3]][i])
    dat <- model2.GetSegment(data,seg)
    rm(seg)
    
    # Extract systolic and diastolic parameters:
    sys <- par2[3]
    #dias <- across_beat_params[2] + par2[3]
    dias <- par2[3] + dias_param
    start <- which(abs(dat[, 1]-sys) == min(abs(dat[, 1] - sys)))   
    end <- which(abs(dat[, 1]-dias) == min(abs(dat[, 1] - dias)))
    
    # Find W:
    w. <- w[i]
    w. <- which(abs(dat[, 1] - w.) == min(abs(dat[, 1] - w.))) 
    
    # Get intially estimated systolic timing and amplitude:
    sys_t <- sys_time[i]
    
    # Fix parameters and calculate penalty:
    temp <- model2.FIX_PAR3(time = dat[, 1], within_beat_params = par2, across_beat_params = across_beat_params, renal_param = renal_param, sys_t = sys_t)  
    penalty <- temp[1]
    fixedPar <- temp[2:length(temp)]     
    rm(temp)
    
    # Calculate fit, residue and max error:
    fit <- model2.Rebuild2(dat,dat[1,2],params = fixedPar)    
    residue <- dat[ ,2] - fit
    max_error[i] <- max(residue)
    
    # Before weighting the residuals, calculate NMRSE for the region we are interested in:
    if(plot == TRUE){
      plot(dat, ylim = c(-500, 2500))
      lines(dat[, 1], fit)
      lines(dat[, 1], residue, col = "Red") 
    }
    
    # Define region of interest:
    rmse_begin <- floor((1+w.)/2)    # Begin half way from O to W
    rmse_end <- end[1] + 10          # End 10 points after the d-peak (this corresponds to half way down the weighted tail)
    if(sum(is.na(residue[rmse_begin:rmse_end])) > 0){
      residue_roi <- residue[rmse_begin:length(residue)]    # If there are fewer than 10 data points after D, use as many as there are
      ind_resid <- rmse_begin:length(residue)
    }else{
      residue_roi <- residue[rmse_begin:rmse_end]
      ind_resid <- rmse_begin:rmse_end
    }
    if(plot == TRUE){lines(dat[ind_resid, 1],residue_roi, col = "green")} 
    
    # Define null model:
    fit_null <- rep(mean(dat[ind_resid, 2], trim = 0), length(residue_roi))
    if(plot == TRUE){lines(dat[ind_resid, 1], fit_null)}
    # Calculate residuals of the null model:
    residuals_of_null_model <- dat[ind_resid, 2] - fit_null
    
    # Calculate NRMSE:
    rmse_model2 <-  sqrt(mean(residue_roi^2, trim = 0))
    rmse_null <- sqrt(mean(residuals_of_null_model^2))
    NRMSE. <- 1 - (rmse_model2 / rmse_null)
    NRMSE[i] <- NRMSE.
    
    # Alternative NRMSE method (Wang et al 2013):
    # SSE / Sum of squared datapoints 
    aNRMSE[i] <- (sum(residue_roi^2) / sum(dat[ind_resid, 2]^2))*100 
    #plot(dat[ind_resid, 1], dat[ind_resid, 2]^2, type = "l")
    #lines(dat[ind_resid, 1], residue_roi^2, col = "red")
    
    # Weighted region is W -> D (with slope)
    residue[w.:end[1]] <-  residue[w.:end[1]]*3
    if(length(residue) > end[1]){
      tail <- (end[1]+1):length(residue)
      for(j in 1:length(tail)){
        wgt <- 3 - (0.1*j)
        if(wgt < 1){wgt <- 1}
        residue[tail[j]] <- residue[tail[j]]*wgt
      }
    }
    
    # Calculate Reduced Chi-Square for the beat:
    nData <- nrow(dat)    
    nPar <- length(par2) 
    if(par2[1] == par2[2]){    # If baselines are the same, consider them as 1 parameter
      nPar <- nPar - 1
    }
    beat_fit[i] <- (sum(residue*residue) / (nData-nPar)) + as.numeric(penalty)
    
  }
  
  # Summate individual beat ChiSq values:
  temp <- c()
  for(i in 1:length(beat_fit)){
    temp[i] <- beat_fit[i]
  }
  ts_fit <- sum(temp)
  rm(temp)
  
  fit <- list(ts_fit, beat_fit, max_error, NRMSE, aNRMSE)
  return(fit)
}


model2.Rebuild2 <- function(xy,offset,params,invert=TRUE){     
  result <- 1:nrow(xy) * 0.0
  # Creating the excess:
  result <- result + model2.Peak(xy[,1],params[3:5])     # Systolic parameters 
  if (length(params)>=8){
    result <- result + model2.Peak(xy[,1],params[6:8]+c(params[3],0,0))   # Diastolic parameters 
  }
  if (length(params)>=11){
    result <- result + model2.Peak(xy[,1],params[9:11]+c(params[3],0,0))  # Renal parameters
  }
  # Adding decay (config.rate + baseline parameters):
  if (invert){
    result <- model2.Excess.Inv2(xy[,1],result,offset,params[1],params[2],params[3]+1*params[6], config.rate = params[12])   
  }
  return(as.double(result))
}


model2.Excess.Inv2 <- function(time,excess,offset,baselineStart,baselineEnd,timeBase,config.rate){   
  nX <- length(excess)
  if (nX == 0){
    print("Help")
  }
  result <- 1:nX * 0.0
  baseline <- time * 0 + baselineStart
  baseline[which(time > timeBase)] = baselineEnd
  
  # If excess has NAs it will interrupt the reconstruction, remove them:
  temp <- which(is.nan(excess))
  if(length(temp) > 0){
    for(i in 1:length(temp)){
      excess[temp][i] <- 0 
    }
  }
  
  # Adding the decay element to the excess (one value at a time):
  result[1] = excess[1] + (baselineStart + config.rate*(offset-baselineStart)) 
  for (j in 2:nX){  
    result[j] = excess[j] + (baseline[j] + config.rate*(result[j-1]-baseline[j]))  
  }
  return(result)
}


model2.FIX_PAR3 <- function(time, within_beat_params, across_beat_params, debug=FALSE, renal_param, sys_t){
  
  # params: {Baseline, {baseline 2}, t_sys, H_sys, W_sys, {dt_1, H_1, W_1, {dt_2, H_2, W_2}}}
  # across_beat_params: { w[1], t[2], w[2], t[3], w[3] }
  
  nPar <- length(within_beat_params) + length(across_beat_params)
  nData <- length(time)
  
  # Transcribe parameters
  nBase <- 1
  baseline <- c( within_beat_params[1], within_beat_params[1] )
  if (nPar == 6 | nPar == 9 | nPar == 12){     
    baseline[2] = within_beat_params[2]
    nBase <- 2
  }
  
  # par: {base1, {base2}, t[1], h[1], #, #, { h[2], #, #, { h[3], ..., ... }}}
  
  t <- c( within_beat_params[nBase + 1], across_beat_params[2], across_beat_params[4])           # time (systolic = within, diastolic / renal = across)
  h <- c( within_beat_params[nBase + 2], 0, 0 )                                                  # height (systolic = within, diastolic / renal default to 0 unless peaks supplied (see below))
  w <- c( across_beat_params[1], across_beat_params[3], across_beat_params[5])    # width (systolic / diastolic / renal = across)
  hasPeak <- c( TRUE, FALSE, FALSE )
  
  # Calculate penalty only if a peak has been supplied       
  
  if (nPar >= nBase + 7){                                    # Assign diastolic values if peak present
    hasPeak[2] = TRUE
    t[2] <- across_beat_params[2]
    h[2] <- within_beat_params[nBase + 3]
    w[2] <- across_beat_params[3]
  }
  
  if (nPar >= nBase + 10){                                   # Assign renal values if peak present  
    hasPeak[3] = TRUE
    t[3] <- across_beat_params[4]
    h[3] <- within_beat_params[nBase + 4]
    w[3] <- across_beat_params[5]
  }
  
  # Clamp and/or penalize parameters
  penalty <- 0
  
  # 
  tMin <- time[1]   
  tMax <- time[length(time)]       
  
  META_BASELINE_SHIFT <- 1.0    # penalty for how big the gap is between the two baselines
  META_MIN_PEAK_DELAY <- 0.1    # peaks cannot be following one another by less than 0.1ms
  MIN_WIDTH <- c(0.05, 0.05, 0.1) 
  MAX_WIDTH <- c(0.5, 0.45, 0.25)
  
  p <- 1:12*0    # One penalty value for each parameter
  # p: { #, #, t[1], h[1], w[1], t[2], h[2], w[2], t[3], h[3], w[3], across_beat_params[6] }
  
  # Fix height and width for each of the three waves (1:3)      
  for (i in 1:3){
    if (h[i] < 0){                                           # Heights should not be negative
      penalty <- penalty + h[i]*h[i] 
      p[3*i+1] <- h[i]*h[i]                
      h[i] <- 0
    }
    
    if (w[i] < MIN_WIDTH[i] | w[i] > MAX_WIDTH[i]){           # Correct widths as per MIN/MAX_WIDTH
      fixed <- max(MIN_WIDTH[i], min( w[i], MAX_WIDTH[i]))             
      diff <- fixed - w[i]
      penalty <- penalty + diff*diff
      p[3*i+2] <- diff*diff            
      w[i] <- fixed
    }
    
    if(i==3){                                           # Renal peak should be penalized as its amplitude increases
      if( h[3] > (h[1]/50)){                            # As soon as it's amplitude exceeds 2% of the systolic amp. 
        diff <- h[3] - (h[1]/50)
        penalty <- penalty + 2*diff*diff
        p[10] <- p[10] + 2*diff*diff
      }
    }
  }
  
  # Fix time
  
  # Systolic:
  fixed <- max((sys_t - 0.04) , min( t[1], (sys_t + 0.04 )))      # Making sure S peak sits within 40ms of the the peak of the data
  if (debug){
    print(paste("time S: ",tMin," < ",t[1]," < min( ",tMin+1,",",tMax," )"))
  }
  if (t[1] != fixed){
    diff <- fixed - t[1]
    penalty <- penalty + 10^8*diff*diff   
    p[3] <- 10^8*diff*diff
    t[1] <- fixed
  }
  
  # Diastolic:
  fixed <- max( 2 * META_MIN_PEAK_DELAY, min( t[2], tMax - tMin + 0.4 * w[2] ) )   # This stops diastolic time being < 0.2
  if (debug){
    print(paste("time D: ",2 * META_MIN_PEAK_DELAY," < ",t[2]," < ",tMax - tMin + 0.4 * w[2]," )"))
  }
  if (t[2] != fixed){
    # Two peak delays between S and D
    diff <- fixed - t[2]
    if (hasPeak[2]){
      penalty <- penalty + diff*diff  
      p[6] <- diff*diff
    }
    t[2] <- fixed
  }
  
  # Renal:
  fixed <- max( max(META_MIN_PEAK_DELAY, renal_param - 0.02), min( t[3], t[2] - META_MIN_PEAK_DELAY, renal_param + 0.02 ) )   # Stops renal peak being < 0.1 after systolic or < 0.1 before diastolic, and within 20ms of renal param
  if (debug){
    print(paste("time R: ",META_MIN_PEAK_DELAY," < ",t[3]," < ",t[2] - META_MIN_PEAK_DELAY," )"))
  }
  if (t[3] != fixed){
    diff <- fixed - t[3]
    if (hasPeak[3]){
      penalty <- penalty + 5*10^7*diff*diff 
      p[9] <- diff*diff
    }
    t[3] <- renal_param
  }
  
  # Config.rate
  if(across_beat_params[6] > 0.95){
    diff <- across_beat_params[6] - 0.95
    penalty <- penalty + 10^7*diff*diff
    p[12] <- 10^7*diff*diff
    across_beat_params[6] <- 0.95
  }
  
  # Baseline1 shouldn't be above 0:
  if(baseline[1] > 0){
    penalty <- penalty + baseline[1]*baseline[1]
    p[1] <- baseline[1]*baseline[1]
    baseline[1] <- 0
  }
  
  fixedPar <- c( baseline[1:2], t[1], h[1], w[1], t[2], h[2], w[2], t[3], h[3], w[3], across_beat_params[6])
  
  if (debug){
    print(p)
  }
  
  return( c( penalty, fixedPar ) )
}


model2.FixParams3 <- function(data,params, across_beat_params = NULL, debug=FALSE, rp = renal_param, sys_t){
  
  # If across_beat_params have not been provided, extract them from params
  if(is.null(across_beat_params)){    
    across_beat_params <- params[c(5, 6, 8, 9, 11, 12)]
  }
  
  temp <- model2.FIX_PAR3(time = data[, 1], within_beat_params = params[c(1:4, 7, 10)], across_beat_params, debug = F, renal_param = rp, sys_t)  
  return( temp[2:length(temp)] )     # first value of temp is penalty
} 



# Make Simplex 2 (for across beat parameters only):
simplex.MakeSimplex2 <- function(data,param,f,inScale,directions=NULL,inTol=-1, optional=NULL,debug=FALSE, beat_vector = beat_vector,
                                 beat = beat, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w){
  ########################################################################################################################################
  # MakeSimplex2 iterates on the across-beat starting parameters, refining them to give to give the simplex a good starting position. 
  
  # Inputs: 
  # data ()
  # param
  # f
  # inScale
  # directions
  # inTol
  # optional
  # debug
  # beat_vector
  # beat
  # renal_param
  # dias_param
  # sys_time
  # w
  
  
  # Outputs:
  # ?? ()
  ########################################################################################################################################
  
  if(debug){print("MakeSimplex -- debug")}                # this seems unecessary now..         
  
  nPar <- length(param)                  # define number of parameters
  nScale <- length(inScale)               # what is inscale?
  
  
  if (nScale == 0)       # Not sure what this is doing.. 
  {
    scale <- 1:nPar * 0 + 1
  } else if (nScale == 1){
    scale <- 1:nPar * 0 + inScale
  } else if (length(inScale) == nPar){
    scale <- inScale
  } else {
    #print("Invalid scale vector length")
    return("Error: Invalid scale vector length")
  }
  
  
  if (length(inTol) == 1 & inTol > 0){   # If no tolerance is provided, the tol is defined as the value of the point / vertex provided
    tol <- inTol[1]
  } else {
    tol <- min(1,f(data,param, beats = beat_vector, beat = beat, renal_param = renal_param,
                   dias_param = dias_param, sys_time = sys_time, w = w))  
  }
  
  
  chiSq <- 1:(nPar+1) * 0.0    # Creating a vector of ChiSq values... 
  
  chiSq[1] <- f(data,param,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param,    # ChiSq[1 is the fit when no parameters are changed... 
                dias_param = dias_param, sys_time = sys_time, w = w)                                            
  
  if (debug){ print(paste("Root chi-squared:",chiSq[1]))}   # is this line still necessary?
  
  result <- matrix(nrow=nPar+1,ncol=nPar)      # Create a matrix of parameters, with as many columns as parameters, and nPar + 1 rows (each will be one point of the simplex)
  result[1,] <- as.double(param)    # Fill the first row with the inputted parameters (our best guess so far)... 
  
  useDirections = !is.null(directions)         # I don't think we ever do use directions...?
  
  if (useDirections){ useDirections <- nrow(directions) == nPar & ncol(directions) == nPar}     # ditto, we are not using directions.. but is it necessary?
  
  
  for (i in c(5, 6, 8, 9, 11, 12)){                  # Create a for loop for each of the across beat parameters (the ones we are refining)
    
    if (debug){ print(paste("Parameter",i)) }  # We are not debugging now so ?necessary
    
    tParam <- param    # create a new vector of parameters to be tweaked / tested, tParam
    
    # Pick a direction
    delta <- 1:nPar * 0              # delta is a finite increment by which to increase or decrease a given parameter value. We create a vector for it here. 
    
    if (useDirections){
      delta <- scale[i] * directions[i,]      # this if is redunant since we are not using directions...
    } else {
      delta[i] <- scale[i]       # scale relates to inScale, somehow.. and is used to set the delta level
    }
    
    tParam <- param - delta                                   # This part tries tweaking each parameter up or down
    
    chiSqMinus <- f(data,tParam,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param,  # we run the test parameters through chisq3 to get a chisq value for them
                    dias_param = dias_param, sys_time = sys_time, w = w)         
    
    tParam <- param + delta                                   # This time delta is added rather than subtracted, and again the resulting Chisq change is found
    
    chiSq[i+1] <- f(data,tParam,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param,
                    dias_param = dias_param, sys_time = sys_time, w = w)            
    
    
    if (debug){                       # redundant now? need to go through specifically looking for debug functions and get a sense of what it is doing
      print("Select direction:")
      print(paste("chi^2(",param[i] - delta[i],") =",chiSqMinus))
      print(paste("chi^2(",param[i],") =",chiSq[1]))
      print(paste("chi^2(",param[i] + delta[i],") =",chiSq[i+1]))
      print("---")
    }
    
    if (chiSqMinus < chiSq[i+1]){      # If going down by delta is better than going up by delta, 
      delta <- -delta                  # then replace chiSq[i+1] with the lower score (ChiSqMinus)
      tParam <- param + delta
      chiSq[i+1] <- chiSqMinus
    }
    
    iKill <- 10       # defining the number of iterations that a given parameter will be refined over
    
    if (chiSq[i+1] < chiSq[1]){         # If the new fit is better than the old fit (with no parameters changed), continue to go in the direction that improved the fit
      if (debug){ print("Extending as best point") }
      while (chiSq[i+1] < chiSq[1] + tol){                 # Chisquare keeps getting iterated here (for 10 iterations)
        delta <- 2*delta   
        tParam <- param + delta
        oldScore <- chiSq[i+1]          # The current best fit is called 'old score'
        chiSq[i+1] <- f(data,tParam,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w)   # The new fit is now designated ChiSq[i+1]
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }
        if (chiSq[i+1] > oldScore){   # Check if the new fit is worse than current fit
          tParam <- param + 0.5*delta   # If so, make delta what it was one iteration previous (undoing the *2)
          chiSq[i+1] <- oldScore        # and redesignate old score to chiSq[i+1]
          break
        }
        #print(paste(i,"-",delta,":",chiSq[i+1]))
        iKill <- iKill - 1             # If the new fit is not worse than the current fit, knock of one on the interation count,
        if (iKill < 0){                # and keep iterating until either the next fit is worse (and break is called), or iKill ends (and break is also called)
          break
        }
      }
    } else if (chiSq[i+1] < chiSq[1] + tol){    # If the new fit is not better than the old fit, is it at least better than the old fit + tol?
      if (debug){ print("Extending below tolerance") }
      while (chiSq[i+1] < chiSq[1] + tol){
        delta <- 2*delta
        tParam <- param + delta
        oldScore <- chiSq[i+1]
        chiSq[i+1] <- f(data,tParam,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w)
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }        
        if (chiSq[i+1] - oldScore < oldScore - chiSq[1]){  
          tParam <- param + 0.5*delta
          chiSq[i+1] <- oldScore
          break
        }
        iKill <- iKill - 1
        if (iKill < 0){
          if(i == 9){
            tParam[9] <- renal_param  # Ignore renal times that can't optomize
            break
          } 
          print(c("simplex constructed as per original parameter"))
          break
          #print("Failed to construct simplex")
          #return(paste("Error: param[",i,"]",sep=""))
        }
      }
    } else {
      if (debug){ print("Shrinking above tolerance") }
      while (chiSq[i+1] > chiSq[1] + tol){              # If the new fit is much worse than the original, reduce the size of delta
        delta <- 0.5*delta
        tParam <- param + delta
        lastChiSq <- chiSq[i+1]
        chiSq[i+1] <- f(data,tParam,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w)
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }
        #print(paste(i,"-",delta,":",chiSq[i+1]))
        if (iKill < 0 & (chiSq[i+1]-chiSq[1]) > 0.75 * (lastChiSq-chiSq[1])){
          if(i == 9){
            tParam[9] <- renal_param  # Ignore renal times that can't optomize
            next
          } 
          print(c("simplex constructed as per original parameter"))
          next
          #print("Failed to construct simplex")
          #return(paste("Error: param[",i,"]",sep=""))   
        }
        iKill <- iKill - 1
      }
      tParam <- param + 0.5 * delta
    }
    
    if(debug){ print(paste("Param[",i,"] =",tParam[i]))}
    result[i+1,] = as.double(tParam) 
    
  }
  
  if (debug){ print("/MakeSimplex") }
  return(result)
}

# [20:15, 19/06/21] Simon Williamson
# Thanks Craig, sorry it's taken me a while to get back round to this. 
# I think I understand how things are working from what you've explained. So I suppose the question is why have tolerance at all? Why not have each vertex extend or shrink until it is below the original vertex, rather than below the original vertex + tol?
# Trying to answer this myself, I suppose there's the possibility that some parameters cannot be improved upon, in which case the test parameter can at least shrink until it is minimally worse than the original? And that minimum is defined by tol?



# Make simplex 3 (for within-beat parameters only)
simplex.MakeSimplex3 <- function(ppg, param,f,inScale, directions=NULL, inTol=-1, optional=NULL, debug=FALSE, beat_vector = beat_vector, renal_param, dias_param = dias_param, sys_time, w){
  
  if(debug){print("MakeSimplex -- debug")}
  nPar <- length(param)
  nScale <- length(inScale)
  if (nScale == 0)
  {
    scale <- 1:nPar * 0 + 1
  } else if (nScale == 1){
    scale <- 1:nPar * 0 + inScale
  } else if (length(inScale) == nPar){
    scale <- inScale
  } else {
    #print("Invalid scale vector length")
    return("Error: Invalid scale vector length")
  }
  if (length(inTol) == 1 & inTol > 0){
    tol <- inTol[1]
  } else {
    tol <- min(1,f(data = ppg, params = param, beats = beat_vector, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w))    
  }    
  
  
  chiSq <- 1:(nPar+1) * 0.0
  chiSq[1] <- f(data = ppg, param, beats = beat_vector, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w)
  if (debug){ print(paste("Root chi-squared:",chiSq[1]))}
  
  result <- matrix(nrow=nPar+1,ncol=nPar)
  result[1,] <- as.double(param)
  
  useDirections = !is.null(directions)
  if (useDirections){ useDirections <- nrow(directions) == nPar & ncol(directions) == nPar }
  
  for(i in c(1:4, 7, 10)){    # within-beat parameters only  
    if (debug){ print(paste("Parameter",i)) }
    tParam <- param
    
    # Pick a direction
    delta <- 1:nPar * 0
    if (useDirections){
      delta <- scale[i] * directions[i,]
    } else {
      delta[i] <- scale[i]
    }
    
    if(i == 3){
      delta <- delta/4
    }
    
    tParam <- param - delta                                   # This part tries tweaking each parameter up or down, 
    chiSqMinus <- f(data = ppg, params = tParam, beats = beat_vector, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w)            # tParam = test parameter. 
    tParam <- param + delta                                   # The chisquare (goodness of fit) is calculated for each direction, 
    chiSq[i+1] <- f(data = ppg, params = tParam, beats = beat_vector, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w)            # the direction with the smaller value is chosen. 
    
    if (debug){
      print("Select direction:")
      print(paste("chi^2(",param[i] - delta[i],") =",chiSqMinus))
      print(paste("chi^2(",param[i],") =",chiSq[1]))
      print(paste("chi^2(",param[i] + delta[i],") =",chiSq[i+1]))
      print("---")
    }
    
    if (chiSqMinus < chiSq[i+1]){
      delta <- -delta
      tParam <- param + delta
      chiSq[i+1] <- chiSqMinus
    }
    
    iKill <- 10    
    
    if (chiSq[i+1] < chiSq[1]){
      if (debug){ print("Extending as best point") }
      while (chiSq[i+1] < chiSq[1] + tol){                 # Chisquare keeps getting iterated here (for 10 iterations)
        delta <- 2*delta
        tParam <- param + delta
        oldScore <- chiSq[i+1]
        chiSq[i+1] <- f(data = ppg,tParam, beats = beat_vector, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w)
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }
        if (chiSq[i+1] > oldScore){
          tParam <- param + 0.5*delta
          chiSq[i+1] <- oldScore
          break
        }
        #print(paste(i,"-",delta,":",chiSq[i+1]))
        iKill <- iKill - 1
        if (iKill < 0){
          break
        }
      }
    } else if (chiSq[i+1] < chiSq[1] + tol){
      if (debug){ print("Extending below tolerance") }
      while (chiSq[i+1] < chiSq[1] + tol){
        delta <- 2*delta
        tParam <- param + delta
        oldScore <- chiSq[i+1]
        chiSq[i+1] <- f(data = ppg, tParam, beats = beat_vector, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w)
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }
        if (chiSq[i+1] - oldScore < oldScore - chiSq[1]){
          tParam <- param + 0.5*delta
          chiSq[i+1] <- oldScore
          break
        }
        iKill <- iKill - 1
        if (iKill < 0){
          #print("Failed to construct simplex")
          #return(paste("Error: param[",i,"]",sep=""))
          print(c("Failed to construct simplex within 10 iterations for parameter", i, "defaulting to inputted value"))
          tParam[i] <- param[i]
          break  # this was next 
        }
      }
    } else {
      if (debug){ print("Shrinking above tolerance") }
      while (chiSq[i+1] > chiSq[1] + tol){
        delta <- 0.5*delta
        tParam <- param + delta
        lastChiSq <- chiSq[i+1]
        chiSq[i+1] <- f(data = ppg,tParam, beats = beat_vector, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w)
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }
        #print(paste(i,"-",delta,":",chiSq[i+1]))
        if (iKill < 0 & (chiSq[i+1]-chiSq[1]) > 0.75 * (lastChiSq-chiSq[1])){
          print(c("Failed to construct simplex within 10 iterations for parameter", i, "defaulting to inputted value"))
          #return(paste("Error: param[",i,"]",sep=""))
          tParam[i] <- param[i]
          break # this was next
        }
        iKill <- iKill - 1
      }
      tParam <- param + 0.5 * delta
    }
    
    if(debug){ print(paste("Param[",i,"] =",tParam[i]))}
    result[i+1,] = as.double(tParam)
  }
  
  if (debug){ print("/MakeSimplex") }
  return(result)
}



simplex.Run2 <- function(data = ppg,simplexParam = mat, f = model2.ChiSq3, optional=NULL, beat_vector = beat_vector, ms = simplex_iterations, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w, run = NULL){
  ########################################################################################################################################
  # Simplex.Run2 is the function that intiates the process of simplex relfection through the multi-dimensional landscape. 
  
  # Inputs: 
  # data ()
  # simplexParam (a matrix... )
  # f (the function used to assess goodness of fit by comparing the fitted wave to the original data wave)
  # optional 
  # beat_vector (an index of beats to be modeled and their x-coordinates to extract from the PPG time series)
  # renal_param (the intially identified renal timing, used as an anchor point for the model)
  # dias_param (the initially identified diastolic timing, used as an achor point for the model / about which the model is constrained)
  # sys_time (the timing of systole (about which the model is constrained))
  # w ()
  # run (indicates which simplex run is being conducted when printed)
  
  
  # Outputs:
  # Basically the set of parameters which best optimizes the function model2.ChiSq
  ########################################################################################################################################
  MAX_STEP <- ms                                               # Maximum number of allowed function evaluations (numerical recipes)
  FTOL <- 1e-5                                  # Tolerance, or fraction of tolerance possibly related to machine precision - need to clarify
  
  debugRtol <- 1:(MAX_STEP+1) * 0.0             # Not sure if we use any of these at any point... 
  debugMin <- 1:(MAX_STEP+1) * 0.0
  debugMax <- 1:(MAX_STEP+1) * 0.0
  
  result <- simplexParam                         # Passing in the 66*66 matrix, which will be the output once iterated on
  nPar <- ncol(result)                      
  chiSq <- 0:nPar * 0.0             # A separate chisq value exists for each parameter, so 66... 
  
  for (i in 1:(nPar+1)){                                    # Find out the ChiSq value for each row from result, there are 66 + 1 rows
    chiSq[i] <- f(data, params = NULL, optional=NULL, a = result[i, ], beats = beat_vector,
                  renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w)
  }                         # Having this vector gives us the values of each point in the simplex
  
  for (iStep in 1:MAX_STEP){                             # beginning of downhill simplex
    extrema <- simplex.SortHighLow(chiSq)                # Finds the results which give the highest, 2nd highest and lowest ChiSq
    low <- extrema[1]
    nHigh <- extrema[2] 
    high <- extrema[3]
    
    if(!is.null(run)){        # just prints run number and iteration number, so you have an idea of it's progress when running it
      print(run)
    }
    print(iStep)
    
    chiSqMax <- chiSq[high]                     # redefine highest and lowest points...
    chiSqMin <- chiSq[low]                      # These should be the ChiSq values for individual rows (the best and worst points), so how can it be compared to the 'score' that resresents the entire new simplex?
    
    print(chiSqMax)
    
    #print(paste("chi^2_min =",chiSqMin))
    #print(paste("argMax = ",high,"[",chiSqMax,"]",sep=""))
    
    rtol <- 2 * (chiSqMax - chiSqMin)/(chiSqMax + chiSqMin + 1e-10)   # measure of how much better high is from low...
    if (rtol < FTOL){
      bestParam <- result[low,]                     # Presumably if the difference in ChiSq (max vs min) is significant, 
      result[low,] <- result[1,]                    # the result that was changed to give the lowest ChiSq gets designated 'best Param'
      result[1,] <- bestParam                       # The best performing row gets upgraded to first row (swapped with what is there currently)
      return(result) 
    }
    
    # So at this point you have input the matrix, identified chiSq values for each row, and moved the row 
    # that generates the best Chisq to the top (at least for iteration 1). 
    
    debugRtol[iStep] <- rtol                         # ? necessary lines?
    debugMin[iStep] <- chiSqMin
    debugMax[iStep] <- chiSqMax
    
    factor <- -1       # not sure what this does
    node <- simplex.HypoCentre(result,high)        # Hypocentre outputs all the parameters that are not the worst  - IS THIS A SINGLE ROW OR A MATRIX?
    apex <- result[high,]                          # Apex is the worst row / point on the simplex
    test <- node - (apex - node)                   # This represents the flipping of the triangle; the whole parameter set is reversed in the direction away from the worst ChiSq point (literally, we are subtracting the difference between each row and the worst row from each row, causing all rows to be less similar to the worst ro)
    score <- f(data, params = rep(0, 12),optional=optional, a = test, beats = beat_vector,
               renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w)   # QUESTION: is this a single value? Does it evaluate every row of the simplex?
    
    if (score < chiSqMin){                          # If flipping improves the ChiSq, try extending further in the same direction (increase the scale by 2) i.e reflection and expansion
      test2 <- node - 2 * (apex - node)
      score2 <- f(data, params = rep(0, 12),optional=optional, a = test2, beats = beat_vector, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w)
      if (score2 >= score){                       # If reflecting a further distance is better than reflecting alone, do that
        # Reflect
        #print(paste("Reflecting",high,": chi^2 ",chiSqMax,"->",score,sep=""))
        result[high,] <- test       # test is replacing one row, so is test just one row?
        chiSq[high] <- score
      } else {
        # Reflect and grow
        #print(paste("Reflect-stretching",high,": chi^2 ",chiSqMax,"->",score2,sep=""))
        result[high,] <- test2
        chiSq[high] <- score2
      }
    } else if (score >= chiSq[nHigh]) {              # If reflecting is not beneficial, try shrinking instead of reflecting
      # Test for shrink with optional reflection
      factor <- 0.5
      if (score < chiSqMax)
      {
        factor <- -0.5
      }
      test2 <- node + factor * (apex - node)
      score2 <- f(data, params = rep(0, 12),optional=optional, a = test2, beats = beat_vector, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w)
      if (score2 < chiSq[nHigh]){
        # Shrink (possibly reflecting)
        #print(paste("Shrinking",high,": chi^2 ",chiSqMax,"->",score2,sep=""))
        result[high,] <- test2
        chiSq[high] <- score2
      } else {
        # Shrink all
        for (i in 1:(nPar+1)){
          if (i != low){
            result[i,] <- 0.5 * (result[i,] + result[low,])
            chiSq[i] <- f(data, params = rep(0, 12),optional=optional, a = result[i, ], beats = beat_vector, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w)
          }
        }
        #print(paste("General contraction: chi^2 ",chiSqMax,"->",max(chiSq),sep=""))
      }
    } else {
      # Reflect
      #print(paste("Reflecting*",high,": chi^2 ",chiSqMax,"->",score,sep=""))
      result[high,] <- test
      chiSq[high] <- score
    }
  }
  
  extrema <- simplex.SortHighLow(chiSq)
  low <- extrema[1]
  bestParam <- result[low,]
  result[low,] <- result[1,]
  result[1,] <- bestParam
  
  chiSqMax <- chiSq[extrema[3]]
  chiSqMin <- chiSq[low]
  rtol <- 2 * (chiSqMax - chiSqMin)/(chiSqMax + chiSqMin + 1e-10)   # Why does rtol need to be redefined here?
  debugRtol[MAX_STEP+1] <- rtol
  debugMin[MAX_STEP+1] <- chiSqMin      # Are these debug lines needed?
  debugMax[MAX_STEP+1] <- chiSqMax
  # plot(debugMax,type='l')
  # lines(debugMin)
  
  
  print(paste("Terminated downhill simplex after",MAX_STEP,"iterations."))
  print(paste("rtol =",rtol))
  return(result)
}


simplex.HypoCentre <- function(mat_Param,index){
  nPar <- ncol(mat_Param)
  
  result <- 1:nPar * 0.0
  for (i in 1:(nPar+1)){
    if (i != index){
      result <- result + mat_Param[i,]
    }
  }
  return( result / nPar )
}


simplex.SortHighLow <- function(vec_ChiSq){
  nPar <- length(vec_ChiSq)
  
  low <- 1
  high <- 1
  nHigh <- 2
  if (vec_ChiSq[2] > vec_ChiSq[1]){
    high <- 2
    nHigh <- 1
  }
  
  for (i in 2:nPar){
    if (vec_ChiSq[i] < vec_ChiSq[low]){
      low <- i
    }
    if (vec_ChiSq[i] > vec_ChiSq[high]){
      nHigh <- high
      high <- i
    }
    if (i != high & vec_ChiSq[i] > vec_ChiSq[nHigh]){
      nHigh <- i
    }
  }
  
  return(c(low,nHigh,high))  
}


# This is an exceptionally inefficient use of space, and needs to be revisited:

PlotRejects <- function(rejected_waves_list1, rejected_waves_list3){
  
  if(length(rejected_waves_list3) > 0){
    len <- 2
  }else{
    len <- 1
  }
  
  participant_extra_long_waves <- c()
  for(k in 1:len){
    
    if(k == 1){
      new_vec <- rejected_waves_list1
    }else{
      new_vec <- rejected_waves_list3
    }
    
    extra_long_waves <- c()
    for(i in 1:length(new_vec)){
      if(!is.null(new_vec[[i]][[1]])){
        extra_long_waves[[i]] <- new_vec[[i]][[1]]
      }
    }
    non_null_names <- which(!sapply(extra_long_waves, is.null))
    extra_long_waves <- extra_long_waves[non_null_names]
    if(length(non_null_names) > 0){names(extra_long_waves) <- non_null_names}
    
    if(length(extra_long_waves) > 0){
      for(j in 1:length(extra_long_waves)){
        
        no_of_rejected_waves <- length(extra_long_waves[[j]])
        
        paticip <- as.numeric(rownames(summary(extra_long_waves[j])))
        
        participant_extra_long_waves[paticip] <- length(participant_extra_long_waves[paticip]) + no_of_rejected_waves 
      }
    }
  }
  if(length(participant_extra_long_waves) > 0){
    for(i in 1:length(participant_extra_long_waves)){
      if(is.na(participant_extra_long_waves[i])){
        participant_extra_long_waves[i] <- 0
      }
    }
  }
  
  
  participant_extra_short_waves <- c()
  for(k in 1:len){
    
    if(k == 1){
      new_vec <- rejected_waves_list1
    }else{
      new_vec <- rejected_waves_list3
    }
    
    extra_short_waves <- list()
    for(i in 1:length(new_vec)){
      if(!is.null(new_vec[[i]][[2]])){
        extra_short_waves[[i]] <- new_vec[[i]][[2]]
      }
    }
    non_null_names <- which(!sapply(extra_short_waves, is.null))
    extra_short_waves <- extra_short_waves[non_null_names]
    if(length(non_null_names) > 0){names(extra_short_waves) <- non_null_names}
    
    if(length(extra_short_waves) > 0){
      for(j in 1:length(extra_short_waves)){
        
        no_of_rejected_waves <- length(extra_short_waves[[j]])
        
        paticip <- as.numeric(rownames(summary(extra_short_waves[j])))
        
        participant_extra_short_waves[paticip] <- length(participant_extra_short_waves[paticip]) + no_of_rejected_waves 
      }
    }
  }
  if(length(participant_extra_short_waves) > 0){
    for(i in 1:length(participant_extra_short_waves)){
      if(is.na(participant_extra_short_waves[i])){
        participant_extra_short_waves[i] <- 0
      }
    }
  }
  
  participant_double_segments <- c()
  for(k in 1:len){
    
    if(k == 1){
      new_vec <- rejected_waves_list1
    }else{
      new_vec <- rejected_waves_list3
    }
    
    double_segments <- list()
    for(i in 1:length(new_vec)){
      if(!is.null(new_vec[[i]][[3]])){
        double_segments[[i]] <- new_vec[[i]][[3]]
      }
    }
    non_null_names <- which(!sapply(double_segments, is.null))
    double_segments <- double_segments[non_null_names]
    if(length(non_null_names) > 0){names(double_segments) <- non_null_names}
    
    if(length(double_segments) > 0){
      for(j in 1:length(double_segments)){
        
        no_of_rejected_waves <- length(double_segments[[j]])
        
        paticip <- as.numeric(rownames(summary(double_segments[j])))
        
        participant_double_segments[paticip] <- length(participant_double_segments[paticip]) + no_of_rejected_waves 
      }
    }
  }
  if(length(participant_double_segments) > 0){
    for(i in 1:length(participant_double_segments)){
      if(is.na(participant_double_segments[i])){
        participant_double_segments[i] <- 0
      }
    }
  }
  
  
  
  participant_systolic_endings <- c()
  for(k in 1:len){
    
    if(k == 1){
      new_vec <- rejected_waves_list1
    }else{
      new_vec <- rejected_waves_list3
    }
    
    systolic_endings <- list()
    for(i in 1:length(new_vec)){
      if(!is.null(new_vec[[i]][[4]])){
        systolic_endings[[i]] <- new_vec[[i]][[4]]
      }
    }
    non_null_names <- which(!sapply(systolic_endings, is.null))
    systolic_endings <- systolic_endings[non_null_names]
    if(length(non_null_names) > 0){names(systolic_endings) <- non_null_names}
    
    if(length(systolic_endings) > 0){
      for(j in 1:length(systolic_endings)){
        
        no_of_rejected_waves <- length(systolic_endings[[j]])
        
        paticip <- as.numeric(rownames(summary(systolic_endings[j])))
        
        participant_systolic_endings[paticip] <- length(participant_systolic_endings[paticip]) + no_of_rejected_waves 
      }
    }
    
  }
  if(length(participant_systolic_endings) > 0){
    for(i in 1:length(participant_systolic_endings)){
      if(is.na(participant_systolic_endings[i])){
        participant_systolic_endings[i] <- 0
      }
    }
  }
  
  
  
  participant_drops_below_o <- c()
  for(k in 1:len){
    
    if(k == 1){
      new_vec <- rejected_waves_list1
    }else{
      new_vec <- rejected_waves_list3
    }
    
    drops_below_o <- list()
    for(i in 1:length(new_vec)){
      if(!is.null(new_vec[[i]][[5]])){
        drops_below_o[[i]] <- new_vec[[i]][[5]]
      }
    }
    non_null_names <- which(!sapply(drops_below_o, is.null))
    drops_below_o <- drops_below_o[non_null_names]
    if(length(non_null_names) > 0){names(drops_below_o) <- non_null_names}
    
    if(length(drops_below_o) > 0){
      for(j in 1:length(drops_below_o)){
        
        no_of_rejected_waves <- length(drops_below_o[[j]])
        
        paticip <- as.numeric(rownames(summary(drops_below_o[j])))
        
        participant_drops_below_o[paticip] <- length(participant_drops_below_o[paticip]) + no_of_rejected_waves 
      }
    }
    
  }
  if(length(participant_drops_below_o) > 0){
    for(i in 1:length(participant_drops_below_o)){
      if(is.na(participant_drops_below_o[i])){
        participant_drops_below_o[i] <- 0
      }
    }
  }
  
  
  participant_hrsd_waves <- c()
  for(k in 1:len){
    
    if(k == 1){
      new_vec <- rejected_waves_list1
    }else{
      new_vec <- rejected_waves_list3
    }
    
    hrsd_waves <- list()
    for(i in 1:length(new_vec)){
      if(!is.null(new_vec[[i]][[6]])){
        hrsd_waves[[i]] <- new_vec[[i]][[6]]
      }
    }
    non_null_names <- which(!sapply(hrsd_waves, is.null))
    hrsd_waves <- hrsd_waves[non_null_names]
    if(length(non_null_names) > 0){names(hrsd_waves) <- non_null_names}
    
    if(length(hrsd_waves) > 0){
      for(j in 1:length(hrsd_waves)){
        
        no_of_rejected_waves <- length(hrsd_waves[[j]])
        
        paticip <- as.numeric(rownames(summary(hrsd_waves[j])))
        
        participant_hrsd_waves[paticip] <- length(participant_hrsd_waves[paticip]) + no_of_rejected_waves 
      }
    }
    
  }
  if(length(participant_hrsd_waves) > 0){
    for(i in 1:length(participant_hrsd_waves)){
      if(is.na(participant_hrsd_waves[i])){
        participant_hrsd_waves[i] <- 0
      }
    }
  }
  
  
  participant_outlier_waves <- c()
  for(k in 1:len){
    
    if(k == 1){
      new_vec <- rejected_waves_list1
    }else{
      new_vec <- rejected_waves_list3
    }
    
    outlier_waves <- list()
    for(i in 1:length(new_vec)){
      if(!is.null(new_vec[[i]][[7]])){
        outlier_waves[[i]] <- new_vec[[i]][[7]]
      }
    }
    non_null_names <- which(!sapply(outlier_waves, is.null))
    outlier_waves <- outlier_waves[non_null_names]
    if(length(non_null_names) > 0){names(outlier_waves) <- non_null_names}
    
    if(length(outlier_waves) > 0){
      for(j in 1:length(outlier_waves)){
        
        no_of_rejected_waves <- length(outlier_waves[[j]])
        
        paticip <- as.numeric(rownames(summary(outlier_waves[j])))
        
        participant_outlier_waves[paticip] <- length(participant_outlier_waves[paticip]) + no_of_rejected_waves 
      }
    }
    
  }
  if(length(participant_outlier_waves) > 0){
    for(i in 1:length(participant_outlier_waves)){
      if(is.na(participant_outlier_waves[i])){
        participant_outlier_waves[i] <- 0
      }
    }
  }
  
  # Make them all the same length:
  if(length(participant_extra_long_waves) < length(Participants)){
    diff <- length(Participants) - length(participant_extra_long_waves)
    participant_extra_long_waves <- c(participant_extra_long_waves, rep(0, diff))
  }
  if(length(participant_extra_short_waves) < length(Participants)){
    diff <- length(Participants) - length(participant_extra_short_waves)
    participant_extra_short_waves <- c(participant_extra_short_waves, rep(0, diff))
  }
  if(length(participant_double_segments) < length(Participants)){
    diff <- length(Participants) - length(participant_double_segments)
    participant_double_segments <- c(participant_double_segments, rep(0, diff))
  }
  if(length(participant_systolic_endings) < length(Participants)){
    diff <- length(Participants) - length(participant_systolic_endings)
    participant_systolic_endings <- c(participant_systolic_endings, rep(0, diff))
  }
  if(length(participant_drops_below_o) < length(Participants)){
    diff <- length(Participants) - length(participant_drops_below_o)
    participant_drops_below_o <- c(participant_drops_below_o, rep(0, diff))
  }
  if(length(participant_hrsd_waves) < length(Participants)){
    diff <- length(Participants) - length(participant_hrsd_waves)
    participant_hrsd_waves <- c(participant_hrsd_waves, rep(0, diff))
  }
  if(length(participant_outlier_waves) < length(Participants)){
    diff <- length(Participants) - length(participant_outlier_waves)
    participant_outlier_waves <- c(participant_outlier_waves, rep(0, diff))
  }
  # Find total of all rejected beats:
  total_rejected_beats <- participant_extra_long_waves + participant_extra_short_waves + participant_double_segments + 
    participant_systolic_endings + participant_drops_below_o + participant_hrsd_waves + participant_outlier_waves
  
  
  # Plot:
  plot(participant_extra_long_waves, t = "l", ylim = c(0, 30), xlim = c(1, 112), ylab = "rejected beats (absolute)", xlab = "participants", lty = "dotted", lwd = 1.5)
  lines(participant_extra_short_waves, col = "red", lty = "dotted", lwd = 1.5)
  lines(participant_double_segments, col = "blue", lty = "dotted", lwd = 1.5)
  lines(participant_systolic_endings, col = "green", lty = "dotted", lwd = 1.5)
  lines(participant_drops_below_o, col = "orange", lty = "dotted", lwd = 1.5)
  lines(participant_hrsd_waves, col = "brown", lty = "dotted", lwd = 1.5)
  lines(participant_outlier_waves, col = "purple", lty = "dotted", lwd = 1.5)
  lines(total_rejected_beats, lwd = 1)
  
}


PlotWavesCarriedForward <- function(waves_carried_forward1, waves_carried_forward3){
  
  test_vec1 <- waves_carried_forward1
  waves_carried_over <- c()
  for(i in 1:length(test_vec1)){
    if(!is.null(test_vec1[[i]])){
      waves_carried_over[i] <- test_vec1[[i]]
    }
  }
  
  if(length(waves_carried_forward3) > 0){
    test_vec2 <- waves_carried_forward3
    waves_carried_over2 <- c()
    for(i in 1:length(test_vec2)){
      if(!is.null(test_vec2[[i]])){
        waves_carried_over2[i] <- test_vec2[[i]]
      }
    }
    waves_carried_over <- c(waves_carried_over, waves_carried_over2)
  }
  
  
  # Currently representing all 2mg times series:
  hist(waves_carried_over, breaks = 30, xlim = c(0, 200))
  
}