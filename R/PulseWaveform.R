lab.time = "time (s)"
lab.ppg = "Detrended"
config.rate = 0.89  # 0.99 for archetypal waves
const.pi = 3.1415926535897932384626433

# PPG function list:

# 1. Preproc  (note for bioradio data only)
# 2. find_w
# 3. find_u_v
# 4. find_o
# 5. preclean_wuv
# 6. Baseline
# 7. clean_wuv
# 8. sep_beats
# 9. find_average
# 10. find_sd
# 11. diast_pk
# 12. osnd_of_average
# 13. feature_extract


preproc <- function(data){
  ########################################################################################################################################
  # preproc undertakes a pre-processing routine specific to Bioradio raw data. This entails downsampling the data and the 'undetrending'
  # it. Downsampling removes repeated values (the BioRadio device provides 250 samples per second, but the PPG is only sampled 75 times
  # per second), whilst undetrending reverses a detrending function inherent to Bioradio hardware. To expand on the latter, analysis of
  # device output indicates that the PPG signal is detrended by application of the following formula:
  #
  #               OUT[i] = 80 + (OUT[i-1]-80) * 0.96875 + (IN[i] - [IN[i-1])
  #
  # where the constant 0.96875 is an approximation fitted to the data. Individual pulse events are more comprehensible if the detrending
  # is not used.
  #
  # Inputs:
  # data (the raw Bioradio data)
  #
  # Outputs:
  # undetrended (downsampled and undetrended data)
  ########################################################################################################################################

  dat <- data[!(data$PPG.PulseOx1=='NaN'),]                                                   # removing all rows with NaN values
  list <- rle(dat$PPG.PulseOx1)
  ID <- rep(1:length(list$values), times = list$lengths)                                      # Assigning an ID to each new unique value
  data_downsampled <- c()

  nSrc <- nrow(dat)
  iDst <- 1
  iSrc <- 1
  iVal <- 1
  print("Deduplicating data...")
  if(TRUE){                                                                                   # Latest method
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
      if (list$lengths[iVal] <= 4){                                                           # If up to 4 repeated values, save the first only
        iSrc <- iSrc + list$lengths[iVal]
        iVal <- iVal + 1
      } else {                                                                                # If more than 4 repeated values, save two (assume genuine repeat)
        iSrc <- iSrc + 4
        list$lengths[iVal] <- list$lengths[iVal] - 4
      }
    }
    data_downsampled <- data_downsampled[1:(iDst-1),]                                         # Trim downsampled data to size
  }else{                                                                                      # Alternative (older) method may be of use (though less computationally efficient)
    data2 <- cbind(dat, ID)
    data_downsampled <-c()
    for (i in 1:max(ID)){
      sub.data <- dplyr::filter(data2, ID == i)
      if(nrow(sub.data) <= 4){
        data_downsampled <- rbind(data_downsampled, sub.data[1,])
      }else if(nrow(sub.data) > 4 ){data_downsampled <- rbind(data_downsampled,
                                                              sub.data[1,], sub.data[5,])}
    }
  }

  print("Removing DC blocker...")

  undetrended <- replicate(length(data_downsampled$PPG.PulseOx1)-1, 0)                        # Undetrending
  undetrended <- c(data_downsampled$PPG.PulseOx1[1],undetrended)
  for(i in 2:length(data_downsampled$PPG.PulseOx1)){
    undetrended[i]<-((data_downsampled$PPG.PulseOx1[i]-80) -
                       ((data_downsampled$PPG.PulseOx1[i-1]-80) * 0.96875) +
                       (undetrended[i-1]))
  }
  print("Done")
  return(undetrended)
}



find_w <- function(d1p, deriv1, sp, sr, pk_thrshd){
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

    if(length(windowPks) == 2 & windowInflxX[windowPks[2]] - windowInflxX[windowPks[1]] > sr){                                                  # If the first two peaks indicate a heart rate < 60 bpm,
      confirmWindow <- data.frame(d1InflxX[windowPks[1]+2]:d1InflxX[windowPks[2]-2],                                                            # ensure a peak between them has not been missed with the current threshold
                                  deriv1[d1InflxX[windowPks[1]+2]:d1InflxX[windowPks[2]-2]])
      inflxConformY <- d1InflxY[which(d1InflxX > confirmWindow[1, 1] & d1InflxX < confirmWindow[, 1][length(confirmWindow[, 1])])]
      inflxConformX <- d1InflxX[which(d1InflxX > confirmWindow[1, 1]  & d1InflxX < confirmWindow[, 1][length(confirmWindow[, 1])])]             # Find inflection points within the window (defined as between peaks 1 and 2)
      threshold <- max(quantile(confirmWindow[, 2], probs=c(.95)), windowInflxY[windowPks[1]]/2)
      missedPks <- inflxConformX[which(inflxConformY > threshold)]                                                                              # Identify peaks above threshold (and greater than 50% of the height of the first peak)
      missedPks <- which(windowInflxX == missedPks[1])                                                                                          # Identify the inflection point in the original window which corresponds to the missed peak
      if(length(missedPks) > 0){windowPks[2] <- missedPks[1]}                                                                                   # Assign the missed peak as the second peak
    }

    if(length(windowPks) == 2 & windowInflxX[windowPks[2]] - windowInflxX[windowPks[1]] < (sr/10)){                                             # Rarely, the first peak contains two inflection points above threshold,
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
      if(windowPksY[1] > (m + (1.5*std)) & windowPksY[1] > pk_thrshd){                                                                          # inflection points, and that the second peak is greater than half the height of the first
        wX[1] <-  windowPks[1]                                                                                                                  # peak (unless the first peak is an artefact).
        if(windowPksY[2] > (m + (1.5*std)) &
           windowPksY[2] > (windowPksY[1]/2) | windowPksY[1] > (mean(deriv1) + (5*std(deriv1))) &
           windowPksY[2] < pk_thrshd){
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

  gt <- 1                                                                                                                                       # Define Global threshold factor for peaks (this may need altering depending on the variability in amplitude in the time series):

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
    #if(length(artefacts) > 0){                                                                                                                  # Recalculate mean and standard deviation after removing any artefacts identified
    #  remove <- c()
    #  for(j in artefacts[which(artefacts < length(wX))]){
    #    remove[j] <- which(abs(d1InflxX - wX[j]) == min(abs(d1InflxX - wX[j])))
    #  }
    #  remove <- remove[!is.na(remove)]
    #  newRem <- c()
    #  for(j in 1:length(remove)){
    #    newRem[(length(newRem)+1):(length(newRem)+21)] <- (remove[j] -10): (remove[j] + 10)                                                     # Remove inflection points around artefact peaks
    #  }
    #  if(sum(newRem < 1) > 0){
    #    newRem <- newRem[-(which(newRem < 1))]
    #  }
    #  m <- mean(d1InflxY[-newRem])
    #  std <- sd(d1InflxY[-newRem])
    #}

    while(length(windowPks) < 1){                                                                                                               # Each subsequent window will adjust (if required) until the next peak is detected

      if(windowExtnd > 10){
        windowStart <- 2
        windowExtnd <- 2.5
      }

      window[[i]] <- data.frame((wX[length(wX)] + windowStart*prevPkDist[i]):(wX[length(wX)] + windowExtnd*prevPkDist[i]),                      # Define the window initially as from (the previous peak + (peak to peak distance / 2))
                                deriv1[(wX[length(wX)] + windowStart*prevPkDist[i]):(wX[length(wX)] + windowExtnd*prevPkDist[i])])              # to (previous peak + (peak to peak distance * 1.35))
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
        if(max(windowPksY) > (m+(gt*std)) |  max(windowPksY) > (mean(deriv1[windowPks[1]:(windowPks[1]                                           # 1. a window does not include a genuine peak, and there are multiple secondary ones of similar height
                                                                                          + (3*prevPkDist[i]))]) + gt*std(deriv1[windowPks[1]:(windowPks[1] + (3*prevPkDist[i]))]))){                           # 2. there is an artefact with more than two peaks
          wX[length(wX)+1] <- windowPks[which(windowPksY == max(windowPksY))]                                           # 3. there are significantly large secondary peaks that also exceed the threshold for identification
          lowPks <- windowPksY[order(windowPksY)[1:2]]                                                                  # 4. a genuine peak has multiple inflection points
          if(lowPks[1] > m+(gt*std) & lowPks[2] > m+(gt*std) & (windowPks[3] - windowPks[1]) > (prevPkDist[i]/10)){
            # Check if the maximum peak is above a threshold relative to the time series (no 1)
            cat('\n','Potential artefact',  ', plot(', (wX[i-1]-100), ':', (wX[i-1]+300), ', deriv1[', (wX[i-1]-100), ':', (wX[i-1]+300),       # If it is, check if both lower peaks also exceed the threshold,
                '], type = "l") ,', 'wave', i, '+/- 2 removed because the non-max peaks were high')                                             # if so mark them as artefactual (no 2), if not assume they are secondary peaks (no 3)
            artefacts[length(artefacts) + c(1, 2, 3, 4, 5)] <- c(i-2, i-1, i, i+1, i+2)                                                         # An additional condition of the above is that the 'peaks' are not too close together so as to be inflection points (no 4)
          }
        }else{
          windowExtnd <- windowExtnd + 0.5                                                                                                      # If no 1 is the case, extend the window to look for peaks again
          windowPks <- c()
        }
      }

      if(length(windowPks) == 2){                                                                                                               # If two peaks are identified, they should be confirmed as genuine:
        if(max(windowPksY) > (m+(gt*std)) | max(windowPksY) > (mean(deriv1[windowPks[1]:(windowPks[1]                                           # Check if the maximum peak exceeds the global or local threshold
                                                                                         + (3*prevPkDist[i]))]) + gt*std(deriv1[windowPks[1]:(windowPks[1] + (3*prevPkDist[i]))]))){                                          # (local defined relative to peak to peak distance)
          if((windowPks[2] - windowPks[1]) < (prevPkDist[i]/3) & (windowPks[2] - windowPks[1]) > (prevPkDist[i]/10)){                           # If it does, check if the two peaks are close together in time
            wX[length(wX)+1] <- windowPks[which(windowPksY == max(windowPksY))]                                                                 # (within 1/3rd of the peak to peak distance, but not so close as to be a case of 4. (see above))
            cat('\n','Potential artefact',  ', plot(', (wX[i-1]-100), ':', (wX[i-1]+300),                                                       # If they are, mark the highest peak as artefactual
                ', deriv1[', (wX[i-1]-100), ':', (wX[i-1]+300), '], type = "l") ,', 'wave',
                i, '+/- 2 removed because two peaks found too close together')
          }else{                                                                                                                                # If they are not, go through both peaks in turn to confirm if genuine:
            if((windowPks[2] - windowPks[1]) < (prevPkDist[i]/10)){
              windowPks <- windowPks[which.max(windowPksY)]
              windowPksY <- windowPksY[which.max(windowPksY)]
            }else{
              if((windowPksY[1] > (m+(gt*std))) | windowPksY[1] > (mean(deriv1[windowPks[1]:(windowPks[1]                                         # Check if the first peak can be marked as genuine by:
                                                                                             + (3*prevPkDist[i]))]) + gt*std(deriv1[windowPks[1]:(windowPks[1] + (3*prevPkDist[i]))])) &                                        # 1. comparing to time series threshold
                 windowPksY[1] > (windowPksY[2]/2)){                                                                                                # 2. comparing to local threshold
                wX[length(wX)+1] <- windowPks[1]                                                                                                 # 3. comparing to height of the other peak (should exceed it's half maximum)
                if(windowPksY[2] > (m+(gt*std)) | windowPksY[2] > (mean(deriv1[windowPks[1]:(windowPks[1]                                        # If it can, the second peak can be marked as genuine with the same criteria
                                                                                             + (3*prevPkDist[i]))]) + gt*std(deriv1[windowPks[1]:(windowPks[1] + (3*prevPkDist[i]))])) &
                   windowPksY[2] > (windowPksY[1]/2)){
                  wX[length(wX)+1] <- windowPks[2]
                }
              }else{                                                                                                                              # If the first peak is not genuine, check the second and identify only the
                if(windowPksY[2] > (m+(gt*std)) | windowPksY[2] > (mean(deriv1[windowPks[1]:(windowPks[1] +                                       # second peak as genuine if appropriate.
                                                                                             (3*prevPkDist[i]))]) + gt*std(deriv1[windowPks[1]:(windowPks[1] + (3*prevPkDist[i]))])) &
                   windowPksY[2] > (windowPksY[1]/2)){
                  wX[length(wX)+1] <- windowPks[2]
                }
              }
            }
          }
        }else{                                                                                                                                  # If the maximum peak did not exceed threshold, extend the window and look again.
          windowExtnd <- windowExtnd + 0.5
          windowPks <- c()
        }
      }


      if(length(windowPks) == 1){                                                                                                               # If one peak is identified, confirm it is genuine by:
        if(windowPksY > (m+(gt*std)) | windowPksY > (mean(deriv1[windowPks[1]:(windowPks[1] +                                                    # 1. comparing to time series threshold
                                                                               (3*prevPkDist[i]))]) + gt*std(deriv1[windowPks[1]:(windowPks[1] + (3*prevPkDist[i]))])) |                                              # 2. comparing to local threshold
           windowPksY > (predict(d1p, wX[i-1])*0.9)){                                                                                             # 3. comparing it to the height of the previous peak
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

      ############# Window Artefact Identification: #############

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
      d1PeakSub <- CubicInterpSplineAsPiecePoly((round(wx[i])-(sr/4)):(round(wx[i])+(sr/4)),  # and search again
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
    o2 <- o2[-dup]
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
  ########################################################################################################################################
  # sep_beats segments the time series into individual pulse waveforms, allowing for an average morphology to be elucidated as well as
  # the later extraction of individual waveform features. Segmented beats undergo a cleaning process to exclude ones with aberrant
  # morphology (thresholds have been determined empirically). There is also an optional subsetting routine which allows for identification
  # of periods of increased heart rate (or lower IBI) in a time series and subsequent selection of the corresponding waveforms - for more
  # details on this see supplementary material.

  # Inputs:
  # odiff (vector of intervals between successive O (trough) points)
  # bc (the baseline corrected time series)
  # samp (sampling rate)
  # wuv (dataframe of w, u and v values for each identified beat)
  # wvlen (the median of O-diff +/- empirically determined value)
  # inx (vector of inflection points (x/time coordinates))
  # o (vector of O points)
  # ibi (vector of interbeat intervals (as determined by w points (peaks in 1st derivative)))
  # scale (logical, determines if pulse waveforms are scaled / normalized for amplitude)
  # q (logical, if TRUE number of beats excluded for each reason will be reported (useful for determining when thresholds may be to liberal/conservative))
  # subset (logical, determines if subsetting (as described above) is carried out)
  # boundaries (ISO analysis specific argument (for carrying over subsetting constants between time series))

  # Outputs:
  # average_wave (vector of values for the averaged waveform)
  # pulse (dataframe of all individual pulse waveforms in segmented form)
  # wuv (dataframe of w, u and v values (updated for excluded beats) for all waveforms)
  # rejects (list of vectors, each comrprised of waves excluded for a particular reason (see comments below))
  ########################################################################################################################################

  #################################### Subsetting #############################################

  if(subset == T){                                                                            # Optional subsetting (see function description)

    ab_ibi <- which(ibi > mean(ibi) + (4*sd(ibi)) | ibi < mean(ibi) - (4*sd(ibi)))            # Excluding points 4 sds above mean from the IBI time series
    ibi[ab_ibi] <- NA
    ibi <- ibi[!is.na(ibi)]
    meds <- rollmedian(ibi, k = 19)                                                           # Apply a rolling median to the IBI time series
    basemed <- mean(meds[1:50])                                                               # Baseline defined as the average of the first 50 rolling median points
    halfs <- (basemed - min(meds))/2                                                          # Find the amplitude between minimum and baseline
    post <- min(meds) + halfs

    meds_roi <- (which(meds == min(meds))[1] - 10) : (which(meds == min(meds))[1] + 10)       # For robustness in finding the minimum, ensure there are 10 beats either side of the minimum

    if( sum(which(meds_roi < 1)) > 0  | sum(which(meds_roi > length(ibi))) > 0 ){             # Minima without 10 beats either side trigger a warning
      message("Warning: No discernable minimum in 2 mcg trace")
      Sys.sleep(10)
      meds_roi <- (round(length(ibi)/2) - 10) : (round(length(ibi)/2) + 10)
    }

    indx <- meds_roi[which(ibi[meds_roi] == min(ibi[meds_roi]))]                              # The minimum is again identified
    if(indx < 15 | indx > (length(ibi) - 15)){                                                # and another check to ensure it is surrounded by adequate beats
      message("Warning: No discernable minimum in 2 mcg trace")
      Sys.sleep(10)
      meds_roi <- (round(length(ibi)/2) - 10) : (round(length(ibi)/2) + 10)
    }

    pre_indx <- which(abs(ibi[1:indx] - post) < 0.5) + 0                                      # Find the values that are closest to the half amplitude point before the minimum, within a certain range
    a <- 0.5
    while(length(pre_indx) < 1){                                                              # If there are no values within the initial y-axis range (0.5), extend the range to 1.5 and so on
      a <- a + 1
      pre_indx <- which(abs(ibi[1:indx] - post) < a) + 0
    }
    in2 <- which(abs(pre_indx - indx) == min(abs(pre_indx - indx)))
    pre_indx <- pre_indx[in2]

    post_indx <- indx + which(abs(ibi[indx:length(ibi)] - post) < 0.5)                        # Find the values that are closest to the half amplitude point after the minimum, same search process as for before the minimum
    a <- 0.5
    while(length(post_indx) < 1){
      a <- a + 0.25
      post_indx <- indx + which(abs(ibi[indx:length(ibi)] - post) < a)
    }
    while(post_indx[1] < indx + 15){
      post_indx <- post_indx[-1]
    }
    post_indx <- post_indx[1]

    rm(meds, basemed, halfs, post, indx, in2)
    ppg_pre <- (which((ppg[,1]) == beat[,1][pre_indx - 1]))                                   # Finding half-minimum timepoint in PPG time series (before minimum)
    ppg_post <- which((ppg[,1]) == beat[,1][post_indx + 1])                                   # Finding half-minimum timepoint in PPG time series (after minimum)

    if(length(ppg_pre) < 1 | length(ppg_post) < 1){                                           # If a half minimum point is not identified, trigger warning
      message("Subsetting failed: review of time series suggested")
      Sys.sleep(10)
      ppg_pre <- which((ppg[,1]) == beat[1,1])
      ppg_post <- which((ppg[,1]) == beat[nrow(beat),1])
    }

    subs <- which(wuv$wX > ppg_pre & wuv$wX < ppg_post)                                       # Identify beats within the subset
    wuv <- wuv[subs, ]
    ibi <- ibi[subs]
    cat(length(subs), "beats in subset (pre-cleaning)")
  }

  if(subset == "rep"){                                                                        # If previous subsetting parameters are being carried forward to a new time series, the following ensues
    subs <- which(wuv$wX > boundaries[1])
    subs_pre <- subs[1:boundaries[2]]
    tmp <- sum(is.na(subs_pre))                                                               # If there aren't enough beats remaining in the new time series, take as many as possible:
    if(tmp > 0){
      subs <- subs[1:(length(subs_pre)-tmp)]
      boundaries[2] <- length(subs)
    }else{
      subs <- subs[1:boundaries[2]]
    }
    wuv <- wuv[subs, ]
    ibi <- ibi[subs]
  }

  #################################### Segmentation ###########################################

  sourcedata <- bc[1:length(undetrended)]                                                     # Redefine baseline corrected data
  pulse <- data.frame(seq((-141/(samp*10)), ((wvlen*15 -9)-142)/(samp*10), by = 1/(samp*10))) # Define a dataframe to contain individual waves (with first column as the x-axis (in seconds)):

  afterO <- list()                                                                            # Lists for defining the number of values before and after the beginning and end of each waveform
  beforeO <- list()
  for(i in 1:length(wuv$wX)){

    splPolySub <- CubicInterpSplineAsPiecePoly((round(wuv$uX[i])-15):(round(wuv$uX[i]) +      # For each peak (w point) a rough segment around it is demarcated
                                                                        (wvlen+10)), sourcedata[(round(wuv$uX[i])-15):(round(wuv$uX[i]) +
                                                                                                                         (wvlen+10))], "natural")

    splSub <- predict(splPolySub, c(seq((wuv$uX[i]-14), (wuv$uX[i]+(wvlen+4)), 0.1)))         # Segment is discretized
    splSub <-  as.data.frame(splSub)
    splSub <- cbind(splSub, c(seq((wuv$uX[i]-14), (wuv$uX[i]+(wvlen+4)), 0.1)))               # Each segment is converted to a dataframe where time is relative to W
    colnames(splSub) <- c('y', 'x')
    if(scale == TRUE){                                                                        # Scaling amplitude such that v-u = 1
      splSub$y <- splSub$y/(wuv$diffVU[i])
    }
    yDiff <- splSub$y[141]                                                                    # Adjust such that u = 0, v = 1 (on y-axis)
    splSub$y <- splSub$y - yDiff

    if(scale == TRUE){
      splPolySub2 <- CubicInterpSplineAsPiecePoly(splSub$x, splSub$y, "natural")              # Find the x-value for each wave that corresponds to when its amplitude is 0.5
      halfCross <- solve(splPolySub2, b = 0.5, deriv = 0)
      halfCross <- halfCross[which(abs(halfCross - wuv$wX[i]) ==                              # Waves are aligned relative to this value (empirically this method has lead to better alignments than aligning by w alone)
                                     min(abs(halfCross - wuv$wX[i])))]
      if(halfCross-14 < splPolySub$knots[1] |
         halfCross+(wvlen+4) > splPolySub$knots[length(splPolySub$knots)]){
        splSub2 <- predict(splPolySub,
                           c(seq((splPolySub$knots[1]),
                                 (splPolySub$knots[length(splPolySub$knots)]), 0.1)))
        splSub2 <- as.data.frame(splSub2)
        splSub2 <- cbind(splSub2, c(seq((splPolySub$knots[1]),
                                        (splPolySub$knots[length(splPolySub$knots)]), 0.1)))
      }else{
        splSub2 <- predict(splPolySub,
                           c(seq((halfCross-14), (halfCross+(wvlen+4)), 0.1)))
        splSub2 <- as.data.frame(splSub2)
        splSub2 <- cbind(splSub2, c(seq((halfCross-14), (halfCross+(wvlen+4)), 0.1)))
      }
      colnames(splSub2) <- c('y', 'x')

      splSub2$y <- splSub2$y/(wuv$diffVU[i])                                                  # Scale again, and again adjust y-axis such that u = 0, v = 1
      yDiff <- wuv$uY[i] / wuv$diffVU[i]
      splSub2$y <- splSub2$y - yDiff
    }else{
      splSub2 <- splSub
    }

    afterO[[i]] <- which(splSub2$x > inx[o][min(which(inx[o] > wuv$wX[i]))])                  # Find values in the segment after the end of the waveform
    beforeO[[i]] <- which(splSub2$x < inx[o][max(which(inx[o] < wuv$wX[i]))])                 # Find values in the segment before the beginning of the waveform

    splSub3 <- c()                                                                            # Correct dataframe of segment
    for(j in 1:nrow(splSub2)){
      splSub3[j+1] <- splSub2$y[j]
    }

    if(length(splSub3) == nrow(pulse)){                                                       # Correct dataframe of segment to ensure is has the same row length as the main (pulse) dataframe
      if(length(afterO[[i]]) > 0){
        diff2 <- abs(length(splSub3) - max(afterO[[i]]))
        for(j in 1:diff2){
          afterO[[i]] <- c(afterO[[i]], (max(afterO[[i]]) + 1))
        }
      }
    }
    if(length(splSub3) > nrow(pulse)){
      diff <- length(splSub3) - nrow(pulse)
      len <- length(splSub3)
      splSub3 <- splSub3[-((len - (diff-1)):len)]
      if(diff > 1){
        if(length(afterO[[i]]) > 1 ){
          afterO[[i]] <- afterO[[i]][-(which(afterO[[i]] > length(splSub3)))]
        }else{
          afterO[[i]] <- afterO[[i]][-(which(afterO[[i]] > length(splSub3)))]
        }
      }
    }
    if(length(splSub3) < nrow(pulse)){                                                        # If there are too few rows in the segment, NA values are added to match that of pulse
      diff <- nrow(pulse) - length(splSub3)
      splSub3 <- c(splSub3, rep(NA, diff))
      if(length(afterO[[i]]) > 0){
        diff2 <- length(splSub3) - max(afterO[[i]])
        for(j in 1:diff2){
          afterO[[i]] <- c(afterO[[i]], (max(afterO[[i]]) + 1))
        }
      }
    }

    pulse <- cbind(pulse, splSub3)                                                            # Combine each segment as an additional column to the main dataframe
  }

  for(i in 1:(ncol(pulse) -1)){                                                               # Assign column names
    colnames(pulse)[i+1] <- paste("wave", i, sep = "_")
  }
  colnames(pulse)[1] <- "x"

  for(i in 2:(ncol(pulse))){                                                                  # For each segment, replace any values outside of the waveform with NA values
    pulse[, i][afterO[[(i-1)]][-1]] <- NA
  }
  for(i in 2:(ncol(pulse))){
    pulse[, i][beforeO[[(i-1)]][-1]] <- NA
  }

  wavelengths <- c()                                                                          # Find length of each waveform in pulse dataframe (excluding NA values)
  for(i in 2:ncol(pulse)){
    wavelengths[i] <- length(pulse[, i][!is.na(pulse[, i])])
  }
  wavelengths <- wavelengths[!is.na(wavelengths)]


  ####################################### Cleaning ############################################

  extra_long_wave <- c()
  for(i in 2:length(wavelengths)){                                                            # Identify waveforms that are excessively long (greater than the mean + sd wavelength,
    if(wavelengths[i] > (mean(wavelengths) +                                                  # and greater than 1.8 times the length of the previous waveform)
                         sd(wavelengths)) & wavelengths[i] > 1.8*(wavelengths[i-1])){
      extra_long_wave[i] <- i
    }
  }
  extra_long_wave <- extra_long_wave[!is.na(extra_long_wave)]
  if(length(extra_long_wave) > 0){                                                            # Exclude identified waveforms, report rejection with option to plot if q = TRUE
    cat("\n", length(extra_long_wave), "/", (ncol(pulse)-1),
        "waves removed for being abnormally long")
    if(q == TRUE){
      plotyyy <- 0
      while(plotyyy == 0){
        plotyy <- readline(prompt = "Would you like to view? (enter yes or no)")
        if(plotyy == "yes"){
          for(i in 1:length(extra_long_wave)){
            plot((wuv$wX[extra_long_wave[i]]-samp*2):(wuv$wX[extra_long_wave[i]]+samp*2),
                 bc[(wuv$wX[extra_long_wave[i]]-samp*2):(wuv$wX[extra_long_wave[i]]+samp*2)],
                 type = "l")
            points(wuv$wX[extra_long_wave[i]],
                   wuv$wY[extra_long_wave[i]], pch =19, col = 'red')
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


  extra_short_waves <- c()                                                                    # Identify waveforms that are excessively short in length (less than the mean - sd*3 wavelength)
  for(i in 2:length(wavelengths)){
    if(wavelengths[i] < mean(wavelengths[-1]) - 3*sd(wavelengths[-1])){
      extra_short_waves[i] <- i
    }
  }
  extra_short_waves <- extra_short_waves[!is.na(extra_short_waves)]                           # Exclude identified waveforms, report rejection with option to plot if q = TRUE
  if(length(extra_short_waves) > 0){
    cat("\n", length(extra_short_waves), "/", (ncol(pulse)-1),
        "waves removed for being abnormally short in duration")
    if(q == TRUE){
      plotyyy <- 0
      while(plotyyy == 0){
        plotyy <- readline(prompt = "Would you like to view? (enter yes or no)")
        if(plotyy == "yes"){
          for(i in 1:length(extra_short_waves)){
            plot((wuv$wX[extra_short_waves[i]]-samp*2):(wuv$wX[extra_short_waves[i]]+samp*2),
                 bc[(wuv$wX[extra_short_waves[i]]-samp*2):(wuv$wX[extra_short_waves[i]]+samp*2)],
                 type = "l")
            points(wuv$wX[extra_short_waves[i]], wuv$wY[extra_short_waves[i]],
                   pch =19, col = 'red')
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


  double_segments <- c()                                                                      # Identify waveforms that are likely to have two systolic peaks
  double_peaks_wave <- list()                                                                 # If more than one peak is identified with a substantial time delay, assume two systolic peaks
  for(i in 2:ncol(pulse)){
    # Make a spline to identify inflection points:
    wave <- pulse[, i][!is.na(pulse[, i])]
    sfunction <- splinefun(1:length(wave), wave, method = "natural")
    spline1 <- sfunction(seq(1, length(wave)), deriv = 0)
    splinePoly <- CubicInterpSplineAsPiecePoly(1:length(wave), spline1, "natural")
    inflexX <- solve(splinePoly, b = 0, deriv = 1)
    inflexY <- predict(splinePoly, inflexX)
    peaks <- which(inflexY > mean(inflexY) + 1.7*sd(inflexY))
    if(length(peaks) > 0){
      if( length(peaks) > 1 & ((inflexX[peaks][length(inflexX[peaks])] - inflexX[peaks[1]]) > 100) ){
        double_segments[i] <- i
        double_peaks_wave[[i]] <- wave
      }
    }
  }
  double_segments <- double_segments[!is.na(double_segments)]
  if(length(double_peaks_wave) > 1){double_peaks_wave <-                                      # Exclude identified waveforms, report rejection with option to plot if q = TRUE
    double_peaks_wave[-(which(sapply(double_peaks_wave, is.null)))]}
  if(length(double_segments) > 0){
    cat("\n", length(double_segments), "/", (ncol(pulse)-1),
        "waves removed for containing two systolic peaks")
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


  systolic_endings <- c()                                                                     # Identify waveforms that enter a second systolic upstroke
  for(i in 2:ncol(pulse)){
    wave <- pulse[, i]
    wave <- wave[!is.na(wave)]
    if(wave[length(wave)] > 0.25 | (length(wave) > mean(wavelengths) &
                                    (max(wave[round(0.75*length(wave)):length(wave)]) > 0.8))){                     # If the wave is above average length and the last quarter has a value above 0.8, consider it a second systolic peak
      systolic_endings[i] <- i
    }
  }
  systolic_endings <- systolic_endings[!is.na(systolic_endings)]
  if(length(systolic_endings) > 0){                                                           # Exclude identified waveforms, report rejection with option to plot if q = TRUE
    cat("\n", length(systolic_endings), "/", (ncol(pulse)-1),
        "waves removed for including the following systolic upstroke")
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

  average_wave <- find_average(p = pulse, ao = afterO)                                        # Calculate average wave

  drops_below_o <- c()                                                                        # Identify waveforms that fall significantly below O, relative to the average wave
  for(i in 2:ncol(pulse)){
    wave <- pulse[, i][!is.na(pulse[, i])]
    thd <- 4
    if((min(wave) < wave[1]*thd | min(wave) < -0.3 | min(wave) <
        min(average_wave[!is.na(average_wave)])*thd) & min(wave) < 0){
      drops_below_o[i] <- i
    }
  }
  drops_below_o <- drops_below_o[!is.na(drops_below_o)]                                       # Exclude identified waveforms, report rejection with option to plot if q = TRUE
  if(length(drops_below_o) > 0){
    cat("\n", length(drops_below_o), "/", (ncol(pulse)-1),
        "waves removed for having values significantly below baseline")
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
    wuv <- wuv[-(drops_below_o-1), ]
    pulse <- pulse[, -drops_below_o]
    ibi <- ibi[-(drops_below_o-1)]
  }

  average_wave <- find_average(p = pulse, ao = afterO)                                        # Recalculate the average after removing any waves in the previous section
  sd_wave <- find_sd(p = pulse, ao = afterO)                                                  # Find the SD wave

  resid_sd <- c()                                                                             # Identify waveforms that have a high SD of residuals
  for(i in 2:ncol(pulse)){
    residuals <- average_wave[142:length(average_wave)] - pulse[, i][142:length(pulse[, i])]
    resid_sd[i] <- sd(residuals[!is.na(residuals)][-c(1:100)])                                # The first 100 values are excluded from this assessment, since at this point
  }                                                                                           # all values are close to the average on the systolic upstroke
  resid_sd <- sqrt(resid_sd)
  resid_sd <- resid_sd[-1]
  thld <-  0.35                                                                               # threshold empirically determined
  hrsd_waves <- which(resid_sd > thld) + 1
  if(length(hrsd_waves) > 0){                                                                 # Exclude identified waveforms, report rejection with option to plot if q = TRUE
    cat("\n", length(hrsd_waves), "/", (ncol(pulse)-1),
        "waves removed for having high residual SD")
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


  outlier_waves <- c()                                                                        # Identify waveforms that are beyond a certain number of SDs from the average
  for(i in 2:ncol(pulse)){
    breaches <- c()
    wave <- pulse[, i]
    for(j in 142:length(wave)){   # look only after w
      if(is.na(pulse[j, i]) == FALSE & is.na(average_wave[j] + 4*sd_wave[j]) == FALSE){
        if(pulse[j, i] > (average_wave[j] + 4*sd_wave[j]) |
           pulse[j, i] < (average_wave[j] - 4*sd_wave[j])){
          breaches[j] <- 1
        }
      }
    }
    if(sum(breaches[!is.na(breaches)]) > 0){
      outlier_waves[i] <- i
    }
  }
  outlier_waves <- outlier_waves[!is.na(outlier_waves)]                                       # Exclude identified waveforms, report rejection with option to plot if q = TRUE
  if(length(outlier_waves) > 0){
    cat("\n", length(outlier_waves), "/", (ncol(pulse)-1),
        "waves removed for having values beyond 5 SDs from the average")
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

  average_wave <- find_average(p = pulse, ao = afterO)                                        # Final recalculation of the average waveform

  cat(nrow(wuv), "beats carried over for 2mg subsetting (post-cleaning)")                     # Print number of beats post subsetting and cleaning (only informative for ISO analysis)

  rejects <- list(extra_long_wave, extra_short_waves, double_segments,                        # Save rejected waves
                  systolic_endings, drops_below_o, hrsd_waves, outlier_waves)

  dat <- list(average_wave, pulse, wuv, rejects)
  if(subset == T){dat <- list(average_wave, pulse, wuv, rejects, c(ppg_pre, nrow(wuv)), ibi)}

  return(dat)
}




find_average <- function(p, ao){
  ########################################################################################################################################
  # find_average generates an average (mean) wave from a sample of segmented waveforms. Due to variance in waveform length, the tail
  # section of the average wave is less robustly calculated, and relies on a degree of inference as to the likely trajectories of
  # waveforms that terminate before the average length is reached.

  # Inputs:
  # p (dataframe consisting of all segmented waveforms)
  # ao (a vector indicating the end point of each waveform)

  # Outputs:
  # avWave (a vector consisting of values for the average waveform)
  ########################################################################################################################################

  last10 <- list()
  for(i in 2:(ncol(p))){
    last10[[i-1]] <- p[, i][!is.na(p[, i])][(length(p[, i][!is.na(p[, i])])-10):(length(p[, i][!is.na(p[, i])]))]      # Find the last 10 values of each wave
  }

  avLast10 <- c()
  for(i in 1:10){
    rowVec <- c()
    for(j in 1:length(last10)){
      rowVec[j] <- last10[[c(j, i)]]
    }
    avLast10[i] <- mean(rowVec[!is.na(rowVec)])                                                                        # Find the average across waveforms for each of these 10 values
  }

  avDiff <- c()
  for(i in 1:9){                                                                                                       # Find the average gradient of tail decay based on these averages
    avDiff[i] <- avLast10[i+1] - avLast10[i]
  }

  p2 <- p
  for(i in 2:(ncol(p2))){
    for(j in 1:length(p2[, i][ao[[(i-1)]][-1]])){                                                                      # Generate a clone dataframe of waveforms, where trajectories of waves that end before
      p2[, i][ao[[(i-1)]][-1]][j] <- p2[, i][ao[[(i-1)]][1]] + j*mean(avDiff[1:5])                                     # the average wavelength are continued according to the average tail gradient calculated above.
    }
  }

  avWav <- c()                                                                                                         # Average wave is generated by calculating the mean value for each row of the dataframe.
  sdWav <- c()                                                                                                         # Note that the standard deviation and median can also be calculated for each row.
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


  end <- c()                                                                                                           # Determine the end of the average waveform as the average (mode) y-value of terminations points on individual waveforms.
  for(i in 2:ncol(p)){
    end[i-1] <- p2[, i][ao[[(i-1)]][1]]
  }

  mode <- function(x){                                                                                                 # define mode function
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }

  if(length(end) > 1){
    avEnd <- mode(round(end[!is.na(end)], digits = 2))
    belowO <- which(avWav < avEnd)                                                                                     # Remove elements of the average waveform after where its end should be
    if(sum(which(belowO > (length(avWav)/(4/3)))) > 0){                                                                # If the average waveform drops below baseline before it's specified endpoint, consider this the termination point instead.
      b. <- belowO[min(which(belowO > (length(avWav)/(4/3))))]                                                         # This conditional is limited to drops below baseline within the final 25% of the average waveform (to prevent premature terminations).
    }else{
      b. <- NA
    }
    if(is.na(b.) == FALSE){
      avWav[b.:length(avWav)] <- NA
    }
  }

  first <- c()                                                                                                         # Define the starting point of the average waveform, determined by the median x-value of all other waveforms' first elements
  for(i in 2:ncol(p)){
    stpqw <- p[, i][1:200]
    stq <- stpqw[is.na(stpqw)]
    first[i-1] <- length(stq) + 1
  }
  start <- median(first)-2
  avWav[1:start] <- NA

  avWave <- as.vector(avWav)                                                                                           # Ensure continuity of baseline in the average waveform - if the previously calculated termination point falls below the
  avWave_tmp <- avWave[!is.na(avWave)]                                                                                 # y-value of the starting point, shorten the wave.
  while(avWave_tmp[length(avWave_tmp)] < avWave_tmp[1]){
    avWave_tmp <- avWave_tmp[-length(avWave_tmp)]
    avWave[length(avWave_tmp) + which(!is.na(avWave))[1]] <- NA
  }

  return(avWave)
}



find_sd <- function(p, ao){
  ########################################################################################################################################
  # find_sd is a clone of find_average; only the output is changed (from mean to SD)

  # Inputs:
  # p (dataframe consisting of all segmented waveforms)
  # ao (a vector indicating the end point of each waveform)

  # Outputs:
  # sd_wave (a vector consisting of standard deviation values (can also be plotted to indicate variance))
  ########################################################################################################################################

  last10 <- list()
  for(i in 2:(ncol(p))){
    last10[[i-1]] <- p[, i][!is.na(p[, i])][(length(p[, i][!is.na(p[, i])])-10):(length(p[, i][!is.na(p[, i])]))]
  }

  mean_last10 <- c()
  for(i in 1:10){
    row_vector <- c()
    for(j in 1:length(last10)){
      row_vector[j] <- last10[[c(j, i)]]
    }
    mean_last10[i] <- mean(row_vector[!is.na(row_vector)])
  }

  mean_diff_last10 <- c()
  for(i in 1:9){
    mean_diff_last10[i] <- mean_last10[i+1] - mean_last10[i]
  }

  p_for_finding_average <- p
  for(i in 2:(ncol(p_for_finding_average))){
    for(j in 1:length(p_for_finding_average[, i][ao[[(i-1)]][-1]])){
      p_for_finding_average[, i][ao[[(i-1)]][-1]][j] <- p_for_finding_average[, i][ao[[(i-1)]][1]] + j*mean(mean_diff_last10[1:5])
    }
  }

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
  ########################################################################################################################################
  # diast_pk finds the position of the diastolic peak (D) on the average wave so as to inform identification of OSND points on all
  # individual waveforms. Since the average waveform being evaluated is also free from NA values, the resultant shift along the
  # x-axis in terms of index points is also calculated.

  # Inputs:
  # avw (average waveform)
  # sr (sampling rate)
  # scale (logical - to establish bounds of y-axis)
  # dias_param (optional - to inform search for D if the waveform is a HED modelled one)

  # Outputs:
  # diastPk (diastolic peak position on the average waveform)
  # xShift (difference in x-axis values of the average waveform due to removal of NA values)
  ########################################################################################################################################

  avw <- avw[!is.na(avw)]                                                                                    # Remove NA values from the average wave

  if(scale == TRUE){                                                                                         # Find the new W peak position (0.5) after removing NAs
    xShift <- which(abs(avw-0.5) == min(abs(avw - 0.5)))
  }else{
    xShift <- which.min(abs(avw))
  }
  avWavePoly <- CubicInterpSplineAsPiecePoly(1:length(avw), avw, "natural")
  avInflexX <- solve(avWavePoly, b = 0, deriv = 1)
  avInflexY <- predict(avWavePoly, avInflexX)

  # plot(avWavePoly)                                                                                         # Optional plots (largely for debugging purposes)
  # points(avInflexX, avInflexY)

  peaks_above_0 <- which(avInflexY > 0)                                                                      # Find any peaks above baseline

  if(!is.null(dias_param)){                                                                                  # This section is for fitted waves only:
    if(avInflexY[peaks_above_0[1]] <  mean(avw)){                                                            # Ensure peaks are above the mean y-value of the waveform
      peaks_above_0 <- peaks_above_0[-1]
    }
    dias_param <- dias_param + avInflexX[peaks_above_0[1]]                                                   # Since model parameters for D will be relative to S, find the objective value by adding S
    diastPk <-  avInflexX[peaks_above_0[which( abs(avInflexX[peaks_above_0] - dias_param) ==                 # The diastolic peak is labelled as the the second peak following S
                                                 min(abs(avInflexX[peaks_above_0] - dias_param)))]]
    if(avInflexY[peaks_above_0[which( abs(avInflexX[peaks_above_0] - dias_param) ==
                                      min(abs(avInflexX[peaks_above_0] - dias_param)))]] - max(avw) < 2){    # If the value assigned is just the max of the wave (i.e systolic peak, default to sr*500):
      diastPk <- 5*sr
    }
    peaks <- peaks_above_0
    fitted <- 1
  }else{
    fitted <- 0
  }

  if(fitted != 1){                                                                                           # This section is for non-fitted waveforms (i.e original data)
    peaks <- order(avInflexY[peaks_above_0], decreasing = TRUE)                                              # Order peaks above 0 in descending order
    peaks <- peaks_above_0[peaks]                                                                            # Ensures peaks refers to the correct avInflexX points
    dif <- peaks[1] - 1                                                                                      # Define the first peak as peak 1
    peaks <- peaks - dif
    if(dif > 0){                                                                                             # Correct avInflexX / avInflexY as a result of the peak assignment
      avInflexX <- avInflexX[-c(1:dif)]
      avInflexY <- avInflexY[-c(1:dif)]
    }
    if(sum(peaks < 1) > 0){                                                                                  # Remove 0 or negative values if they have resulted from the above adjustment
      peaks <- peaks[-which(peaks < 1)]
    }

    if(peaks[1] != 1 &
       length(peaks) > 1){diastPk <- avInflexX[peaks[1]]}else{diastPk <- avInflexX[peaks[2]]}                # Identify the second highest peak as the diastolic peak, unless it is unreasonably
    if(!is.na(diastPk) & diastPk == avInflexX[peaks[2]]){                                                    # high (i.e a diastolic to systolic peak ratio of > 0.75)
      ai. <- avInflexY[peaks[2]] / avInflexY[peaks[1]]                                                       # (the threshold for unreasonably high must be dependent on how close the peak is
      sd_time <- avInflexX[peaks[2]] - avInflexX[peaks[1]]                                                   # to S on the x-axis (lower threshold for closer))

      while( (ai. > 0.45 & sd_time < 75) | ai. > 0.7 | ai. < 0.05 ){                                         # Move forward through the peaks until finding one that is not unreasonably high (or unreasonably low)
        old <- peaks[2]                                                                                      # Assign the previous Dpeak as old, and remove it from the peaks list
        peaks <- peaks[-2]
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

        if(is.na(ai.)){                                                                                      # If there are no valid peaks, choose the default estimation (5 * sampling rate)
          diastPk <- 5*sr
          ai. <- 0.5
          sd_time <- 100
        }
      }
    }
  }

  if(is.na(diastPk) | diastPk < avInflexX[peaks[1]]){                                                        # diastPk will be NA for class 3 waveforms, in which case set the default estimation also
    diastPk <- 5*sr
  }
  return(c(diastPk, xShift))
}



osnd_of_average <- function(aw, dp, diff, sr, plot = TRUE){
  ########################################################################################################################################
  # osnd_of_average identifies OSND points on the average waveform and individual waves.

  # Inputs:
  # aw (average waveform, or an individual waveform)
  # dp (diastolic peak (as established by the diast_pk function))
  # diff (difference in x-axis values of the average waveform due to removal of NA values (outputted from diast_pk as xShift))
  # sr (sampling rate)
  # plot (logical - for plotting the inputted waveform with identified OSND points)

  # Outputs:
  # osnd (OSND values for the inputted waveform)
  ########################################################################################################################################

  switch <- 0                                                                                                # Switch is defined according to whether a D peak is found (switch <- 1) or not found (switch <- 0) in the waveform
  aw <- aw[!is.na(aw)]

  avWavPoly <- CubicInterpSplineAsPiecePoly(1:length(aw), aw, "natural")                                     # Convert inputted waveform to polynomial spline
  sfunction <- splinefun(1:length(aw), aw, method = "natural")
  d1Wav <- sfunction(1:length(aw), deriv = 1)
  d1WavPoly <- CubicInterpSplineAsPiecePoly(1:length(aw), d1Wav, "natural")                                  # Create spline of first derivative

  # plot(avWavPoly)

  d1InflxX <- solve(d1WavPoly, b = 0, deriv = 1)                                                             # Identify inflection points on first derivative
  d1InflxY <- predict(d1WavPoly, d1InflxX)

  # plot(d1WavPoly)
  # points(d1InflxX, d1InflxY)

  wavInflxX <- solve(avWavPoly, b = 0, deriv = 1)                                                            # Identify inflection on original waveform
  wavInflxY <- predict(avWavPoly, wavInflxX)

  notchRange <- which(d1InflxX > (3.104572*sr - diff) & d1InflxX < dp)                                       # Identify a range on the x-axis within which the notch is likely to be located


  if(length(notchRange) < 1 | (length(notchRange) == 1 & d1InflxY[notchRange][1] < -0.02)){                  # If there is no inflexion point detected within the notch range,
    new.n <- ((3.104572*sr) + dp)/2                                                                          # this could be because there is a plateu rather than a peak.
    # In this case, taking the mean value of the notch range boundaries gives a reasonable approximation
    if(new.n > length(aw) | length(aw) - new.n < (5*(length(aw)/100)) ){                                     # In case new.n is greater than the number of datapoints (or very close to the end),
      new.n <- d1InflxX[length(d1InflxX)]                                                                    # set D to last inflection point on derivative (assuming the wave is very short)
    }
  }else{
    if(length(notchRange) != 1){
      a. <- which(d1InflxY[notchRange] == max(d1InflxY[notchRange]))
      while(d1InflxX[notchRange[a.]] < sr*3){                                                                # In cases where the renal peak is higher on 1st deriv than the notch peak,
        b. <- 2                                                                                              # ensure the notch peak is limited by x-axis
        a. <- order(d1InflxY[notchRange], decreasing = TRUE)[b.]
        b. <- 3
      }
      if(which(d1InflxY == max(d1InflxY)) == notchRange[a.]){                                                # Ensure the 1st peak is not the notch
        if(length(notchRange[which(notchRange > notchRange[a.])]) > 0){
          notchRange <- notchRange[which(notchRange > notchRange[a.])]
        }
      }
      if(is.na(d1InflxX[notchRange[a.]])){
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

  if(sum(which(wavInflxX < new.n)) < 1){                                                                     # If there are no inflection points before new.n, assume it is
    if (sum(which(wavInflxX > new.n & wavInflxY < new.ny)) > 1){                                             # incorrect and move to the next inflection point that has a lower y-value than it
      new.n <- wavInflxX[min(which(wavInflxX > new.n & wavInflxY < new.ny))]                                 # If, however, there are no inflection points that fit this criterion, keep N as it is.
      new.ny <- predict(avWavPoly, new.n)
    }
  }

  if (sum(which(wavInflxX < new.n)) > 0){
    if(length(wavInflxX) > 1 & wavInflxY[max(which(wavInflxX < new.n))] < new.ny){                           # After having found the notch on all waves, see if there are inflexion points either side
      new.n <- wavInflxX[max(which(wavInflxX < new.n))]                                                      # If inflexion point before is lower and inflexion point after is higher (on y axis), this must mean a second peak aka a canonical / class 1 wave
      new.ny <-  predict(avWavPoly, new.n)                                                                   # Thus if this criterion is fulfilled, consider N and D as distinct.

      if(sum(wavInflxX > new.n) > 0){                                                                        # Take the inflection point after the notch as the D peak.
        d. <- wavInflxX[min(which(wavInflxX > new.n))]                                                       # If there is no inflection point after the notch, take instead the next inflection point on deriv1
        d.y <- predict(avWavPoly, d.)
      }else{
        d. <- d1InflxX[min(which(d1InflxX > new.n))]
        d.y <- predict(avWavPoly, d.)
      }
      switch <- 1
    }
  }


  if(plot == TRUE){
    plot(avWavPoly)
    points(wavInflxX, wavInflxY)
    points(new.n, new.ny, col = "red")
    if(switch == 1){
      points(d., d.y, col = "blue")
    }
  }

  w. <- d1InflxX[which( d1InflxY ==  max(d1InflxY) & d1InflxX[which(d1InflxY ==  max(d1InflxY))] < new.n)]   # Find W: the max inflection point on first deriv
  if(length(w.) < 1){
    w. <- 1
  }
  w.y <- predict(avWavPoly, w.)
  if(plot ==TRUE){points(w., w.y)}

  hhaw <- max(d1InflxY)/2                                                                                    # Find U and V, by finding half the height of w (on derivative y-axis)
  halfHeights <- solve(d1WavPoly, b = hhaw)

  if(length(halfHeights) < 2){                                                                               # If only one half height value is detected, assume the segment
    halfHeights[2] <- halfHeights[1]                                                                         # had quite a high O value such that O was higher than U
    halfHeights[1] <- solve(d1WavPoly, b = d1Wav[1])[1]
  }

  if(length(halfHeights) > 2){                                                                               # If more than two points at half height are detected:
    a <- halfHeights - w.                                                                                    # Find the distance between each detected point and w and consider the two that are equidistant to w as u and v
    postW <- which(a > 0)
    preW <- which(a < 0)
    if(length(preW) < 1){                                                                                    # If a point is not found before W despite >2 points, take the first value as U
      halfHeights <- c(1, min(halfHeights[postW]))
    }else{
      halfHeights <- c(max(halfHeights[preW]), min(halfHeights[postW]))
    }
  }
  u <- halfHeights[1]
  v <- halfHeights[2]

  uvY <- predict(avWavPoly, halfHeights)                                                                     # Find the corresponding y-values for U and V on the original wave
  uY <- uvY[1]
  vY <- uvY[2]
  if(plot == TRUE){
    points(u, uY)
    points(v, vY)
  }

  if(length(which(wavInflxX < w.)) > 0){                                                                     # Find O:
    o. <- wavInflxX[max(which(wavInflxX < w.))]                                                              # If there are inflection points before w, use the maximum one as 0
  }else{
    if(length(which(d1InflxX < w.)) < 1){                                                                    # If there are no inflection points before w, check if there are any in the first derivative before w
      o. <- 1                                                                                                # If not, assign O to the first value
    }else{                                                                                                   # If there are, make sure the max inflection point is not greater than 0 on the original y-axis (which would be higher than U)
      if((predict(avWavPoly, d1InflxX[max(which(d1InflxX < w.))]) > 0)){
        o. <- 1                                                                                              # If the max inflection point is above 0, assign O to the first value
      }else{
        o. <- d1InflxX[max(which(d1InflxX < w.))]                                                            # If the max inflection point is not above 0, define it as O
      }
    }
  }
  o._yval <- predict(avWavPoly, o.)
  if(plot == TRUE){points(o., o._yval, pch = 19)}

  if( (v - w.) > 100 ){                                                                                      # Find S:
    s2Y <- max(wavInflxY)                                                                                    # First define S as the point on the waveform that is twice the distance of w to v (x-axis)
    s2 <- wavInflxX[which(wavInflxY == max(wavInflxY))]                                                      # (if v is erroneous (making the (v - w.) > 100 statement true), then the maximum point of the waveform is chosen instead)
  }else{
    s2 <- w. + 2*(abs(v - w.))
    s2Y <- predict(avWavPoly, s2)
  }
  # points(s2, s2Y)
  s1Y <- max(wavInflxY)                                                                                      # Define an alternative S as the maximum point on the waveform
  s1 <- wavInflxX[which(wavInflxY == max(wavInflxY))]
  if((s1 - w.) < (s2 - w.)){                                                                                 # Define S definitvely as the closer of the above two points to w
    s. <- s1                                                                                                 # (This routine is necessary in class 3 and above waveforms, where the maximum of the waveform does
    s.y <- s1Y                                                                                               # not correspond to the maximum of the systolic peak, or when the systolic peak is more of a shoulder)
  }else{
    s. <- s2
    s.y <- s2Y
  }
  if(plot == TRUE){points(s., s.y, pch = 19)}

  if(switch == 0){                                                                                           # Define the final OSND values to be outputted
    x <- c(o., s., new.n, new.n)                                                                             # If no D was found, consider D to have the same value as N
  }else{                                                                                                     # (This is the equivalent to defining an 'inflection point' in diastole and is common in higher class waveforms)
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
  ####################################################################################################
  # feature_extract extracts descriptive morphological features for each waveform inputted.

  # Inputs:
  # oa (a list of osnd values for each wave (x and y coordinates))
  # p (dataframe of all individual waveforms (discrete form))
  # pw (list of all individual waveforms (polynomial spline form))

  # Outputs:
  # features (a dataframe of all extracted features for each waveform)
  ####################################################################################################

  s_vals <- c()                                                                          # Systolic peak values
  for(i in 1:length(oa)){
    s_vals[i] <- oa[[i]]$y[2]
  }

  n_vals <- c()                                                                          # Notch values
  for(i in 1:length(oa)){
    n_vals[i] <- oa[[i]]$y[3]
  }

  d_vals <- c()                                                                          # Diastolic peak values
  for(i in 1:length(oa)){
    d_vals[i] <- oa[[i]]$y[4]
  }

  np_ratio <- c()                                                                        # Notch to peak ratio
  for(i in 1:length(oa)){
    np_ratio[i] <- oa[[i]]$y[3] / oa[[i]]$y[2]
  }

  ppt <- c()                                                                             # Peak to peak time (systolic to diastolic)
  for(i in 1:length(oa)){
    ppt[i] <- oa[[i]]$x[4] - oa[[i]]$x[2]
  }

  max_amp <- c()                                                                         # Maximum amplitude
  for(i in 1:length(oa)){
    inx <- solve(pw[[i]], b = 0, deriv = 1)
    iny <- predict(pw[[i]], inx)
    max_amp[i] <- max(iny)
  }

  auc <- c()                                                                             # Total area under waveform
  for(i in 1:length(oa)){
    wave <- p[, (i+1)]
    v <- which(!is.na(wave))
    auc[i] <- AUC(x = v, y = wave[v], method = "spline")
  }

  auc_s <- c()                                                                           # Area under the waveform after systolic peak
  for(i in 1:length(oa)){
    wave <- p[, (i+1)]
    wave <- wave[!is.na(wave)]
    s <- oa[[i]]$x[2]
    v <- 1:length(wave)
    v <- which(v > s)
    auc_s[i] <- AUC(x = v, y = wave[v], method = "spline")
  }

  l <- c()                                                                               # Length of waveform
  for(i in 1:length(oa)){
    l[i] <- length(p[, (i+1)][!is.na(p[, (i+1)])])
  }

  ipa_ratio <- c()                                                                       # Inflexion point area ratio (For canonical waveform uses x[3], if using inflection point then uses x[4])
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

  pn_time <- c()                                                                         # Systolic peak to Notch time (relative to length of wave)
  for(i in 1:length(oa)){
    pn_time[i] <- (oa[[i]]$x[3] - oa[[i]]$x[2])/l[i]
  }

  nt_ratio <- c()                                                                        # Notch-time ratio = time interval from notch to end of waveform / time interval from notch to beginning of waveform
  for(i in 1:length(oa)){
    wave <- p[, (i+1)][!is.na(p[, (i+1)])]
    next_o <- max(which(is.na(wave) ==  0))                                              # Find next_o i.e the last value of the wave:
    nt_ratio[i] <- (next_o - oa[[i]]$x[3]) / (oa[[i]]$x[3] -  oa[[i]]$x[1])
  }

  ai <- c()                                                                              # Reflectance peak to forward peak ratio (augmentation index)
  for(i in 1:length(oa)){
    ai[i] <- oa[[i]]$y[4] / oa[[i]]$y[2]
  }

  aai <- c()                                                                             # Alternative augmentation index
  for(i in 1:length(oa)){
    aai[i] <- (oa[[i]]$y[2] - oa[[i]]$y[4]) / oa[[i]]$y[2]
  }

  ct <- c()                                                                              # Crest time (O to S time)
  for(i in 1:length(oa)){
    ct[i] <- oa[[i]]$x[2] - oa[[i]]$x[1]
  }

  features <- data.frame(s_vals, n_vals, d_vals, np_ratio,
                         ppt, max_amp, auc, auc_s, l, ipa_ratio,
                         pn_time, nt_ratio, ai, aai, ct)
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
  ##################################################################################################################
  # GetParticipants creates a list of all useable files in the ISO dataset

  # Inputs:
  # direc (file directory for the ISO dataset)

  # Output:
  # string_list (list of ISO participant IDs)
  ###################################################################################################################

  string_list <- list.files(path = direc)                                      # Create list of participant IDs
  string_list <- string_list[-119]
  string_list <- substr(string_list, 0, 5)
  rejected_ts <- c("AL826", "AO139", "AP493", "AU602", "AZ985", "AZ883")       # Remove IDs known to have unusable data
  to_remove <- c()
  for(i in 1:length(rejected_ts)){
    tmp <- which(string_list == rejected_ts[i])
    if(length(tmp) > 0){to_remove <- c(to_remove, tmp)}
  }
  string_list <- string_list[-to_remove]
  return(string_list)
}


GetDirec <- function(run, Participants, dir){
  ##################################################################################################################
  # GetDirec finds the file directory for an individual participant

  # Inputs:
  # run (participant number)
  # Participants (list of participants)
  # dir (file directory for the ISO dataset)

  # Output:
  # direc (file directory for an individual participant)
  ###################################################################################################################
  subjectID <- Participants[run]
  participant_number <- run
  direc <- paste(dir, subjectID, sep = "")
  scan_no <- list.files(direc)[grep(list.files(direc), pattern = "scan")]

  if(length(scan_no) < 1){                                                     # If an 'ISO_ONLY' file another step is needed
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
  ##################################################################################################################
  # GetPairs is an ISO analysis specific function. It is used to identify and order the four data files per
  # participant in the ISO dataset, arranding them into two pairs of 'isoprenaline' time series and 'saline'
  # time series.

  # Inputs:
  # direc (the directory for an individual participant in the ISO dataset)
  # run_order (the order of infusions given during the ISO study protocol)
  # participant_number (the individual participant number for a given participant)
  # subjectID (the participant ID as coded in the ISO study data files)

  # Output:
  # pairs (list of two vectors, each vector containing one ISO time series and one Saline time series)
  ###################################################################################################################
  str <- direc
  direc <- list.files(path = direc)                                            # Identify all files within the participant's directory

  direc <- direc[-grep("REST", direc)]                                         # Refine the directory by removing irrelevant files
  direc <- direc[-grep(".fig", direc)]
  direc <- direc[-grep(".ecgpk", direc)]
  direc <- direc[-grep(".resppk", direc)]
  direc <- direc[-grep(".hrrv", direc)]
  direc <- direc[grep("ECG", direc)]

  dose_numbers <- parse_number(direc)                                          # Identify the infusion order number for each file
  direc <- direc[rev(order(dose_numbers, decreasing = T))]
  if(length(direc) < 6){                                                       # Each directory should now contain 6 files for 6 infusions given
    message("Warning: Some time series are missing from this folder")
  }

  dose_order <- run_order[which(run_order$subj_ID == subjectID), 3]            # By cross-referencing with the infusion run order, identify
  dose_order <- as.numeric(unlist(strsplit(dose_order,",")))                   # which four files correspond to either the ISO 2mcg dose or saline
  direc <- direc[rev(order(dose_order, decreasing = T))]


  set.seed(32)                                                                 # Randomize which files are paired together (each combination of ISO
  # and saline pairs has a probability of 0.5)
  if(rnorm(participant_number)[participant_number] > 0){                       # (set.seed(32) ensures reproducibility)
    pair1 <- c(direc[5], direc[1])
    pair2 <- c(direc[6], direc[2])
  }else{
    pair1 <- c(direc[6], direc[1])
    pair2 <- c(direc[5], direc[2])
  }

  pairs <- list(pair1, pair2)

  return(pairs)
}



UnDetrend <- function(ppg,factor=0,offset=1){
  ##################################################################################################################
  # UnDetrend is a function used as part of ISO data preprocessing. It reverses a detrending effect apparent when
  # viewing the data and presumably originating from the hardware source, according to the inputted factor and
  # offset values.

  # Inputs:
  # ppg (ppg time series to be pre-processed  / undetrended)
  # factor (determines the degree of y-axis decay between data points in the time series)
  # offset (determines the overall trend of the time series)

  # Output:
  # result (an undetrended ppg time series)
  ###################################################################################################################
  k <- offset * (1-factor)                                                     # Define k as the product of offset and inverse factor
  n <- nrow(ppg)
  result <- (1:n)*0                                                            # Create a result vector of same length as inputted time series
  result[1] = ppg[1,2]

  for (i in 2:n){                                                              # Populate the vector with ppg values adjusted for factor and offset
    result[i] = ppg[i,2] - ppg[i-1,2] * factor - k + result[i-1]
  }

  return(result)
}


FactorAdjust <- function(data, factorCutoff, ppg, u, beat, a., test, gs=gs, beatTime, nextTime, plot = T){
  ##################################################################################################################
  # FactorAdjust identifies the most appropriate factor value to be used in undetrending an ISO study ppg time
  # series (for inputting into the Undetrend function). It focuses on a single beat segment of the ppg time
  # series and adjusts the factor value to correct morphological distortions common to detrending in the ISO
  # data. These are 1. abnormal tail decay (e.g positive gradient to tail) 2. values significantly below the
  # approximated baseline (such that notch features are negative and therefore meaningless). The function
  # operates within the FindUndetrendingParams function.

  # Inputs:
  # data (an individual beat segment of raw ppg time series data)
  # factorCutoff (the maximum acceptable tail gradient (default 0))
  # ppg (the ppg time series)
  # u (Undetrend function)
  # beat (dataframe of all detected peaks in the ppg time series)
  # a. (index number indicating which beat (relative to the calculated minimum IBI of the time series) to sample)
  # test (index number indicating the beat identified as the minimum IBI of the time series)
  # gs (model2.GetSegment function)
  # beatTime (index value for the start of the beat to be assessed)
  # nextTime (index value for the end of the beat to be assesssed)
  # plot (logical, if TRUE, plots the beat segment for each iteration of the while loop, demonstrating the change
  # in morphology as it is corrected)

  # Output:
  # factor_value (the appropriate factor value for restoring proper morphology to the ppg time series, to be inputted
  # into the Undetrend function to complete ISO study data preprocessing)

  # Notes:
  # A minimum IBI section of the time series is relevant to this function given its purpose in preprocessing ISO
  # study data. The morpological change associated with isoprenaline results in a lowered notch (and IBI), and thus
  # any alterations to morphology to ensure notches above 0 should be targeted at the waves with the lowest notches.

  # Within the ISO data four types of waveform endings are possible:

  # 1. prolonged negative slope (morphologically normal)
  # 2. prolonged positive slope (due to heavy detrending)
  # 3. ~5 positive slope coming out of the notch, followed by ~5 negative slope as a rather short tail i.e an n shape
  # (these tend to be ISO waves (they are short due to high HR))
  # 4. ~5 positive slope at the end of the tail, preceded by a prolonged negative slope i.e a v shape (due to tail noise)

  # # 1. and 2. will be processed the same regardless of which of the last 10 datapoints the gradient is taken from
  # 3. In the initial 'last 5' datapoint assessment, these will be found to be negative and thus will continue being
  # iterated on based on the last 5 datapoints.
  # 4. In the initial 'last 5' datapoint assessement, these will be found to be positive, and so will be iterated on
  # based on NOT the last 5 datapoints.
  ###################################################################################################################

  tail <- c(data[nrow(data), 2], data[nrow(data)-1, 2],                        # Calculate the gradient of the tail of the beat,
            data[nrow(data)-2, 2], data[nrow(data)-3, 2],                      # from the final 5 values
            data[nrow(data)-4, 2])
  tail <- rev(tail)
  xx. <- 1:length(tail)
  y. <- lm(tail~xx.)

  factor_value <- 1                                                            # Define the factor value as 1 initially (equivalent to no change)

  while(y.[[1]][2] > factorCutoff |                                            # Adjust the factor value until the gradient of the tail reaches an
        (which.min(data[, 2]) > quantile(1:nrow(data))[[2]] &                  # appropriate threshold, or until the minimum value is not the notch
         which.min(data[, 2]) < quantile(1:nrow(data))[[4]])){

    if(factor_value < 0.7){                                                    # If 30 iterations have occurred without the above criteria being met,
      break                                                                    # assume the data is artefactual and end the while loop
    }

    factor_value <- factor_value - 0.01                                        # Reduce the factor value marginally with each iteration

    ppg2 <- ppg                                                                # Create a test data segment with the new factor value implemented
    ppg2[, 2] <- u(ppg,factor=factor_value,offset=1)
    seg <- c(which(ppg$`time (s)` ==  beatTime),
             0, which(ppg$`time (s)` == nextTime))                             # [DEBUG NOTE]: The following lines may be useful if inputs not provided:
    data <- gs(ppg2,seg)                                                       # #beatTime <- beat[test + a., 1] # nextTime <- beat[test + (a. + 1), 1]
    if(plot == TRUE){plot(data, pch = 19)}

    if(y.[[1]][2] > 0){                                                        # If the tail gradient is positive, use values from before the last 5.
      tail <- c(data[nrow(data)-5, 2], data[nrow(data)-6, 2],                  # (If the last 5 were positive due to noise, this will make sure it is
                data[nrow(data)-7, 2], data[nrow(data)-8, 2],                  # counted as negative from that point on)
                data[nrow(data)-9, 2], data[nrow(data)-10, 2],                 # (If the last 5 were genuinely positive, then the positivity should
                data[nrow(data)-11, 2], data[nrow(data)-12, 2],                # extend beyond the last 5)
                data[nrow(data)-13, 2], data[nrow(data)-14, 2])
    }else{
      tail <- c(data[nrow(data), 2], data[nrow(data)-1, 2],                    # If the gradient is negative, use the last 5 values.
                data[nrow(data)-2, 2], data[nrow(data)-3, 2],                  # (Misleadingly negative tails do not tend to occur)
                data[nrow(data)-4, 2], data[nrow(data)-5, 2],
                data[nrow(data)-6, 2], data[nrow(data)-7, 2],
                data[nrow(data)-8, 2], data[nrow(data)-9, 2])
    }
    tail <- rev(tail)
    xx. <- 1:length(tail)
    y. <- lm(tail~xx.)
  }

  return(factor_value)
}


OffsetAdjust <- function(ppg3, ppg, u = UnDetrend, factor_value, plot = F){
  ##################################################################################################################
  # OffsetAdjust identifies the most appropriate offset value to be used in undetrending an ISO study ppg time
  # series (for inputting into the Undetrend function). It is used within the FindUndetrendingParams function
  # once a factor value has been determined.

  # Inputs:
  # ppg3 (the ppg time series adjusted for a given factor value)
  # ppg (the original ppg time series)
  # u (Undetrend function)
  # factor_value (the appropriate factor value for restoring proper morphology to the ppg time series)
  # plot (logical, if TRUE, plots the entire time series for each iteration of offset adjustment, thereby
  # visualising the adjustment as it occurs)

  # Output:
  # offset_value (determines the overall trend of the time series)
  ###################################################################################################################

  if(factor_value == 1){                                                       # If factor value is 1, adjustment is unnecessary
    print("no adjustment needed")
    return(1)
  }

  vv. <- ppg3[, 1]                                                             # Identify the gradient of the time series
  yv. <- lm(ppg3[, 2]~vv.)

  offset_value <- 1                                                            # Begin adjustment by defining offset as 1

  while (yv.[[1]][2] > 0){                                                     # While the time series gradient is > 0, continue to adjust the offset value

    if(yv.[[1]][2] > 5){                                                       # Make more minor adjustments to the offset value as it approaches 0,
      offset_value <- offset_value + 1                                         # to avoid overshooting
    }else{
      if(yv.[[1]][2] > 1){
        offset_value <- offset_value + 0.5
      }else{
        offset_value <- offset_value + 0.05
      }
    }

    ppg3 <- data.frame(ppg[,1], u(ppg,factor=factor_value,                     # Once an adjustment is made, redefine the test time series and recalculate the gradient
                                  offset=offset_value))
    vv. <- ppg3[, 1]
    yv. <- lm(ppg3[, 2]~vv.)
    if(plot == TRUE){
      if(yv.[[1]][2]>0){plot(ppg[,1],
                             u(ppg,factor=factor_value,
                               offset=offset_value), type = "l")}
      if(yv.[[1]][2]>0){abline(a = yv.[[1]][1],
                               b = yv.[[1]][2], col = "red")}
      if(yv.[[1]][2]>0){print(yv.[[1]][2])}
    }
  }

  return(offset_value)
}


FindUndetrendingParams <- function(direc, gs = model2.GetSegment, u = UnDetrend, oa = OffsetAdjust, fa = FactorAdjust, factorCutoff = 0, sr = samplingRate, pk_thrshd, pairs, plot = T){
  ##################################################################################################################
  # FindUndetrendingParams is specific to the ISO study analysis and preprocesses the data therein. However, it
  # may also be useful for preprocessing any PPG time series data in need of morphological correction due to
  # hardware detrending algorithms. It functions by optimising factor and offset undetrending parameters across
  # a number of sample waveforms from the four PPG time series per participant in the ISO study. From this
  # reasonable values to apply to all four time series are inferred.

  # Inputs:
  # direc (the directory for an individual participant in the ISO study)
  # gs (model2.GetSegment function)
  # u (UnDetrend function)
  # oa (OffsetAdjust function)
  # fa (FactorAdjust function)
  # factorCutoff (the maximum acceptable tail gradient (default 0))
  # sr (sample rate)
  # pk_thrshd (the objective threshold for initial identification of peaks in the 1st derivative (see ISO main script))
  # pairs (file names of the two pairs of time series used for each participant in the ISO study analysis)
  # plot (logical, if TRUE plots data during both factor and offset adjustment)

  # Output:
  # values (a vector containing the offset and factor values required to undetrend the raw PPG data)
  ###################################################################################################################
  if(plot == TRUE){p <- TRUE}else{p <- FALSE}

  factor_value_vec <- c()                                                      # First optimise factor value across time series
  for(ps in 1:2){                                                              # (Four waveforms are sampled from each time series)

    if(ps == 1){pair <- pairs[[1]]}else{pair <- pairs[[2]]}

    for(pr in 1:2){
      new_direc <- paste(direc, "/", pair[pr], sep = "")
      ppg <- read.csv(new_direc, sep = "")
      ppg <- data.frame(                                                       # Identify a single time series of the four as 'ppg'
        time = (0:(nrow(ppg)-1)) / samplingRate,
        ppg = ppg[,1]
      )
      names(ppg)[1] <- "time (s)"
      names(ppg)[2] <- "Detrended"

      n <- dim(ppg)[1]                                                         # Identify beats in the time series using peak detection
      vpg <- ppg[2:n,2] - ppg[1:(n-1),2]
      beat <- data.frame(ppg[which(vpg[1:(n-1)] < pk_thrshd &
                                     vpg[2:n] >= pk_thrshd),1])                # If the number of beats found suggests a HR of < 30bpm, trigger warning
      if(nrow(beat) < length(vpg)/(sr*2)){
        message("Warning: Minimal peaks found - consider resetting vpg peak threshold")
        Sys.sleep(10)
      }

      t_value <- c()
      for(i in 1:4){
        print(i)

        pre_ibi <- abs(beat[1:(nrow(beat)-1), 1] - beat[2:nrow(beat), 1])      # Create a vector of inter-beat intervals (IBI) inferred from the peaks detected
        meds <- rollmedian(pre_ibi, k = 19)
        min <- which.min(meds)                                                 # The minimum IBI interval likely represents the timepoint of maximum isoprenaline effect
        if((min - 10) < 1){min = 11}
        if((min + 10) > length(meds)){min = length(meds) - 11}                 # Failsafes in case the minimum is close to the end / beginning of the time series are implemented
        test <- round(quantile((min-10):(min+30)))[[i]]                        # Define the approximation of the minimum of the IBI time series as 'test'

        beatTime <- beat[test, 1]                                              # Extract beat corresponding to said timepoint
        nextTime <- beat[(test + 1), 1]
        seg <- c(which(ppg$`time (s)` ==  beatTime), 0,
                 which(ppg$`time (s)` == nextTime))
        data <- gs(ppg,seg)

        a. <- 0                                                                # If a segment with multiple beats is extracted, move on to the next segment instead
        while (nrow(data) > sr*1.5 | nrow(data) < (0.375*sr) |
               sum(c(order(data[, 2], decreasing = T)[1:5][-1],
                     order(data[, 2], decreasing = T)[1:5][1]) > 20) > 0 ){    # Multiple beats are detected by calculating the time difference between maximum y-value datapoints in the segment
          beatTime <- beat[test + a., 1]
          nextTime <- beat[test + (a. + 1), 1]
          seg <- c(which(ppg$`time (s)` ==  beatTime), 0,
                   which(ppg$`time (s)` == nextTime))
          data <- gs(ppg,seg)
          a. <- a. + 1
        }
        if(a. != 0){
          a. <- a. - 1
        }

        if(plot == TRUE){plot(data, pch = 19)}
        t_value[i] <- fa(data, factorCutoff, ppg, u, beat, a.,                 # Once a valid beat is extracted, adjust the factor value for the beat and store
                         test, gs, beatTime, nextTime, plot = p)               # in vector defined 't_value'
      }

      factor_value_vec <- c(factor_value_vec, t_value)                         # Combine the vectors of optimised factor values across time series into a single vector
    }
  }

  min2 <- order(factor_value_vec)[2]                                           # Arrange the vector in ascending order and select the 2nd from minimum value
  if(sum(min2 == 5:8) > 0 | sum(min2 == 13:16) > 0){
    message("Warning: factor value taken from 0mg time series")                # Check which dose level time series the selected factor value arose from
    Sys.sleep(8)                                                               # (It is expected that the 2nd to minimum value should arise from an isoprenaline beat, which have lower notches)
    min1 <- order(factor_value_vec)[1]
    if(sum(min1 == 5:8) > 0 | sum(min1 == 13:16) > 0){
      message("though minimum factor value found in 2mg time series")
      Sys.sleep(8)
    }
  }
  factor_value <- factor_value_vec[min2]                                       # Final factor value defined


  offset_value <- c()                                                          # Offset Value Adjustment:
  for(ps in 1:2){
    if(ps == 1){pair <- pairs[[1]]}else{pair <- pairs[[2]]}
    for(pr in 1:2){
      new_direc <- paste(direc, "/", pair[pr], sep = "")
      ppg <- read.csv(new_direc, sep = "")
      ppg <- data.frame(                                                       # The four time series are used to optimise the offset value in turn
        time = (0:(nrow(ppg)-1)) / samplingRate,                               # (once adjusted for the now determined factor value)
        ppg = ppg[,1]
      )
      names(ppg)[1] <- "time (s)"
      names(ppg)[2] <- "Detrended"

      ppg3 <- data.frame(ppg[,1],UnDetrend(ppg,
                                           factor=factor_value,
                                           offset=1))
      offset_value <- c(offset_value, oa(ppg3,                                 # A vector of offset values (one per time series) is added to
                                         ppg, u = UnDetrend,
                                         factor_value, plot = p))
    }
  }

  offset_value <- median(offset_value)                                         # The median of offset values is selected as a reasonable approximation for all
  values <- c(factor_value, offset_value)
  return(values)
}



AddOutput <- function(beat){
  ##################################################################################################################
  # AddOutput adds output columns to the beat dataframe (see general purpose script) to be populated with model
  # parameter outputs.

  # Inputs:
  # beat (dataframe where number of rows = number of waveforms in sample)

  # Output:
  # beat (beat dataframe with output columns added)
  ###################################################################################################################
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
  ##################################################################################################################
  # FindStartParams generates starting parameters for the HED model, using estimates derived from the data itself.
  # Since the HED model assumes both decay and excess elements, an approximation of the decay is removed from
  # the original waveform data to infer the excess. Once the excess is defined, peaks in the excess are
  # searched for within empirically defined ranges. Values for width, timing and amplitude for each peak can then
  # be inferred.

  # Inputs:
  # batch_number (the number of batches of waves to process)
  # beats_in (the number of beats in a batch, over which the width, reflectance timing, and decay rate parameters will be fixed)
  # beat (empty dataframe to be filled with starting parameters)
  # ppg (ppg time series (unsegmented))
  # gs (model2.GetSegment function)
  # e (model2.Excess function)
  # sep (model2.SubtractExcessPeak function)
  # o_points (vector of locations of the o (origin) points for each waveform in the ppg time series)
  # wuv (dataframe of locations of the w, u and v points for each waveform in the ppg time series)
  # inflexX (vector of x coordinates for each inflection point on the ppg time series)
  # all_beats (logical, if true ensures all available beats in the time series are modeled)
  # plot (logical, if true plots each beat for which parameters are being estimated)

  # Output:
  # temp (list consisting of the following:)
  # beat (updated dataframe populated with estimated starting parameters, to be used as input for downhill simplex)
  # ppg (additional columns providing parallel time series of inferred baseline, excess and residuals of ppg data)
  ###################################################################################################################
  nBeats <- nrow(beat)
  seg <- c(0,0,0)
  if((batch_number*beats_in) > nBeats){                                                  # Establish number of beats to estimate parameters for
    print("Batch and beat values request more beats than
          are in time series, defaulting to max number of beats")
    maxn <- nBeats
  }else{
    maxn <-(batch_number*beats_in)
    if(all_beats == TRUE){maxn <- maxn + (nrow(beat) - maxn)}                            # [DEBUG NOTE] If concerns that O2 and O points are not identical, can check with:
  }                                                                                      # plot(5900:6500, ppg$Detrended[5900:6500], type = "l")
  # points(o_points, rep(0, length(o_points)))
  # points(inflexX[wuv$o2], rep(0, length(wuv$o2)), col = "red")


  for (i in 1:maxn){                                                                     # For each beat, estimate model parameters:

    beatTime <- beat[i,1]                                                                # Find the position of the ith segmented beat in the ppg time series
    current_o <- which(ppg[, 1] == beatTime)
    next_o <- min(which(o_points > current_o))                                           # Find the minimum o_point that is after the current o_point (to define the limits of the beat within the time series)
    next_o <- o_points[next_o]
    if ((next_o - current_o) < 5){                                                       # Ensure that the distance between o points is not unreasonably small (indicating rounding error)
      next_o <- min(which(o_points > current_o)) + 1
      next_o <- o_points[next_o]
    }
    nextTime <- ppg[round(next_o), 1]                                                    # Find the closest ppg time series x-value to the next o_point
    seg <- c(which(ppg$`time (s)` ==  beatTime), 0,
             which(abs(ppg[, 1] - nextTime) == min(abs(ppg[, 1] - nextTime))))
    data <- gs(ppg,seg)
    if(plot == TRUE){plot(data)}

    tStart <- ppg[seg[1],1]                                                              # Define the first value of the beat
    yPrev <- ppg[max(seg[1]-1,1),2]                                                      # Define the value in the ppg time series immediately preceding the first value of the beat

    amp <- max(data[, 2]) - min(data[, 2])
    constant <- 0.1092254*amp
    baseline <- min(data[,2]) - constant                                                 # Estimate baseline parameter and create a residue time series
    residue <- e(data[,2], ppg[seg[1]-1,2], -0)

    count <- nrow(data)
    excess <- 1:count * 0.0
    excess[1] = data[1,2] - (baseline + config.rate*(yPrev-baseline))                    # Substract estimated decay element from the data to estimate an excess
    for (j in 2:count){
      excess[j] = data[j,2] - (baseline + config.rate*(data[j-1,2]-baseline))
    }
    rm(count)
    rm(j)                                                                                # plot(data[,1],excess, type = "l")
    par <- 1:10 * 0.
    par[1] = baseline
    residue <- excess

    peak.w <- which(data[,1] > beat[i,1]-0.2 & data[,1] < beat[i,1]+0.2)                 # Estimate systolic peak parameters from the excess
    peak.t <- data[peak.w,1]
    peak.y <- residue[peak.w]                                                            # plot(data[,1],excess) # lines(peak.t,peak.y)
    par[3] <- max(peak.y)                                                                # Amplitude of S peak
    par[2] <- peak.t[which(peak.y==par[3])]                                              # Timing of S peak
    par[4] <- 0.25                                                                       # Width of S peak (empirically determined)
    rm(peak.w,peak.t,peak.y)
    residue <- sep(data[,1],residue,par[2:4])                                            # Determine residue (excess subtracted by S peak)
    # plot(data[,1],excess)
    # lines(data[,1],residue)


    peak.w <- which(data[,1] > beat[i,1]+0.3 & data[,1] < beat[i,1]+0.6)                 # Estimate diastolic (2nd reflectance) peak parameters from the excess (minus the S peak)
    peak.t <- data[peak.w,1]                                                             # Find a range within which to locate the D peak
    peak.y <- residue[peak.w]                                                            # plot(data[,1],excess)  # lines(peak.t,peak.y)
    par[6] <- max(peak.y)
    par[5] <- peak.t[which(peak.y==par[6])]
    par[7] <- 0.25
    rm(peak.w,peak.t,peak.y)
    residue <- sep(data[,1],residue,par[5:7])                                            # plot(data[,1],excess)  # lines(data[,1],residue)


    t <- par[2] + c(0.25,0.75) * (par[5]-par[2])                                         # Estimate notch (first reflectance) peak parameters from the excess (minus S and D peaks)
    peak.w <- which(data[,1] > t[1] & data[,1] < t[2])
    peak.t <- data[peak.w,1]
    peak.y <- residue[peak.w]                                                            # plot(data[,1],excess)  # lines(peak.t,peak.y)
    par[9] <- max(peak.y)
    par[8] <- peak.t[which(peak.y==par[9])]
    par[10] <- 0.25
    rm(peak.w,peak.t,peak.y,t)
    residue <- sep(data[,1],residue,par[8:10])                                           # plot(data[,1],excess)  # lines(data[,1],residue)


    w <- seg[1]:seg[3]                                                                   # Store estimated parameters, as well as baseline, excess and residue time series
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
  ##################################################################################################################
  # FindWithinParams generates a 13x12 matrix for each waveform in an inputted batch of waveforms. Each row is a
  # unique parameter set for the waveform in question, while columns 1:12 correspond to each of the 12 parameters
  # (baseline1, baseline2, systolic time, systolic amplitude, systolic width, 2nd reflectance wave time, 2nd
  # reflectance wave amplitude, 2nd reflectance wave width, 1st reflectance wave time, 1st reflectance wave amplitude,
  # 1st reflectance wave width, and decay rate).

  # Each populated row differs from the top row by a single parameter. The specific value that differs in each row
  # represents the optimised value for that parameter when all other parameters are held constant. Rows corresponding
  # to changes in parameters which are to be fixed across beats are unpopulated (NA). The matrices outputted by
  # FindWithinParams are the building blocks for the final matrix that will be inputted into the downhill simplex routine.

  # Inputs:
  # beats_in (number of beats in a batch, over which the width, reflectance timing, and decay rate parameters will be fixed)
  # ppg (the preprocessed ppg time series)
  # beat (the dataframe of model parameter outputs populated with initial starting parameters)
  # gs (model2.GetSegment function)
  # fp (model2.FixParams3 function)
  # ms (simplex.MakeSimplex3 function)
  # m2 (model2.ChiSq3 function)
  # beat_vector (a list consisting of 1. the number of beats in the batch 2. the starting points of each beat on the ppg time series 3. the end points of each beat on the ppg time series)
  # renal_param (the starting parameter for 1st reflectance peak timing (inputted to prevent drastic deviations from this value))
  # dias_param (the starting parameter for 2nd reflectance peak timing (inputted to prevent drastic deviations from this value))
  # sys_time (the starting parameter for systolic peak timing (inputted to prevent drastic deviations from this value))
  # w (the timing of the 1st derivative peaks on the ppg time series)

  # Output:
  # a (list of n 13x12 matrices, where n = number of beats in the inputted batch)
  ###################################################################################################################
  a <- list()
  for(i in 1:beats_in){

    par <- as.numeric(beat[i,5:16])                                                      # Define the 12 parameters as initially estimated (see FindStartParams)
    beat_indi <- list(1, beat_vector[[2]][i], beat_vector[[3]][i])                       # [DEBUG NOTE] par <- fp(data[, 1:2], par, rp = renal_param, sys_t = sys_t) may be required as an intermediate step before this line to fix parameters
    a[[i]] <- ms(ppg = ppg, param = par, f = m2, inScale = 0.1, inTol=-1,                # Generate a matrix for each beat using the simplex.MakeSimplex3 function.
                 beat_vector = beat_indi, renal_param = renal_param,
                 dias_param = dias_param, sys_time = sys_time[i], w = w[i])
  }
  return(a)
}

make_matrix <- function(sim, a){
  ##################################################################################################################
  # make_matrix combines individual waveform matrices for non-fixed parameters and the batch-wide matrix for fixed
  # parameters into a singular matrix for inputting into the downhill simplex routine. The dimensions of the matrix
  # are dependent on the number of beats in the inputted batch, such that matrix dimensions are:

  # (6*n + 6) + 1   x    (6*n + 6)

  # Where n = number of beats in the inputted batch. Therefore a batch of 10 waveforms would generate a matrix of
  # 67x66. The number of columns also represents the number of dimensions the simplex will navigate in optimising
  # parameters across the entire batch.

  # Inputs:
  # sim (the matrix representing variations in fixed parameters across the batch)
  # a (the set of matrices representing variations in non-fixed parameters for each waveform in the batch)

  # Output:
  # sim (a combined matrix representative of both fixed and non-fixed parameters across the batch)
  ###################################################################################################################

  top_row_sim <- sim[1, c(5:6, 8:9, 11:12)]                                    # The top row represents the initial starting parameters that are fixed
  sim <- sim[c(6:7, 9:10, 12:13), c(5:6, 8:9, 11:12)]                          # Redundant (NA) rows and columns are removed
  top_row_sim <- matrix(data = top_row_sim, nrow = 6, ncol = 6, byrow = TRUE)  # The top_row is replcated in matrix form

  top_row <- list()                                                            # Remove redundant rows for each waveform matrix
  for(i in 1:beats_in){                                                        # Save the top row of each waveform matrix
    top_row[[i]] <- a[[i]][1, -c(5:6, 8:9, 11:12)]
    a[[i]] <- a[[i]][-c(1, 6:7, 9:10, 12:13), -c(5:6, 8:9, 11:12)]
  }
  for(i in 1:beats_in){                                                        # Replicate the top row in matrix form
    top_row[[i]] <- matrix(data = top_row[[i]],
                           ncol = 6, nrow = 6, byrow = TRUE)
  }

  for(i in 1:beats_in){                                                        # Assemble Matrix
    sim <- cbind(sim, top_row[[i]])                                            # First combine all the respective 'top rows'
  }

  beat_rows <- list()                                                          # Create additional rows for each beat in the batch
  for(i in 1:beats_in){

    beat_rows[[i]] <- top_row_sim                                              # Add values to the left (i.e from previous beats to the ith)
    if(i != 1){
      for(j in 1:(i-1)){
        beat_rows[[i]] <- cbind(beat_rows[[i]], top_row[[j]])
      }
    }

    beat_rows[[i]] <- cbind( beat_rows[[i]], a[[i]])                           # Add rows corresponding to the ith beat

    if(i == beats_in){                                                         # Add values to the right (i.e from beats succeeding the ith)
      break
    }else{
      for(j in (i+1):beats_in){
        beat_rows[[i]] <- cbind(beat_rows[[i]], top_row[[j]])
      }
    }
  }

  for(i in 1:beats_in){                                                        # Bind rows together
    sim <- rbind(sim, beat_rows[[i]])
  }

  final_top_row <- top_row_sim[1, ]                                            # Add the combined top rows as the final top row
  for(i in 1:beats_in){
    final_top_row <- c(final_top_row, top_row[[i]][1, ])
  }
  sim <- rbind(final_top_row, sim)

  return(sim)
}


extractOutput <- function(beats_in, sim){
  ##################################################################################################################
  # extractOutput extracts the simplex-optimised values for fixed and non-fixed parameters across a given batch.

  # Inputs:
  # beats_in (number of beats in the batch)
  # sim (the final matrix outputted by the downhill simplex routine)

  # Output:
  # temp (list consisting of:)
  # across (optimised fixed parameters across the entire batch (6 in total))
  # within (list of optimised non-fixed parameters for each waveform in the batch)
  ###################################################################################################################
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
  ##################################################################################################################
  # FixOutput applies constraints to the optimised model parameter outputs to ensure they remain within feasible
  # limits, by applying the model2.FixParams3 function.

  # Inputs:
  # beats_in (number of beats in an inputted batch)
  # beat (dataframe of model output parameters (not yet updated with optimised parameters))
  # ppg (ppg time series (unsegmented))
  # gs (model2.GetSegment function)
  # fp (model2.FixParams3 function)
  # across (optimised fixed parameters across the batch)
  # within (optimised non-fixed parameters for each waveform in the batch)
  # sys_time (the starting parameter for systolic peak timing (inputted to prevent drastic deviations from this value))

  # Output:
  # fixed (list of vectors, each vector corresponds to the 12 optimised and fixed parameters for one waveform in the
  # batch)
  ###################################################################################################################
  fixed <- list()
  for(i in 1:beats_in){                                                        # Fix parameters for each beat in the batch
    seg <- c(beat[i,3],0,beat[i,4])
    data <- model2.GetSegment(ppg,seg)
    rm(seg)
    fixed[[i]] <- model2.FixParams3(data,
                                    params = as.numeric(within[[i]]),
                                    across_beat_params = across, sys_t = sys_time[i])
  }
  return(fixed)
}


UpdateBeat <- function(beats_in, beat, fixed){
  ##################################################################################################################
  # UpdateBeat generates a new matrix of model parameters by replacing values in the existing dataframe with
  # parameters that are both optimised and constrained.

  # Inputs:
  # beats_in (number of beats in the batch)
  # beat (existing dataframe of model parameter estimates (populated with initial estimations))
  # fixed (list of optimised and constrained parameters for each waveform, to replace initial estimates)

  # Output:
  # new_beat (updated dataframe of model parameter outputs)
  ###################################################################################################################
  new_beat <- data.frame(matrix(0, ncol = 12, nrow = beats_in))
  for(i in 1:beats_in){
    new_beat[i, ] <- fixed[[i]]
  }
  new_beat <- cbind(beat[, 1:4], new_beat)
  return(new_beat)
}


FixBaseline <- function(new_beat, f = model2.ChiSq3, renal_param, dias_param, sys_time, w){
  ################################################################################################################################
  # FixBaseline assesses the first and second fitted baselines for each waveform in a batch. If the baselines are close in value
  # and the reduced ChiSq value (goodness of fit) for a waveform's parameter set is lower when the baselines are equal (i.e
  # baseline2 equal to baseline1), then baseline2 will be made equal to baseline1.

  # Inputs:
  # new_beat (dataframe of model parameter outputs)
  # f (model2.ChiSq3 function)
  # renal_param (the starting parameter for 1st reflectance peak timing (inputted to prevent drastic deviations from this value))
  # dias_param (the starting parameter for 2nd reflectance peak timing (inputted to prevent drastic deviations from this value))
  # sys_time (the starting parameter for systolic peak timing (inputted to prevent drastic deviations from this value))
  # w (the timing of the 1st derivative peaks on the ppg time series)

  # Output:
  # new_beat (dataframe of model parameter outputs with baseline parameters further corrected)
  #################################################################################################################################
  for(j in 1:nrow(new_beat)){                                                            # Baselines are assessed for each beat in the batch
    if(abs(new_beat[j, 6] - new_beat[j, 5]) < 5){                                        # Assess if baselines 1 and 2 are close in value
      wave_check <- model2.ChiSq3(data = ppg,                                            # If so, assess goodness of fit with unique values
                                  params = as.numeric(new_beat[j, 5:16]),
                                  beats = list(1, new_beat[j, 3], new_beat[j, 4]),
                                  beat = NULL, a = NULL, plot = FALSE,
                                  renal_param = renal_param, dias_param = dias_param,
                                  sys_time = sys_time[j], w = w[j])
      wave_check2 <- model2.ChiSq3(data = ppg,                                           # Then assess goodness of fit with baselines equal
                                   params = c(rep(new_beat[j, 5], 2),
                                              as.numeric(new_beat[j, 7:16])),
                                   beats = list(1, new_beat[j, 3],
                                                new_beat[j, 4]), beat = NULL, a = NULL,
                                   plot = FALSE, renal_param = renal_param,
                                   dias_param = dias_param,
                                   sys_time = sys_time[j], w = w[j])
      if(wave_check2 < wave_check){                                                      # If equal values give a better fit, fix them to be so
        new_beat[j, 6] <- new_beat[j, 5]
      }
    }
  }
  return(new_beat)
}


PlotFits <- function(beats_in, ppg, beat2, gs = model2.GetSegment, rb = model2.Rebuild2){
  ##################################################################################################################
  # PlotFits produces base R plots of each modelled waveform in the batch, including component waves and baselines.

  # Inputs:
  # beats_in (number of beats in the batch)
  # ppg (ppg time series)
  # beat2 (finalised dataframe of model parameter outputs for each waveform in the batch)
  # gs (model2.GetSegment function)
  # rb (model2.Rebuild2 function)

  # Output:
  # Automated plotting of each fitted waveform in the batch
  ###################################################################################################################
  for(i in 1:beats_in){
    seg <- c(beat[i,3],0,beat[i,4])                                            # Identify beat within ppg time series
    data <- model2.GetSegment(ppg,seg)
    yPrev <- ppg[seg[1]-1,2]
    xPrev <- ppg[seg[1]-1, 1]
    xNext <- ppg[seg[3], 1]
    rm(seg)
    temp <- model2.Rebuild2(data, yPrev, as.double(beat2[i,]),TRUE)            # Generate the fitted wave from parameter outputs
    plot(data[, 1], data[, 2], ylim = c(beat2$Baseline[1]*1.5,                 # Plot both original data wave and fitted wave
                                        max(data[, 2]*1.2)), main = paste(c("batch",
                                                                            k, "wave", i), collapse = " "))
    lines(data[,1],temp)
    lines(c(xPrev, (beat2[i, 3]  + (1*beat2[i, 6]))), rep(beat2[i, 1], 2))     # Plot baselines
    lines(c((beat2[i, 3]  + (1*beat2[i, 6])), xNext), rep(beat2[i, 2], 2))

    par <- as.double(beat2[i,])                                                # Plot systolic wave
    par[c(7:8, 10:11)] <- 0
    temp<-model2.Rebuild2(data,yPrev,par,TRUE)
    lines(data[,1],temp, col = "red")

    par <- as.double(beat2[i,])                                                # Plot second reflectance wave
    par[c(4:5, 10:11)] <- 0
    temp<-model2.Rebuild2(data,yPrev,par,TRUE)
    lines(data[,1],temp, col = "blue")

    par <- as.double(beat2[i,])                                                # Plot first reflectance wave
    par[c(4:5, 7:8)] <- 0
    temp <- model2.Rebuild2(data,yPrev,par,TRUE)                               # [DEBUG NOTE] When plotting component waves, 12 elements are always needed;
    lines(data[,1],temp, col = "green")                                        # the amplitude and width for the peaks that aren't being used should be set to 0
  }
}


GGplotFits <- function(beats_in, ppg, beat2, gs = model2.GetSegment, rb = model2.Rebuild2, run, pr, p = F){
  ##################################################################################################################
  # GGplotFits produces GGplot plots of each modelled waveform in the batch, including component waves and baselines.

  # Inputs:
  # beats_in (number of beats in the batch)
  # ppg (ppg time series)
  # beat2 (finalised dataframe of model parameter outputs for each waveform in the batch)
  # gs (model2.GetSegment function)
  # rb (model2.Rebuild2 function)
  # run (participant number)
  # pr (dose level)
  # p (logical, setting to true stores the plot as a variable)

  # Output:
  # Automated plotting of each fitted waveform in the batch
  ###################################################################################################################
  if(pr == 1){pr = 0}

  for(i in 1:beats_in){
    seg <- c(beat[i,3],0,beat[i,4])                                            # Generate each fitted wave and the original
    data <- model2.GetSegment(ppg,seg)
    yPrev <- ppg[seg[1]-1,2]
    xPrev <- ppg[seg[1]-1, 1]
    xNext <- ppg[seg[3], 1]
    rm(seg)
    temp<-model2.Rebuild2(data, yPrev, as.double(beat2[i,]),TRUE)
    fit <- model2.Rebuild2(data, yPrev, as.double(beat2[i,]),TRUE)

    waves <- data.frame(data)                                                  # Create a dataframe of each plot element
    waves <- cbind(data, fit)
    par <- as.double(beat2[i,])
    par[c(7:8, 10:11)] <- 0
    temp <- model2.Rebuild2(data,yPrev,par,TRUE)
    waves <- cbind(waves, temp)
    par <- as.double(beat2[i,])
    par[c(4:5, 7:8)] <- 0
    temp <- model2.Rebuild2(data,yPrev,par,TRUE)
    waves <- cbind(waves, temp)
    par <- as.double(beat2[i,])
    par[c(4:5, 10:11)] <- 0
    temp <- model2.Rebuild2(data,yPrev,par,TRUE)
    waves <- cbind(waves, temp)

    b1y <- beat2[i, 1]                                                         # Create baselines for adding with geom_segment()
    b1x <-  c(xPrev, (beat2[i, 3]  + (1*beat2[i, 6])))
    b2x <- c((beat2[i, 3]  + (1*beat2[i, 6])), xNext)
    b2y <-  beat2[i, 2]

    sfunction <- splinefun(1:nrow(waves), fit, method = "natural")             # Convert component waves and fitted wave into splines
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

    time <- waves[, 1]                                                         # Compile elements into a stacked dataframe
    time <- seq(from = time[1], to = time[length(time)],
                length.out = length(fit2))
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

    data_stacked <- data.frame(data)                                           # Stack data
    data_stacked <- cbind(data_stacked, rep("data"))
    colnames(data_stacked) <- c("x", "values", "Wave")

    if(p == T){                                                                # Plot
      c <- ggplot(data = waves_stacked_final,
                  aes(x = x, y = values, col = Wave)) +
        geom_line(aes(size = Wave, alpha = Wave)) +
        scale_color_manual(values = c("#03fc7b",
                                      "#03b5fc", "black", "black", "black", "#ff4242",
                                      "black")) +
        scale_size_manual(values = c(0.7, 0.7, 1.5, 0.7, 0.7)) +
        scale_alpha_manual(values = c(1, 1, 1, 1, 1)) +
        ylab("PPG Signal") + xlab("Time") +
        geom_point(data = data_stacked) +
        geom_segment(aes(x = data[1, 1], y = b1y, xend = b1x[2],
                         yend = b1y, colour = "black")) +
        geom_segment(aes(x = b2x[1], y = b2y, xend = b2x[2],
                         yend = b2y, colour = "black")) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black")) +
        theme(legend.position = "none") +
        ggtitle(paste0("Participant", " ", run, "\n", pr,
                       " ", "mcg", sep = ""))
    }else{
      ggplot(data = waves_stacked_final, aes(x = x, y = values,
                                             col = Wave)) +
        geom_line(aes(size = Wave, alpha = Wave)) +
        scale_color_manual(values = c("green", "blue", "black",
                                      "black", "black", "red", "black")) +
        scale_size_manual(values = c(0.7, 0.7, 1.5, 0.7, 0.7)) +
        scale_alpha_manual(values = c(1, 1, 1, 1, 1)) +
        ylab("PPG Signal") + xlab("Time") +
        geom_point(data = data_stacked) +
        geom_segment(aes(x = b1x[1], y = b1y, xend = b1x[2],
                         yend = b1y, colour = "black")) +
        geom_segment(aes(x = b2x[1], y = b2y, xend = b2x[2],
                         yend = b2y, colour = "black")) +
        theme(legend.position = "none")
    }
  }
}



osnd_fit <- function(bf = beat_final, ppg, gs = model2.GetSegment, r = model2.Rebuild2, sf = splinefun, dp = diast_pk, oa = osnd_of_average, sr = samplingRate, plot = FALSE){
  ##################################################################################################################
  # osnd_fit identifies OSND points on both original waveforms and model-generated waveforms, such that the HED
  # model's ability to recapitulate important morphological features and fiducial points can be assessed.

  # Inputs:
  # bf (finalised model parameter outputs)
  # ppg (ppg time series)
  # gs (model2.GetSegment function)
  # r (model2.Rebuild2 function)
  # sf (splinefun function)
  # dp (diast_pk function)
  # oa (osnd_of_average function)
  # sr (sampling rate)
  # plot (logical, plots each waveform with modelled version superimposed, as well as OSND points for each)

  # Output:
  # osnd_diff (list of values representing the error between model-generated OSND values and actual OSND values)
  ###################################################################################################################
  osnd_diff <- list()
  for(i in 1:nrow(bf)){
    seg <- c(bf[i,3],0,bf[i,4])                                                # Identify data segments and corresponding model fit for each waveform.
    data <- gs(ppg,seg)
    yPrev <- ppg[seg[1]-1,2]
    xPrev <- ppg[seg[1]-1, 1]
    xNext <- ppg[seg[3], 1]
    rm(seg)
    temp <- r(data, yPrev, as.double(bf[i,-c(1:4)]),TRUE)

    sfunction <- sf(1:length(temp), temp, method = "natural")                  # Upsample modelled waves to match precision of fiducial point finding on original waveforms
    fit <-  sfunction(seq(1, length(temp), 0.1), deriv = 0)
    sfunction <- sf(1:length(data[, 2]), data[, 2], method = "natural")        # Upsample data segment
    dat <-  sfunction(seq(1, length(data[, 2]), 0.1), deriv = 0)               # plot(dat)  # plot(fit)

    tmp <- dp(avw = fit, sr = sr, scale = T, dias_param = bf[i, 10]*sr*10)     # Find OSND points on the modelled waveform
    dPeak <- tmp[1]                                                            # pass in 2nd reflectance parameters to help find D on the modelled waveform
    xShift <- tmp[2]
    rm(tmp)
    osnd_fit <- oa(fit, dp = dPeak, diff = 0, sr = sr, plot = F)

    tmp <- dp(avw = dat, sr = sr, scale = T, dias_param = bf[i, 10]*sr*10)     # Find OSND points on the original waveform
    dPeak <- tmp[1]
    xShift <- tmp[2]
    rm(tmp)
    osnd_dat <- oa(dat, dp = dPeak, diff = 0, sr = sr, plot = F)

    osnd_dat$x <- osnd_dat$x/(10*sr)                                           # Adjust x-axis for upsampling and sampling rate, then find the difference between OSND points
    osnd_fit$x <- osnd_fit$x/(10*sr)                                           # on the fitted / modelled waveform vs those on the original waveform
    osnd_diff[[i]] <- osnd_dat - osnd_fit

    if(plot == TRUE){                                                          # Plot both original and modelled waveforms along with OSND points for each
      plot((1:length(dat))/(10*samplingRate), dat, type = "l",
           xlab = "time", ylab = "")
      lines((1:length(dat))/(10*samplingRate), fit, col = "red")
      points(osnd_dat, pch = 19)
      points(osnd_fit, col = "red", pch = 19)
    }
  }
  return(osnd_diff)
}


ArrangeOutputs <- function(beat_final, beat_orig, features, pulse, fit_check, ps, pr){
  ##################################################################################################################
  # ArrangeOutputs organises lists of outputs from various parts of the processing pipeline into a single output list.

  # Inputs:
  # beat_final (dataframe of finalised model parameter outputs)
  # beat_orig (dataframe of intially estimated model parameter outputs)
  # features (dataframe of morphological features for each waveform)
  # pulse (dataframe of all individual waveforms (discrete form))
  # fit_check (list of goodness of fit measures for each modelled waveform (ChiSq, Max error, NRMSE, aNRMSE))
  # ps (time series number (for ISO data))
  # pr (dose level (for ISO data))

  # Output:
  # tmp (list composed of:)
  # beat_final (as above)
  # features (as above)
  # fit_check (as above)
  ###################################################################################################################

  rownames(beat_final) <- colnames(pulse)[-1][1:nrow(beat_final)]              # Rename rows for consistency across outputs
  rownames(features) <- colnames(pulse)[-1]

  wave_fits <- c()                                                             # Reorganize fit_check output
  for(i in 1:length(fit_check)){
    wave_fits <- c(wave_fits, as.numeric(fit_check[[i]][[2]]))
  }
  max_err <- c()
  for(i in 1:length(fit_check)){
    max_err <- c(max_err, as.numeric(fit_check[[i]][[3]]))
  }
  NRMSE <- c()
  for(i in 1:length(fit_check)){
    NRMSE <- c(NRMSE, as.numeric(fit_check[[i]][[4]]))
  }
  aNRMSE <- c()
  for(i in 1:length(fit_check)){
    aNRMSE <- c(aNRMSE, as.numeric(fit_check[[i]][[5]]))
  }
  fit_check <- list(wave_fits, max_err, NRMSE, aNRMSE)

  tmp <- list(beat_final, features, fit_check)                                 # Combine output lists
  return(tmp)
}


model2.GetSegment <- function(ppg, limits){
  ##################################################################################################################
  # model2.GetSegment extracts an original waveform from the preprocessed ppg time series. It is used by a number
  # of other functions.

  # Inputs:
  # ppg (ppg time series)
  # limits (limits of segment of time series to be extracted)

  # Output:
  # result (segment of ppg time series corresponding to a single waveform)
  ###################################################################################################################
  w <- c(limits[1]:limits[3])
  result <- matrix(nrow=length(w),ncol=2)
  result[,1] <- ppg[[lab.time]][w]
  result[,2] <- ppg[[lab.ppg]][w]

  return(result)
}

model2.Excess <- function(y, offset, baseline){
  ##################################################################################################################
  # model2.Excess converts a segmented waveform into its excess element by subtracting an approximated decay
  # element.

  # Inputs:
  # y (the y values of the waveform, as a vector)
  # offset (the y value of the sample in the time series directly preceding the first value of the inputted waveform)
  # baseline (anticpated baseline towards which the decay element decays)

  # Output:
  # result (the approximated excess element of the waveform in question)
  ###################################################################################################################
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

model2.Peak <- function(time, peakParams){
  ##################################################################################################################
  # model2.Peak generates a sine wave from the three parameters used to define each component wave (timing, width,
  # and amplitude).

  # Inputs:
  # time (x-axis values of the sine wave)
  # peakParams (vector of three parameters defining a given peak)

  # Output:
  # result (sine wave constructed from model peak parameters)
  ###################################################################################################################
  temp <- 2*const.pi*(time - as.double(peakParams[1]))/as.double(peakParams[3])
  temp[which(temp < -const.pi)] = -const.pi
  temp[which(temp >  const.pi)] =  const.pi
  result <- as.double(peakParams[2]) * (0.5 * (1+cos(temp)))^2

  return(result)
}

model2.SubtractExcessPeak <- function(time, residue, peakParams){
  ##################################################################################################################
  # model2.SubtractExcessPeak removes modelled component sine waves from excess segments of the ppg time series.

  # Inputs:
  # time (x-values of ppg time series)
  # residue (excess component of ppg time series (before or after the extraction of other component waves))
  # peakParams (vector of three parameters to construct the component sine wave to be removed)

  # Output:
  # result (excess time series with modelled sine wave subtracted (residue))
  ###################################################################################################################
  result <- residue - model2.Peak(time,peakParams)

  return(result)
}


model2.ChiSq3 <- function(data, params, debug=FALSE, beats, optional = NULL, beat = NULL, a = NULL, plot = FALSE, renal_param, dias_param, sys_time, w){
  #################################################################################################################################
  # model2.ChiSq3 calculates a reduced ChiSq (goodness of fit) value for a model generated waveform, when provided
  # with the original waveform being modelled and a set of model parameters. The function can either be applied to
  # calculate the goodness of fit of an individual waveform, or batches of waveforms, and is used to drive the
  # downhill simplex optimisation of model parameters.

  # Inputs:
  # data (section of ppg time series)
  # params (model parameters)
  # debug (logical, currently redundant)
  # beats (list of number of inputted beats, there beginnings in the time series, and their endings)
  # optional (logical, now redundant)
  # beat (dataframe of model parameters)
  # a (the combined matrix used in the downhill simplex routine (alternative source of parameters))
  # plot (logical, if set to true plots model generated waveform against original waveform)
  # renal_param (the starting parameter for 1st reflectance peak timing (inputted to prevent drastic deviations from this value))
  # dias_param (the starting parameter for 2nd reflectance peak timing (inputted to prevent drastic deviations from this value))
  # sys_time (the starting parameter for systolic peak timing (inputted to prevent drastic deviations from this value))
  # w (the timing of the 1st derivative peaks on the ppg time series)

  # Output:
  # ts_fit (total reduced ChiSq value summed across all beats inputted)

  # Notes:
  # makesimplex3 inputs a single set of parameters (tParam) to test, whereas run.simplex2 inputs a matrix
  # This determines where within beat parameters should be extracted from, hence beat and a are NULL unless otherwise specified.

  # During the downhill simplex parameter optimization process, movements of parameter values outside of defined contraints
  # will lead to the application of penalties, which model2.ChiSq3 adds to the ChiSq value outputted, driving the
  # optimisation process to occur within constraints. The penalties are proportionate to the degree to which the constraints
  # are exceeded, to prevent mere ignorance of the constraint once violated.
  ################################################################################################################################

  if(!is.null(a)){                                                             # Extract across-beat parameters
    across_beat_params <- a[1:6]
  }else{                                                                       # If a simplex matrix has been supplied, extract from that
    par <- params
    across_beat_params <- par[c(5, 6, 8, 9, 11, 12)]                           # If not, extract from the params input
  }

  beat_fit <- list()                                                           # Calculate ChiSq for all beats:
  for(i in 1:beats[[1]]){

    if(!is.null(a)){                                                           # Extract within-beat parameters (this happens within the for loop since they are specific to each individual waveform)
      par2 <- a[((i*6)+1):((i*6)+6)]
    }else{                                                                     # If a simplex matrix has been supplied, extract from that
      if(!is.null(beat)){                                                      # If not, extract from the beat dataframe
        par2 <- as.numeric(beat[i, c(5:8, 11, 14)])
      }else{                                                                   # If the beat dataframe is also not provided, extract from the params input
        par2 <- par[c(1:4, 7, 10)]
      }
    }

    seg <- c(beats[[2]][i],0,beats[[3]][i])                                    # Identify the individual ith waveform segment
    dat <- model2.GetSegment(data,seg)
    rm(seg)

    sys <- par2[3]                                                             # Extract systolic and diastolic parameters for the ith waveform
    dias <- par2[3] + dias_param                                               # [DEBUG NOTE: # dias <- across_beat_params[2] + par2[3] is an alternative if not wishing to input dias_param]
    end <- which(abs(dat[, 1]-dias) == min(abs(dat[, 1] - dias)))

    w. <- w[i]                                                                 # Find w (first derivative peak) on the ith waveform
    w. <- which(abs(dat[, 1] - w.) == min(abs(dat[, 1] - w.)))

    sys_t <- sys_time[i]
    start <- which(abs(dat[, 1]-sys_t) == min(abs(dat[, 1] - sys_t)))

    temp <- model2.FIX_PAR3(time = dat[, 1], within_beat_params = par2,        # Fix parameters and calculate penalty for those exceeding constraints
                            across_beat_params = across_beat_params,
                            renal_param = renal_param, sys_t = sys_t)
    penalty <- temp[1]
    fixedPar <- temp[2:length(temp)]
    rm(temp)

    fit <- model2.Rebuild2(dat,dat[1,2],params = fixedPar)                     # Generate modelled wave and calculate residual error (original wave - model generated wave)
    residue <- dat[ ,2] - fit

    if(residue[start]*residue[start] > 0){                                     # Add to the penalty if the residual around the time of the systolic peak (sys_t) is high
      penalty <- penalty + residue[start]*residue[start]
    }

    residue[w.:end[1]] <-  residue[w.:end[1]]*3                                # Weight the residual error in the morphologically relevant part of the waveform more strongly (weighted region is W -> D (with slope))
    if(length(residue) > end[1]){
      tail <- (end[1]+1):length(residue)
      for(j in 1:length(tail)){
        wgt <- 3 - (0.1*j)
        if(wgt < 1){wgt <- 1}
        residue[tail[j]] <- residue[tail[j]]*wgt
      }
    }

    nData <- nrow(dat)                                                         # Calculate Reduced Chi-Square for the ith beat, and add the calculated penalty for that beat also
    nPar <- length(par2)
    beat_fit[[i]] <- (sum(residue*residue) / (nData-nPar)) + as.numeric(penalty)

    if(plot == TRUE){
      plot(dat,  ylim = c(-150, 2000))                                         # ylim = c(76, 86) for bioradio data, ylim = c(-150, 1600) for ISO
      lines(dat[, 1], fit)                                                     # lines(dat[, 1], residue + dat[1, 2]) can also be added
    }
  }

  temp <- c()                                                                  # Summate all individual beat ChiSq values
  for(i in 1:length(beat_fit)){
    temp[i] <- beat_fit[[i]][1]
  }
  ts_fit <- sum(temp)
  return(ts_fit)
}


model2.ChiSq4 <- function(data, params, debug=FALSE, beats, beat, a = NULL, plot = FALSE, renal_param, dias_param, sys_time, w){
  #################################################################################################################################
  # model2.ChiSq4 is a clone of model2.ChiSq3. It differs by providing additional measures of goodness of fit besides reduced
  # ChiSq for the inputted batch. These are listed below.

  # Inputs:
  # data (section of ppg time series)
  # params (model parameters)
  # debug (logical, currently redundant)
  # beats (list of number of inputted beats, there beginnings in the time series, and their endings)
  # beat (dataframe of model parameters)
  # a (the combined matrix used in the downhill simplex routine (alternative source of parameters))
  # plot (logical, if set to true plots model generated waveform against original waveform)
  # renal_param (the starting parameter for 1st reflectance peak timing (inputted to prevent drastic deviations from this value))
  # dias_param (the starting parameter for 2nd reflectance peak timing (inputted to prevent drastic deviations from this value))
  # sys_time (the starting parameter for systolic peak timing (inputted to prevent drastic deviations from this value))
  # w (the timing of the 1st derivative peaks on the ppg time series)

  # Output:
  # fit (list consisting of:)
  # ts_fit (total reduced ChiSq value summed across all beats inputted)
  # beat_fit (list of reduced ChiSq values for each individual beat)
  # max_error (residual of greatest value for each individual beat)
  # NRMSE (Normalised root mean square error (see supplementary material))
  # aNRMSE (alternative normalised root mean square error (see supplementary material))
  ################################################################################################################################

  if(!is.null(a)){                                                             # Across-beat parameter extraction
    across_beat_params <- a[1:6]
  }else{                                                                       # If a simplex matrix has been supplied, extract from that
    par <- params                                                              # If not, extract from the params input
    across_beat_params <- par[c(5, 6, 8, 9, 11, 12)]
  }

  beat_fit <- c()                                                              # Calculate ChiSq for all beats:
  max_error <- c()
  NRMSE <- c()
  aNRMSE <- c()
  for(i in 1:beats[[1]]){                                                      # The number of beats is determined by the first object of beats

    # Within-beat parameter extraction:
    if(!is.null(a)){                                                           # If a simplex matrix has been supplied, extract from that
      par2 <- a[((i*6)+1):((i*6)+6)]
    }else{
      if(!is.null(beat)){                                                      # If not, extract from beat.
        par2 <- as.numeric(beat[i, c(5:8, 11, 14)])
      }else{                                                                   # If beat is also not provided, take them from the params input
        par2 <- par[c(1:4, 7, 10)]
      }
    }

    seg <- c(beats[[2]][i],0,beats[[3]][i])                                    # Extract individual beat data
    dat <- model2.GetSegment(data,seg)
    rm(seg)

    sys <- par2[3]                                                             # Extract systolic and diastolic parameters
    dias <- par2[3] + dias_param                                               # [DEBUG NOTE] dias <- across_beat_params[2] + par2[3] can be used if dias_param not inputted
    start <- which(abs(dat[, 1]-sys) == min(abs(dat[, 1] - sys)))
    end <- which(abs(dat[, 1]-dias) == min(abs(dat[, 1] - dias)))

    w. <- w[i]                                                                 # Find W
    w. <- which(abs(dat[, 1] - w.) == min(abs(dat[, 1] - w.)))

    sys_t <- sys_time[i]                                                       # Define intially estimated systolic timing and amplitude

    temp <- model2.FIX_PAR3(time = dat[, 1], within_beat_params = par2,        # Fix parameters and calculate penalty
                            across_beat_params = across_beat_params,
                            renal_param = renal_param, sys_t = sys_t)
    penalty <- temp[1]
    fixedPar <- temp[2:length(temp)]
    rm(temp)

    fit <- model2.Rebuild2(dat,dat[1,2],params = fixedPar)                     # Calculate fit, residue and max error
    residue <- dat[ ,2] - fit
    max_error[i] <- max(residue)

    if(plot == TRUE){
      plot(dat, ylim = c(-500, 2500))
      lines(dat[, 1], fit)
      lines(dat[, 1], residue, col = "Red")
    }

    rmse_begin <- floor((1+w.)/2)                                              # Define region of interest (Begin half way from O to W, End 10 points after the d-peak (this corresponds to half way down the weighted tail))
    rmse_end <- end[1] + 10
    if(sum(is.na(residue[rmse_begin:rmse_end])) > 0){
      residue_roi <- residue[rmse_begin:length(residue)]                       # If there are fewer than 10 data points after D, use as many as there are
      ind_resid <- rmse_begin:length(residue)
    }else{
      residue_roi <- residue[rmse_begin:rmse_end]
      ind_resid <- rmse_begin:rmse_end
    }
    if(plot == TRUE){lines(dat[ind_resid, 1],residue_roi, col = "green")}

    fit_null <- rep(mean(dat[ind_resid, 2], trim = 0), length(residue_roi))    # Define null model (mean of the waveform)
    if(plot == TRUE){lines(dat[ind_resid, 1], fit_null)}

    residuals_of_null_model <- dat[ind_resid, 2] - fit_null                    # Calculate residuals of the null model

    rmse_model2 <-  sqrt(mean(residue_roi^2, trim = 0))                        # Calculate NRMSE
    rmse_null <- sqrt(mean(residuals_of_null_model^2))
    NRMSE. <- 1 - (rmse_model2 / rmse_null)
    NRMSE[i] <- NRMSE.

    aNRMSE[i] <- (sum(residue_roi^2) / sum(dat[ind_resid, 2]^2))*100           # Calculate alternative NRMSE method (Wang et al 2013) (SSE / Sum of squared datapoints)
    # Optional plotting: # plot(dat[ind_resid, 1], dat[ind_resid, 2]^2, type = "l")  # lines(dat[ind_resid, 1], residue_roi^2, col = "red")

    residue[w.:end[1]] <-  residue[w.:end[1]]*3                                # Weighted region is W -> D (with slope)
    if(length(residue) > end[1]){
      tail <- (end[1]+1):length(residue)
      for(j in 1:length(tail)){
        wgt <- 3 - (0.1*j)
        if(wgt < 1){wgt <- 1}
        residue[tail[j]] <- residue[tail[j]]*wgt
      }
    }

    nData <- nrow(dat)                                                         # Calculate Reduced Chi-Square for the beat
    nPar <- length(par2)
    if(par2[1] == par2[2]){                                                    # If baselines are the same, consider them as 1 parameter
      nPar <- nPar - 1
    }
    beat_fit[i] <- (sum(residue*residue) / (nData-nPar)) + as.numeric(penalty)
  }

  temp <- c()                                                                  # Summate individual beat ChiSq values
  for(i in 1:length(beat_fit)){
    temp[i] <- beat_fit[i]
  }
  ts_fit <- sum(temp)
  rm(temp)

  fit <- list(ts_fit, beat_fit, max_error, NRMSE, aNRMSE)
  return(fit)
}


model2.Rebuild2 <- function(xy,offset,params,invert=TRUE){
  #################################################################################################################################
  # model2.Rebuild2 takes HED model parameters and reconstructs a segment of the ppg time series corresponding to one waveform.
  # First the excess is calculated from the 9 peak parameters, then the decay element is added to generate a final waveform.

  # Inputs:
  # xy (x and y values for a particular segment of the ppg time series)
  # offset (the y value of the data point immediately preceding the inputted segment of the ppg time series)
  # params (vector of 12 parameters for generating a modelled wave)
  # invert (logical, if set to true adds decay element as well as excess element)

  # Output:
  # result (dataframe of x and y values defining the model generated waveform)
  ################################################################################################################################
  result <- 1:nrow(xy) * 0.0                                                   # Create output dataframe

  result <- result + model2.Peak(xy[,1],params[3:5])                           # Add systolic peak
  if (length(params)>=8){
    result <- result + model2.Peak(xy[,1],params[6:8]+c(params[3],0,0))        # Add diastolic (2nd reflectance) peak
  }
  if (length(params)>=11){
    result <- result + model2.Peak(xy[,1],params[9:11]+c(params[3],0,0))       # Add renal (1st reflectance) peak
  }

  if (invert){                                                                 # Add decay (config.rate + baseline parameters)
    result <- model2.Excess.Inv2(xy[,1],result,offset,params[1],params[2],
                                 params[3]+1*params[6],
                                 config.rate = params[12])
  }
  return(as.double(result))
}


model2.Excess.Inv2 <- function(time,excess,offset,baselineStart,baselineEnd,timeBase,config.rate){
  #################################################################################################################################
  # model2.Excess.Inv2 adds the model's decay element to its excess element. It is used within model2.Rebuild2 as part of the
  # construction of a modelled waveforms from inputted model parameters.

  # Inputs:
  # time (x coordinates of a particular segment of the PPG time series)
  # excess (segment corresponding to excess element of reconstructed waveform)
  # offset (the y value of the data point immediately preceding the inputted segment of the ppg time series)
  # baselineStart (the y-value towards which the first baseline decays)
  # baselineEnd (the y-value towards which the second baseline decays)
  # timeBase (the point in time (x value) when baseline 1 switches to baseline 2, defined as the D (2nd reflectance) peak timing)
  # config.rate (decay rate parameter)

  # Output:
  # result (a modelled waveform with combined excess and decay elements)
  ################################################################################################################################
  nX <- length(excess)
  if (nX == 0){
    print("Help")
  }
  result <- 1:nX * 0.0                                                         # Define length of output as length of the excess
  baseline <- time * 0 + baselineStart                                         # Define baselines 1 and 2
  baseline[which(time > timeBase)] = baselineEnd

  temp <- which(is.nan(excess))                                                # NA values in the excess will not suffice, replace with 0s
  if(length(temp) > 0){
    for(i in 1:length(temp)){
      excess[temp][i] <- 0
    }
  }

  result[1] = excess[1] + (baselineStart +                                     # Adding the decay element to the first value of the excess (defined in terms of decay rate and baseline parameters)
                             config.rate*(offset-baselineStart))
  for (j in 2:nX){                                                             # Add decay elements to the remainder of excess values (each in turn is dependent on the previously generated value)
    result[j] = excess[j] + (baseline[j] +
                               config.rate*(result[j-1]-baseline[j]))
  }
  return(result)
}


model2.FIX_PAR3 <- function(time, within_beat_params, across_beat_params, debug=FALSE, renal_param, sys_t){
  #################################################################################################################################
  # model2.FIX_PAR3 takes a set of 12 parameters as input and fixes any that exceed constraints defined within the function.
  # For every violation of model constraints a penalty is accrued, the sum of which is outputted along with the fixed parameters.

  # Inputs:
  # time (x values of the ppg segment being modelled)
  # within_beat_params (parameters not fixed across beats)
  # across_beat_params (parameters fixed across beats)
  # debug (logical, if true function prints parameter values relative to constraint values, and the penalties accrued)
  # renal_param (the timing of the 2nd reflectance wave, used as a reference to the initial estimation for this parameter)
  # sys_t (the timing of the systolic wave, used as a reference to the initial estimation for this parameter)

  # Output:
  # result (vector of fixed parameters and total penalty accrued)

  # Notes:
  # model2.FIX_PAR3 can handle parameter inputs that indicate only a single systolic wave, or two component waves only. Likewise
  # it can handle one rather than two baseline parameter inputs. This is to support modelling of waveforms with only as many
  # parameters as are necessary (hence also the use of reduced ChiSq as the chosen measure of goodness of fit). The optionality
  # of parameter inputs is summarized as follows:

  # params: {Baseline, {baseline 2}, t_sys, H_sys, W_sys, {dt_1, H_1, W_1, {dt_2, H_2, W_2}}}
  # across_beat_params: { w[1], t[2], w[2], t[3], w[3] }

  # Or:
  # par: {base1, {base2}, t[1], h[1], #, #, { h[2], #, #, { h[3], ..., ... }}}
  ################################################################################################################################


  nPar <- length(within_beat_params) + length(across_beat_params)              # Establish number of inputted parameters
  nData <- length(time)

  nBase <- 1
  baseline <- c(within_beat_params[1], within_beat_params[1])                  # Establish number of baselines
  if (nPar == 6 | nPar == 9 | nPar == 12){
    baseline[2] = within_beat_params[2]
    nBase <- 2
  }
  # Arrange peak parameters into vectors of timing, height and width:

  t <- c( within_beat_params[nBase + 1],                                       # timing (systolic = within, diastolic / renal (R1, R2) = across)
          across_beat_params[2], across_beat_params[4])
  h <- c( within_beat_params[nBase + 2], 0, 0 )                                # height (systolic = within, diastolic / renal (R1, R2) default to 0 unless peaks supplied (see below))
  w <- c( across_beat_params[1],
          across_beat_params[3], across_beat_params[5])                        # width (systolic / diastolic (R1) / renal (R2) = across)


  hasPeak <- c(TRUE, FALSE, FALSE)                                             # Assume initially that only one peak is supplied



  if (nPar >= nBase + 7){                                                      # If parameters for a second peak are available, add its parameters
    hasPeak[2] = TRUE
    t[2] <- across_beat_params[2]
    h[2] <- within_beat_params[nBase + 3]
    w[2] <- across_beat_params[3]
  }

  if (nPar >= nBase + 10){                                                     # Assign parameters for a third peak are available, add its parameters
    hasPeak[3] = TRUE
    t[3] <- across_beat_params[4]
    h[3] <- within_beat_params[nBase + 4]
    w[3] <- across_beat_params[5]
  }

  penalty <- 0                                                                 # Define penalty (as 0 initially)

  tMin <- time[1]                                                              # Define earliest and latest time points in the PPG segment
  tMax <- time[length(time)]

  META_BASELINE_SHIFT <- 1.0                                                   # Define the constraint for difference in y-value between baselines
  META_MIN_PEAK_DELAY <- 0.1                                                   # Define the constraint for minimum timing between peaks (peaks cannot follow one another by less than 0.1ms)
  MIN_WIDTH <- c(0.05, 0.05, 0.1)                                              # Define the constaints for (minimum) width, for each of the three component waves
  MAX_WIDTH <- c(0.5, 0.45, 0.25)                                              # Define the constaints for (maximum) width, for each of the three component waves

  p <- 1:12*0                                                                  # Define also a penalty vector, each value corresponding to 1 of the 12 parameters as follows:
  # p: { #, #, t[1], h[1], w[1], t[2], h[2], w[2], t[3], h[3], w[3], across_beat_params[6] }

  # Penalizing and fixing of component wave parameter values (amplitude and width):

  for(i in 1:3){                                                               # For each component wave, penalize and fix amplitudes below 0
    if (h[i] < 0){
      penalty <- penalty + h[i]*h[i]
      p[3*i+1] <- h[i]*h[i]
      h[i] <- 0
    }

    if (w[i] < MIN_WIDTH[i] | w[i] > MAX_WIDTH[i]){                            # For each component wave, penalize and fix widths outside of defined constraints
      fixed <- max(MIN_WIDTH[i], min( w[i], MAX_WIDTH[i]))
      diff <- fixed - w[i]
      penalty <- penalty + diff*diff
      p[3*i+2] <- diff*diff
      w[i] <- fixed
    }

    if(i==3){                                                                  # For the 1st reflectance wave, penalize amplitudes exceeding 2% of the systolic amplitude
      if( h[3] > (h[1]/50)){
        diff <- h[3] - (h[1]/50)
        penalty <- penalty + 2*diff*diff
        p[10] <- p[10] + 2*diff*diff
      }
    }
  }

  # Penalizing and fixing of component wave parameter values (timing):

  fixed <- max((sys_t - 0.04) , min( t[1], (sys_t + 0.04 )))                   # Penalise and fix systolic waves occuring outside of 40ms of the peak of the ppg data segment
  if (debug){
    print(paste("time S: ",tMin," < ",t[1],
                " < min( ",tMin+1,",",tMax," )"))
  }
  if (t[1] != fixed){
    diff <- fixed - t[1]
    penalty <- penalty + 10^8*diff*diff
    p[3] <- 10^8*diff*diff
    t[1] <- fixed
  }

  fixed <- max(2 * META_MIN_PEAK_DELAY,                                        # Penalize and fix diastolic (2nd reflectance wave) timing if < 0.2 seconds after
               min( t[2], tMax - tMin + 0.4 * w[2]))                           # systolic wave, or otherwise occurring too early in the segment
  if (debug){
    print(paste("time D: ",2 * META_MIN_PEAK_DELAY," < ",
                t[2]," < ",tMax - tMin + 0.4 * w[2]," )"))
  }
  if (t[2] != fixed){
    diff <- fixed - t[2]
    if (hasPeak[2]){
      penalty <- penalty + diff*diff
      p[6] <- diff*diff
    }
    t[2] <- fixed
  }

  fixed <- max( max(META_MIN_PEAK_DELAY, renal_param - 0.02),                  # Penalise and fix renal (1st reflectance wave) timing if < 0.1 seconds after
                min(t[3], t[2] - META_MIN_PEAK_DELAY, renal_param + 0.02))     # systolic wave, or < 0.1 seconds before the 2nd reflectance wave, or outside 20ms
  if (debug){                                                                  # of the initially estimated value for the same parameter.
    print(paste("time R: ",META_MIN_PEAK_DELAY,
                " < ",t[3]," < ",t[2] - META_MIN_PEAK_DELAY," )"))
  }
  if (t[3] != fixed){
    diff <- fixed - t[3]
    if (hasPeak[3]){
      penalty <- penalty + 5*10^7*diff*diff
      p[9] <- diff*diff
    }
    t[3] <- renal_param
  }

  # Penalizing and fixing of decay element parameters:

  if(across_beat_params[6] > 0.95){                                            # Penalise and fix the decay rate parameter if > 0.95 (indicating less than 5% decay over the time interval of a single data point)
    diff <- across_beat_params[6] - 0.95
    penalty <- penalty + 10^7*diff*diff
    p[12] <- 10^7*diff*diff
    across_beat_params[6] <- 0.95
  }

  if(baseline[1] > 0){                                                         # Penalize and fix Baseline1 if its y-value is above 0
    penalty <- penalty + baseline[1]*baseline[1]
    p[1] <- baseline[1]*baseline[1]
    baseline[1] <- 0
  }

  fixedPar <- c(baseline[1:2], t[1], h[1], w[1], t[2],                         # Combine all fixed parameters into a single vector for outputting
                h[2], w[2], t[3], h[3], w[3], across_beat_params[6])

  if(debug){
    print(p)
  }

  return(c(penalty, fixedPar))
}


model2.FixParams3 <- function(data,params, across_beat_params = NULL, debug=FALSE, rp = renal_param, sys_t){
  #################################################################################################################################
  # model2.FixParams3 organises inputs and outputs of the model2.FIX_PAR3 function.

  # Inputs:
  # data (ppg time series segment relevant to waveform being modelled)
  # params (vector of model parameters)
  # across_beat_params (vector of fixed parameters (set to null if all parameters to be sourced from params input))
  # debug (logical, now redundant)
  # rp (1st reflectance wave timing paramter intiial estimation)
  # sys_t (systolic wave timing paramter intiial estimation)

  # Output:
  # temp (vector of fixed parameters)
  ################################################################################################################################
  if(is.null(across_beat_params)){
    across_beat_params <- params[c(5, 6, 8, 9, 11, 12)]                        # If across_beat_params have not been provided, extract from params
  }

  temp <- model2.FIX_PAR3(time = data[, 1],                                    # Pass inputs into model2.FIX_PAR3 function
                          within_beat_params = params[c(1:4, 7, 10)],
                          across_beat_params, debug = F, renal_param = rp,
                          sys_t)
  return( temp[2:length(temp)] )                                               # Return all outputs from model2.FIX_PAR3 except penalty
}




simplex.MakeSimplex2 <- function(data,param,f,inScale,directions=NULL,inTol=-1, optional=NULL,debug=FALSE, beat_vector = beat_vector,
                                 beat = beat, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w){
  ########################################################################################################################################
  # simplex.MakeSimplex2 iterates on the initially estimated values for fixed parameters, refining them to give the simplex an improved
  # starting position.

  # Inputs:
  # data (ppg time series segment)
  # param (model parameter inputs)
  # f (model2.ChiSq3 function)
  # inScale (the starting value by which to iterate on parameter values (0.1 recommended))
  # directions (the initial direction in which to iterate parameter values (default null))
  # inTol (the minimum value a parameter can change if it cannot be iterated upon to improve goodness of fit)
  # optional (for direct input into model2.ChiSq3 function, default NULL)
  # debug (logical, if true iteration and refinement statuses will be printed to console)
  # beat_vector (an index of beats to be modeled and their x-coordinates to extract from the PPG time series (for direct input to f))
  # beat (matrix of intial parameter estimates)
  # renal_param (the starting parameter for 1st reflectance peak timing (inputted to prevent drastic deviations from this value))
  # dias_param (the starting parameter for 2nd reflectance peak timing (inputted to prevent drastic deviations from this value))
  # sys_time (the starting parameter for systolic peak timing (inputted to prevent drastic deviations from this value))
  # w (the timing of the 1st derivative peaks on the ppg time series)

  # Outputs:
  # result (a matrix of parameter sets, including the initially inputted parameters (row 1) and a parameter set for each refined parameter)

  # Notes:
  # simplex.MakeSimplex2 is for fixed (across beat) parameters only. Thus the outputted matrix will include absent (NA) rows where non-fixed
  # parameters are represented (these are refined separately in the simplex.MakeSimplex3 function).
  ########################################################################################################################################

  if(debug){print("MakeSimplex -- debug")}                                     # Indicate if in debug mode

  nPar <- length(param)                                                        # Define number of parameters
  nScale <- length(inScale)


  if (nScale == 0)                                                             # Set scale values according to number of parameters
  {
    scale <- 1:nPar * 0 + 1
  } else if (nScale == 1){
    scale <- 1:nPar * 0 + inScale
  } else if (length(inScale) == nPar){
    scale <- inScale
  } else {
    return("Error: Invalid scale vector length")                               # [DEBUG NOTE] print("Invalid scale vector length") could be added here
  }


  if (length(inTol) == 1 & inTol > 0){                                         # If no tolerance is provided, tol is defined as the reduced ChiSq value of the inputted parameter set
    tol <- inTol[1]
  } else {
    tol <- min(1,f(data,param, beats = beat_vector, beat = beat,
                   renal_param = renal_param,
                   dias_param = dias_param,
                   sys_time = sys_time, w = w))
  }


  chiSq <- 1:(nPar+1) * 0.0                                                    # Create a vector of ChiSq values (to store goodness of fit values for each parameter)

  chiSq[1] <- f(data,param,optional=optional, beats = beat_vector,             # ChiSq[1] is the goodness of fit when no parameters are changed
                beat = beat, renal_param = renal_param,
                dias_param = dias_param, sys_time = sys_time, w = w)

  if(debug){print(paste("Root chi-squared:",chiSq[1]))}

  result <- matrix(nrow=nPar+1,ncol=nPar)                                      # Create a matrix to be outputted by the function, with as many columns as parameters, and
  # the same number of rows + 1 (row 1 will be the initially estimated parameters, whilst all
  # other rows will represent parameter sets where all parameters remain the same as row 1 with
  # the exception of one paramter that is refined as per the below procedure.

  result[1,] <- as.double(param)                                               # Fill the first row with the inputted parameters

  useDirections = !is.null(directions)                                         # [DEBUG NOTE] The directions input is currently redundant
  if(useDirections){
    useDirections <- nrow(directions) == nPar & ncol(directions) == nPar
  }

  for (i in c(5, 6, 8, 9, 11, 12)){                                            # Refine the values of each fixed parameter in turn:

    if(debug){print(paste("Parameter",i))}

    tParam <- param                                                            # Create a new vector of test parameters to be tested for goodness of fit, tParam

    delta <- 1:nPar * 0                                                        # Create a vector for delta, the finite increment by which to increase or decrease a given parameter value

    if(useDirections){                                                         # Give a value to the element of the delta vector corresponding to the parameter being refined
      delta <- scale[i] * directions[i,]
    }else{
      delta[i] <- scale[i]                                                     # [DEBUG NOTE] scale relates to inScale, and is used to set the delta level
    }

    tParam <- param - delta                                                    # Differentiate test parameters from inputted parameters by subtracting the delta value from the parameter being refined

    chiSqMinus <- f(data,tParam,optional=optional,                             # Calculate the goodness of fit of the new test parameters when the delta increment is subtracted
                    beats = beat_vector, beat = beat,
                    renal_param = renal_param,
                    dias_param = dias_param, sys_time = sys_time, w = w)

    tParam <- param + delta                                                    # Redefine the test parameters by this time adding rather than subtracting the delta value

    chiSq[i+1] <- f(data,tParam,optional=optional,                             # Calculate the goodness of fit of the test parameters when the delta increment is added
                    beats = beat_vector, beat = beat,
                    renal_param = renal_param,
                    dias_param = dias_param, sys_time = sys_time, w = w)


    if (debug){
      print("Select direction:")
      print(paste("chi^2(",param[i] - delta[i],") =",chiSqMinus))
      print(paste("chi^2(",param[i],") =",chiSq[1]))
      print(paste("chi^2(",param[i] + delta[i],") =",chiSq[i+1]))
      print("---")
    }

    if (chiSqMinus < chiSq[i+1]){                                              # If subtracting delta yields a better goodness of fit than adding delta,
      delta <- -delta                                                          # then invert delta, redefine the test parameters to have delta subtracted,
      tParam <- param + delta                                                  # and define the goodness of fit as the lower ChiSq value. The direction in
      chiSq[i+1] <- chiSqMinus                                                 # which to the refine the ith parameter value (increase or decrease) is now decided.
    }

    iKill <- 10                                                                # Define the number of iterations that a given parameter will be refined over

    if (chiSq[i+1] < chiSq[1]){                                                # If the new fit is better than the old fit (with no parameters changed),
      if (debug){ print("Extending as best point") }                           # extend the value in the same direction so long as fit continues to improve
      while (chiSq[i+1] < chiSq[1] + tol){
        delta <- 2*delta
        tParam <- param + delta
        oldScore <- chiSq[i+1]                                                 # During each iteration, the current best fit is defined as 'old score'
        chiSq[i+1] <- f(data,tParam,optional=optional,                         # while the new fit (with the ith parameter value further extended) is designated ChiSq[i+1]
                        beats = beat_vector, beat = beat,
                        renal_param = renal_param,
                        dias_param = dias_param,
                        sys_time = sys_time, w = w)
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }
        if (chiSq[i+1] > oldScore){                                            # During each iteration, check if the iterated parameter value results in an improved or worsened fit
          tParam <- param + 0.5*delta
          chiSq[i+1] <- oldScore                                               # If an iteration results in a worsened fit, redesignate 'old score' to chiSq[i+1]
          break
        }
        iKill <- iKill - 1                                                     # If an iteration results in an improved fit, continue to iterate until further iteration does not
        if (iKill < 0){                                                        # improve fit or until the number of defined iterations is reached
          break
        }
      }
    } else if (chiSq[i+1] < chiSq[1] + tol){                                   # If the new fit is not better than the old fit (with no parameters changed),
      if (debug){ print("Extending below tolerance") }                         # check if it is at least better than the old fit + the defined tolerance level.
      while (chiSq[i+1] < chiSq[1] + tol){
        delta <- 2*delta
        tParam <- param + delta
        oldScore <- chiSq[i+1]
        chiSq[i+1] <- f(data,tParam,optional=optional,                         # If so, continue to iterate, 'extending below tolerance', until fit can not be improved
                        beats = beat_vector, beat = beat,                      # or the number of defined iterations is reached.
                        renal_param = renal_param,
                        dias_param = dias_param,
                        sys_time = sys_time, w = w)
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }
        if (chiSq[i+1] - oldScore < oldScore - chiSq[1]){
          tParam <- param + 0.5*delta
          chiSq[i+1] <- oldScore
          break
        }
        iKill <- iKill - 1
        if (iKill < 0){
          if(i == 9){
            tParam[9] <- renal_param                                           # [DEBUG NOTE]: Ignore reflectance wave one times that can't optomize
            break
          }
          print(c("simplex constructed as per original parameter"))
          break                                                                # [DEBUG NOTE]: May be useful to include # print("Failed to construct simplex") # return(paste("Error: param[",i,"]",sep=""))
        }
      }
    } else {
      if (debug){ print("Shrinking above tolerance") }                         # If the new fit is not better than the old fit (with no parameters changed) nor the old fit + the
      while (chiSq[i+1] > chiSq[1] + tol){                                     # defined tolerance level, iterate in the same direction but with a smaller delta value (shrinking above tolerance)
        delta <- 0.5*delta
        tParam <- param + delta
        lastChiSq <- chiSq[i+1]
        chiSq[i+1] <- f(data,tParam,optional=optional,
                        beats = beat_vector, beat = beat,
                        renal_param = renal_param,
                        dias_param = dias_param,
                        sys_time = sys_time, w = w)
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }
        if (iKill < 0 & (chiSq[i+1]-chiSq[1]) > 0.75 * (lastChiSq-chiSq[1])){
          if(i == 9){
            tParam[9] <- renal_param
            next
          }
          print(c("simplex constructed as per original parameter"))
          next
        }
        iKill <- iKill - 1
      }
      tParam <- param + 0.5 * delta
    }

    if(debug){ print(paste("Param[",i,"] =",tParam[i]))}
    result[i+1,] = as.double(tParam)                                           # The value of the ith parameter providing the best goodness of fit (with all other parameters
    # held constant) defines a unique parameter set that become the (i+1)th row in the output matrix
  }

  if (debug){ print("/MakeSimplex") }
  return(result)
}




simplex.MakeSimplex3 <- function(ppg, param,f,inScale, directions=NULL, inTol=-1, optional=NULL, debug=FALSE, beat_vector = beat_vector, renal_param, dias_param = dias_param, sys_time, w){
  ########################################################################################################################################
  # simplex.MakeSimplex3 iterates on the initially estimated values for non-fixed parameters, refining them to give the simplex an improved
  # starting position. It is used within the FindWithinParams function, such that a matrix of optimised non-fixed parameters can be generated
  # for each beat in an inputted batch.

  # Inputs:
  # ppg (ppg time series segment)
  # param (model parameter inputs)
  # f (model2.ChiSq3 function)
  # inScale (the starting value by which to iterate on parameter values (0.1 recommended))
  # directions (the initial direction in which to iterate parameter values (default null))
  # inTol (the minimum value a parameter can change if it cannot be iterated upon to improve goodness of fit)
  # optional (for direct input into model2.ChiSq3 function, default NULL)
  # debug (logical, if true iteration and refinement statuses will be printed to console)
  # beat_vector (an index of beats to be modeled and their x-coordinates to extract from the PPG time series (for direct input to f))
  # renal_param (the starting parameter for 1st reflectance peak timing (inputted to prevent drastic deviations from this value))
  # dias_param (the starting parameter for 2nd reflectance peak timing (inputted to prevent drastic deviations from this value))
  # sys_time (the starting parameter for systolic peak timing (inputted to prevent drastic deviations from this value))
  # w (the timing of the 1st derivative peaks on the ppg time series)

  # Outputs:
  # result (a matrix of parameter sets, including the initially inputted parameters (row 1) and a parameter set for each refined parameter)

  # Notes:
  # simplex.MakeSimplex3 is for non-fixed (within beat) parameters only. Thus the outputted matrix will include absent (NA) rows where non-fixed
  # parameters are represented (these are refined separately in the simplex.MakeSimplex2 function).
  ########################################################################################################################################

  if(debug){print("MakeSimplex -- debug")}                                     # Indicate if in debug mode

  nPar <- length(param)                                                        # Define number of parameters
  nScale <- length(inScale)

  if (nScale == 0)                                                             # Set scale values according to number of parameters
  {
    scale <- 1:nPar * 0 + 1
  } else if (nScale == 1){
    scale <- 1:nPar * 0 + inScale
  } else if (length(inScale) == nPar){
    scale <- inScale
  } else {
    return("Error: Invalid scale vector length")
  }

  if (length(inTol) == 1 & inTol > 0){                                         # If no tolerance is provided, tol is defined as the reduced ChiSq value of the inputted parameter set
    tol <- inTol[1]
  } else {
    tol <- min(1,f(data = ppg, params = param, beats = beat_vector,
                   renal_param = renal_param, dias_param = dias_param,
                   sys_time = sys_time, w = w))
  }


  chiSq <- 1:(nPar+1) * 0.0                                                    # Create a vector of ChiSq values (to store goodness of fit values for each parameter)

  chiSq[1] <- f(data = ppg, param, beats = beat_vector,                        # ChiSq[1] is the goodness of fit when no parameters are changed
                renal_param = renal_param,
                dias_param = dias_param,
                sys_time = sys_time, w = w)

  if (debug){ print(paste("Root chi-squared:",chiSq[1]))}

  result <- matrix(nrow=nPar+1,ncol=nPar)                                      # Create a matrix to be outputted by the function
  result[1,] <- as.double(param)                                               # Fill the first row with the inputted parameters

  useDirections = !is.null(directions)
  if (useDirections){
    useDirections <- nrow(directions) == nPar & ncol(directions) == nPar
  }

  for(i in c(1:4, 7, 10)){                                                     # Refine the value of each non-fixed parameter in turn:

    if (debug){ print(paste("Parameter",i))}

    tParam <- param                                                            # Create a new vector of test parameters

    delta <- 1:nPar * 0                                                        # Create a vector for delta, the finite increment by which to increase or decrease a given parameter value

    if (useDirections){                                                        # Give a value to the element of the delta vector corresponding to the parameter being refined
      delta <- scale[i] * directions[i,]
    } else {
      delta[i] <- scale[i]
    }

    if(i == 3){                                                                # Specify a unique delta value when refining systolic timing parameter, which often requires smaller adjustments
      delta <- delta/4
    }

    tParam <- param - delta                                                    # Differentiate test parameters from inputted parameters by subtracting the delta value from the parameter being refined

    chiSqMinus <- f(data = ppg, params = tParam,                               # Calculate the goodness of fit of the new test parameters when the delta increment is subtracted
                    beats = beat_vector,
                    renal_param = renal_param,
                    dias_param = dias_param,
                    sys_time = sys_time, w = w)

    tParam <- param + delta                                                    # Redefine the test parameters by this time adding rather than subtracting the delta value

    chiSq[i+1] <- f(data = ppg, params = tParam,                               # Calculate the goodness of fit of the test parameters when the delta increment is added
                    beats = beat_vector,
                    renal_param = renal_param,
                    dias_param = dias_param,
                    sys_time = sys_time, w = w)

    if (debug){
      print("Select direction:")
      print(paste("chi^2(",param[i] - delta[i],") =",chiSqMinus))
      print(paste("chi^2(",param[i],") =",chiSq[1]))
      print(paste("chi^2(",param[i] + delta[i],") =",chiSq[i+1]))
      print("---")
    }

    if (chiSqMinus < chiSq[i+1]){                                              # If subtracting delta yields a better goodness of fit than adding delta,
      delta <- -delta                                                          # then invert delta, redefine the test parameters to have delta subtracted,
      tParam <- param + delta                                                  # and define the goodness of fit as the lower ChiSq value. The direction in
      chiSq[i+1] <- chiSqMinus                                                 # which to the refine the ith parameter value (increase or decrease) is now decided.
    }

    iKill <- 10                                                                # Define the number of iterations that a given parameter will be refined over


    if (chiSq[i+1] < chiSq[1]){                                                # If the new fit is better than the old fit (with no parameters changed),
      if (debug){ print("Extending as best point") }                           # extend the value in the same direction so long as fit continues to improve.
      while (chiSq[i+1] < chiSq[1] + tol){
        delta <- 2*delta
        tParam <- param + delta
        oldScore <- chiSq[i+1]                                                 # During each iteration, the current best fit is defined as 'old score'
        chiSq[i+1] <- f(data = ppg,tParam,                                     # while the new fit (with the ith parameter value further extended) is designated ChiSq[i+1]
                        beats = beat_vector,
                        renal_param = renal_param,
                        dias_param = dias_param,
                        sys_time = sys_time, w = w)
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1]))}
        if (chiSq[i+1] > oldScore){                                            # During each iteration, check if the iterated parameter value results in an improved or worsened fit
          tParam <- param + 0.5*delta
          chiSq[i+1] <- oldScore                                               # If an iteration results in a worsened fit, redesignate 'old score' to chiSq[i+1]
          break
        }
        iKill <- iKill - 1                                                     # If an iteration results in an improved fit, continue to iterate until further iteration does not
        if (iKill < 0){                                                        # improve fit or until the number of defined iterations is reached
          break
        }
      }
    } else if (chiSq[i+1] < chiSq[1] + tol){                                   # If the new fit is not better than the old fit (with no parameters changed),
      if (debug){print("Extending below tolerance")}                           # check if it is at least better than the old fit + the defined tolerance level.
      while (chiSq[i+1] < chiSq[1] + tol){
        delta <- 2*delta
        tParam <- param + delta
        oldScore <- chiSq[i+1]
        chiSq[i+1] <- f(data = ppg, tParam,                                    # If so, continue to iterate, 'extending below tolerance', until fit can not be improved
                        beats = beat_vector,                                   # or the number of defined iterations is reached.
                        renal_param = renal_param,
                        dias_param = dias_param,
                        sys_time = sys_time, w = w)
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }
        if (chiSq[i+1] - oldScore < oldScore - chiSq[1]){
          tParam <- param + 0.5*delta
          chiSq[i+1] <- oldScore
          break
        }
        iKill <- iKill - 1                                                     # If 10 iterations do not improve the parameter, default to initial starting estimation
        if (iKill < 0){
          print(c("Failed to construct simplex within 10 iterations for parameter",
                  i, "defaulting to inputted value"))
          tParam[i] <- param[i]
          break
        }
      }
    } else {
      if (debug){ print("Shrinking above tolerance")}                          # If the new fit is not better than the old fit (with no parameters changed) nor the old fit + the
      while (chiSq[i+1] > chiSq[1] + tol){                                     # defined tolerance level, iterate in the same direction but with a smaller delta value (shrinking above tolerance)
        delta <- 0.5*delta
        tParam <- param + delta
        lastChiSq <- chiSq[i+1]
        chiSq[i+1] <- f(data = ppg,tParam,
                        beats = beat_vector,
                        renal_param = renal_param,
                        dias_param = dias_param,
                        sys_time = sys_time, w = w)
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }
        if (iKill < 0 & (chiSq[i+1]-chiSq[1]) > 0.75 * (lastChiSq-chiSq[1])){
          print(c("Failed to construct simplex within 10 iterations for parameter",
                  i, "defaulting to inputted value"))
          tParam[i] <- param[i]
          break
        }
        iKill <- iKill - 1
      }
      tParam <- param + 0.5 * delta
    }

    if(debug){ print(paste("Param[",i,"] =",tParam[i]))}                       # The value of the ith parameter providing the best goodness of fit (with all other parameters
    result[i+1,] = as.double(tParam)                                           # held constant) defines a unique parameter set that become the (i+1)th row in the output matrix
  }

  if (debug){ print("/MakeSimplex") }
  return(result)
}



simplex.Run2 <- function(data = ppg,simplexParam = mat, f = model2.ChiSq3, optional=NULL, beat_vector = beat_vector, ms = simplex_iterations, renal_param = renal_param, dias_param = dias_param, sys_time = sys_time, w = w, run = NULL){
  ########################################################################################################################################
  # simplex.Run2 intiates the Nelder Mead downhill simplex routine to optimise model parameters and goodness of fit.

  # Inputs:
  # data (ppg time series segment)
  # simplexParam (the matrix outputted by make.matrix, representative of all model parameters in an inputted batch of beats)
  # f (model2.ChiSq3 function)
  # optional (currently redundant)
  # beat_vector (an index of beats to be modeled and their x-coordinates to extract from the PPG time series)
  # ms (number of simplex iterations to conduct)
  # renal_param (the starting parameter for 1st reflectance peak timing (inputted to prevent drastic deviations from this value))
  # dias_param (the starting parameter for 2nd reflectance peak timing (inputted to prevent drastic deviations from this value))
  # sys_time (the starting parameter for systolic peak timing (inputted to prevent drastic deviations from this value))
  # w (the timing of the 1st derivative peaks on the ppg time series)
  # run (indicates which iteration of the simplex routine is being run)

  # Outputs:
  # A matrix with equivalent structure to the one inputted, with optimised parameter values.
  ########################################################################################################################################

  MAX_STEP <- ms                                                               # Maximum number of allowed simplex iterations
  FTOL <- 1e-5                                                                 # The minimum 'height' difference between highest and lowest simplex vertices to warrant further iteration

  debugRtol <- 1:(MAX_STEP+1) * 0.0                                            # [DEBUG NOTE]: These vectors may be redundant
  debugMin <- 1:(MAX_STEP+1) * 0.0
  debugMax <- 1:(MAX_STEP+1) * 0.0

  result <- simplexParam                                                       # Defining the simplex matrix as inputted
  nPar <- ncol(result)
  chiSq <- 0:nPar * 0.0                                                        # Create a vector of goodness of fit values for each parameter set in the matrix

  for (i in 1:(nPar+1)){                                                       # Populate said vector with values by assessing goodness of fit for each row / parameter set
    chiSq[i] <- f(data, params = NULL, optional=NULL,
                  a = result[i, ], beats = beat_vector,
                  renal_param = renal_param,                                   # Each value of the vector represents the 'height' of one vertex of the simplex in multi-dimensional space
                  dias_param = dias_param,
                  sys_time = sys_time, w = w)
  }

  for (iStep in 1:MAX_STEP){                                                   # Begin downhill simplex routine:
    extrema <- simplex.SortHighLow(chiSq)                                      # Calculate the vertices (parameter sets) with the highest, 2nd highest and lowest ChiSq values (expressed as row numbers)
    low <- extrema[1]
    nHigh <- extrema[2]
    high <- extrema[3]

    if(!is.null(run)){                                                         # Print run number and iteration number
      print(run)
    }
    print(iStep)

    chiSqMax <- chiSq[high]                                                    # Redefine highest and lowest vertices of the simplex
    chiSqMin <- chiSq[low]
    print(chiSqMax)                                                            # [DEBUG NOTE]: may be informative to run # print(paste("chi^2_min =",chiSqMin)) # print(paste("argMax = ",high,"[",chiSqMax,"]",sep=""))


    rtol <- 2 * (chiSqMax - chiSqMin)/(chiSqMax + chiSqMin + 1e-10)            # Define the relative heights of lowest and heighest vertices as 'rtol'
    if (rtol < FTOL){                                                          # [DEBUG NOTE]: if issues with simplex not terminating regardless of number of iterations, may be worth reconsidering FTOL / 1e-10 value
      bestParam <- result[low,]
      result[low,] <- result[1,]                                               # If rtol is smaller than FTOL, this indicates that further simplex iterations are likely to yield minimal improvements in goodness of fit.
      result[1,] <- bestParam                                                  # Therefore, the routine is terminated: the minimum vertex is considered the final parameter set and the matrix is outputted.
      return(result)
    }

    debugRtol[iStep] <- rtol                                                   # [DEBUG NOTE]: These vectors may be redundant
    debugMin[iStep] <- chiSqMin
    debugMax[iStep] <- chiSqMax

    factor <- -1                                                               # Define the degree to which the simplex reflects away from it's highest vertex as 'factor'
    node <- simplex.HypoCentre(result,high)                                    # simplex.HypoCentre outputs a single vertex (parameter set) representing the centre of the simplex when excluding the highest vertex
    apex <- result[high,]                                                      # Apex is the highest vertex on the simplex
    test <- node - (apex - node)                                               # The 'test' vertex represents the apex vertex reflected about the node vertex (akin to flipping one point of a triangle through the middle of the base in 2D)
    score <- f(data, params = rep(0, 12),                                      # The newly reflected test vertex is evaluated for goodness of fit (the downhill simplex routine works on the premise that reflecting the highest vertex about
               optional = optional,                                            # the node will result in a new 'lowest' vertex allowing the simplex to reach continually decreasing minimum vertex values.
               a = test, beats = beat_vector,
               renal_param = renal_param,
               dias_param = dias_param,
               sys_time = sys_time, w = w)

    if (score < chiSqMin){                                                     # If reflecting the apex improves goodness of fit, try extending further in the same direction (increase the scale by 2) i.e reflection and expansion
      test2 <- node - 2 * (apex - node)
      score2 <- f(data, params = rep(0, 12),                                   # Calculate the goodness of fit for the reflected and extended vertex
                  optional = optional,
                  a = test2, beats = beat_vector,
                  renal_param = renal_param,
                  dias_param = dias_param,
                  sys_time = sys_time, w = w)
      if (score2 >= score){                                                    # If reflecting a further distance does not result in a better fit, revert to reflection only and redefine the simplex
        result[high,] <- test                                                  # [DEBUG NOTE]: print(paste("Reflecting",high,": chi^2 ",chiSqMax,"->",score,sep="")) may be useful to assess simplex behaviour
        chiSq[high] <- score
      } else {                                                                 # If reflecting a further distance does result in a better fit, redefine the simplex in its reflected and extended form
        result[high,] <- test2                                                 # [DEBUG NOTE]: print(paste("Reflect-stretching",high,": chi^2 ",chiSqMax,"->",score2,sep="")) may be useful
        chiSq[high] <- score2
      }
    } else if (score >= chiSq[nHigh]) {                                        # If reflecting is not beneficial, try shrinking away from the apex instead of reflecting away.
      factor <- 0.5                                                            # Factor is redefined as a fraction, indicating shrinkage
      if (score < chiSqMax)
      {
        factor <- -0.5
      }
      test2 <- node + factor * (apex - node)
      score2 <- f(data, params = rep(0, 12),                                   # Calculate goodness of fit of 'shrunken apex' parameter set
                  optional=optional, a = test2,
                  beats = beat_vector,
                  renal_param = renal_param,
                  dias_param = dias_param,
                  sys_time = sys_time, w = w)
      if (score2 < chiSq[nHigh]){                                              # If shrunken apex represents an improved goodness of fit, redefine the simplex to include it  # print(paste("Shrinking",high,": chi^2 ",chiSqMax,"->",score2,sep=""))
        result[high,] <- test2
        chiSq[high] <- score2
      } else {                                                                 # If shrinking from the apex does not result in an improved fit, contract instead towards the lowest vertex # print(paste("General contraction: chi^2 ",chiSqMax,"->",max(chiSq),sep=""))
        for (i in 1:(nPar+1)){
          if (i != low){
            result[i,] <- 0.5 * (result[i,] + result[low,])
            chiSq[i] <- f(data, params = rep(0, 12),
                          optional = optional,
                          a = result[i, ],
                          beats = beat_vector,
                          renal_param = renal_param,
                          dias_param = dias_param,
                          sys_time = sys_time, w = w)
          }
        }
      }
    } else {                                                                   # If the reflected apex represented an improved goodness of fit but did not become the new minimum of the simplex, reflect anyway
      result[high,] <- test                                                    # print(paste("Reflecting*",high,": chi^2 ",chiSqMax,"->",score,sep=""))
      chiSq[high] <- score
    }
  }                                                                            # Continue to iterate until FTOL or number of iterations is reached

  extrema <- simplex.SortHighLow(chiSq)                                        # Define new lowest vertex and move this to the top row of the matrix, representing the best parameter set
  low <- extrema[1]
  bestParam <- result[low,]
  result[low,] <- result[1,]
  result[1,] <- bestParam

  chiSqMax <- chiSq[extrema[3]]
  chiSqMin <- chiSq[low]

  rtol <- 2 * (chiSqMax - chiSqMin)/(chiSqMax + chiSqMin + 1e-10)              # [DEBUG NOTE]: rtol and debug variabets may be redundant
  debugRtol[MAX_STEP+1] <- rtol                                                # plot(debugMax,type='l') # lines(debugMin)
  debugMin[MAX_STEP+1] <- chiSqMin
  debugMax[MAX_STEP+1] <- chiSqMax

  print(paste("Terminated downhill simplex after",MAX_STEP,"iterations."))     # Print to indicate termination of downhill simplex
  print(paste("rtol =",rtol))
  return(result)
}


simplex.HypoCentre <- function(mat_Param, index){
  #################################################################################################################################
  # simplex.HypoCentre identifies the average vertex of the simplex when the apex (highest point / parameter set yielding worst
  # goodness of fit) is excluded. This is referred to as the 'node' or 'hypocentre' and represents the point about which the
  # worst vertex is reflected in order to generate a new and lower vertex.

  # Inputs:
  # mat_Param (a simplex matrix)
  # index (the row of the inputted matrix corresponding to the worst (poorest fit) parameter set)

  # Output:
  # result (a vertex representing the mean of simplex vertices when the apex is excluded (equivalent to a triangle's base))
  ################################################################################################################################
  nPar <- ncol(mat_Param)                                                      # Number of parameters defined

  result <- 1:nPar * 0.0                                                       # Define a vector of length nPar
  for (i in 1:(nPar+1)){
    if (i != index){                                                           # For each row of the matrix, if the ith row is not the worst row, add the row to the result vector
      result <- result + mat_Param[i,]                                         # The end result is a vector representing the sum of all parameter sets excluding the worst row
    }
  }
  return( result / nPar )                                                      # Return the summed vector once divided by number of parameters, to give the average (mean) vertex (excluding the worst)
}


simplex.SortHighLow <- function(vec_ChiSq){
  #################################################################################################################################
  # simplex.SortHighLow calculates the vertices (parameter sets) of a simplex with the highest, second highest and lowest ChiSq
  # values (expressed as row numbers of the simplex matrix).

  # Inputs:
  # vec_ChiSq (a vector of goodness of fit (ChiSq) values corresponding to each row (/parameter set) of the simplex)

  # Output:
  # low,nHigh,high (vector of the lowest, second highest, and highest ChiSq values in the simplex)
  ################################################################################################################################
  nPar <- length(vec_ChiSq)                                                    # Define number of parameters

  low <- 1                                                                     # Assign initial values to function outputs
  high <- 1                                                                    # These define the first row as the highest (and lowest),
  nHigh <- 2                                                                   # and the second row as the second highest

  if (vec_ChiSq[2] > vec_ChiSq[1]){                                            # If the second row's parameter set is worse fitting than
    high <- 2                                                                  # the first row's, consider the second as the current highest
    nHigh <- 1                                                                 # and swap their rankings.
  }

  for (i in 2:nPar){                                                           # For the remainder of rows in the simplex matrix, determine
    if (vec_ChiSq[i] < vec_ChiSq[low]){                                        # how they score relative to the existing highest, second highest,
      low <- i                                                                 # and lowest ranked rows, and adjust the rankings accordingly.
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




PlotRejects <- function(rejected_waves_list1, rejected_waves_list3){
  #################################################################################################################################
  # PlotRejects takes a list of the output list 'reject' from the sep_beats function and consutructs an R plot indicating the number
  # of rejects within a given samples, colour coded to indicate specific reasons for rejection. The function can be used as a tool
  # to assess how conservative the sep_beats function is behaving for a given sample such that appropriate adjustments to
  # empirical thresholds can be altered if necessary. Initially written to work within the ISO study main script.

  # Inputs:
  # rejected_waves_list1 (a list of number of beats rejected for each reason, for a given sample, nested within a list of samples
  # (participants in ISO study)).
  # rejected_waves_list3 (equivalent data structure to rejected_waves_list1 but for a different category of sample (e.g dose level
  # in ISO study) (optional))

  # Output:
  # Automatically plots when called.

  # Notes:
  # The code for this function is inefficient and, while functioning, will need cleaning up for future releases.

  # Reasons for beat rejection include 1. excessively long waveforms (dotted), 2. excessively short waveforms (red), 3. bi-peaked
  # waveforms (blue), 4. waveforms that begin to enter a second systolic upstroke (green), 5. waveforms that drop significantly
  # below baseline (orange), 6. waveforms with significant variability in relation to the average (brown), 7. waveforms that are
  # significantly different in morphology than the average (purple), and total number of rejected beats (black))
  ################################################################################################################################
  if(length(rejected_waves_list3) > 0){                                        # Establish if one or two lists inputted
    len <- 2                                                                   # (for different categories of data e.g dose levels)
  }else{
    len <- 1
  }

  participant_extra_long_waves <- c()                                          # Establish a unique vector of the number of beats rejected
  for(k in 1:len){                                                             # due to being too long in length for each sample (participant)

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

        participant_extra_long_waves[paticip] <- length(participant_extra_long_waves[paticip]) +
          no_of_rejected_waves
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


  participant_extra_short_waves <- c()                                         # Establish a unique vector of the number of beats rejected for being too short
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

        participant_extra_short_waves[paticip] <- length(participant_extra_short_waves[paticip]) +
          no_of_rejected_waves
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

  participant_double_segments <- c()                                           # Establish a unique vector of the number of beats rejected for having two systolic peaks
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

        participant_double_segments[paticip] <- length(participant_double_segments[paticip]) +
          no_of_rejected_waves
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



  participant_systolic_endings <- c()                                          # Establish a unique vector of the number of beats rejected for entering into a second systolic wave
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

        participant_systolic_endings[paticip] <- length(participant_systolic_endings[paticip]) +
          no_of_rejected_waves
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



  participant_drops_below_o <- c()                                             # Establish a unique vector of the number of beats rejected for having values that drop considerably below baseline
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

        participant_drops_below_o[paticip] <- length(participant_drops_below_o[paticip]) +
          no_of_rejected_waves
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


  participant_hrsd_waves <- c()                                                # Establish a unique vector of the number of beats rejected for significant variability in relation to the average
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

        participant_hrsd_waves[paticip] <- length(participant_hrsd_waves[paticip]) +
          no_of_rejected_waves
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


  participant_outlier_waves <- c()                                             # Establish a unique vector of the number of beats rejected for significantly different morphology to the average
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

        participant_outlier_waves[paticip] <- length(participant_outlier_waves[paticip]) +
          no_of_rejected_waves
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

  if(length(participant_extra_long_waves) < length(Participants)){             # Ensure all created vectors are of equal length
    diff <- length(Participants) - length(participant_extra_long_waves)
    participant_extra_long_waves <- c(participant_extra_long_waves,
                                      rep(0, diff))
  }
  if(length(participant_extra_short_waves) < length(Participants)){
    diff <- length(Participants) - length(participant_extra_short_waves)
    participant_extra_short_waves <- c(participant_extra_short_waves,
                                       rep(0, diff))
  }
  if(length(participant_double_segments) < length(Participants)){
    diff <- length(Participants) - length(participant_double_segments)
    participant_double_segments <- c(participant_double_segments,
                                     rep(0, diff))
  }
  if(length(participant_systolic_endings) < length(Participants)){
    diff <- length(Participants) - length(participant_systolic_endings)
    participant_systolic_endings <- c(participant_systolic_endings,
                                      rep(0, diff))
  }
  if(length(participant_drops_below_o) < length(Participants)){
    diff <- length(Participants) - length(participant_drops_below_o)
    participant_drops_below_o <- c(participant_drops_below_o,
                                   rep(0, diff))
  }
  if(length(participant_hrsd_waves) < length(Participants)){
    diff <- length(Participants) - length(participant_hrsd_waves)
    participant_hrsd_waves <- c(participant_hrsd_waves,
                                rep(0, diff))
  }
  if(length(participant_outlier_waves) < length(Participants)){
    diff <- length(Participants) - length(participant_outlier_waves)
    participant_outlier_waves <- c(participant_outlier_waves,
                                   rep(0, diff))
  }

  total_rejected_beats <- participant_extra_long_waves +                       # Calculate total number of rejected beats (irrespective of reason) for each sample
    participant_extra_short_waves + participant_double_segments +
    participant_systolic_endings + participant_drops_below_o +
    participant_hrsd_waves + participant_outlier_waves


  plot(participant_extra_long_waves, t = "l", ylim = c(0, 30),                 # Plot all vectors (x = 1:number of samples, y = 0: number of rejected beats, col = reason for rejection)
       xlim = c(1, 112), ylab = "rejected beats (absolute)",                   # The x and y limits defined here are specified for ISO study data and may need altering.
       xlab = "participants", lty = "dotted", lwd = 1.5)
  lines(participant_extra_short_waves, col = "red",
        lty = "dotted", lwd = 1.5)
  lines(participant_double_segments, col = "blue",
        lty = "dotted", lwd = 1.5)
  lines(participant_systolic_endings, col = "green",
        lty = "dotted", lwd = 1.5)
  lines(participant_drops_below_o, col = "orange",
        lty = "dotted", lwd = 1.5)
  lines(participant_hrsd_waves, col = "brown",
        lty = "dotted", lwd = 1.5)
  lines(participant_outlier_waves, col = "purple",
        lty = "dotted", lwd = 1.5)
  lines(total_rejected_beats, lwd = 1)

}


PlotWavesCarriedForward <- function(waves_carried_forward1, waves_carried_forward3){
  #################################################################################################################################
  # PlotWavesCarriedForward is a function specific to the ISO study data. The number of waves subsetted in the sep_beats function
  # is plotted for each participant (for each dose level) to allow for assessment of which samples are likely to be more or less
  # robust due to the number of beats included.

  # Inputs:
  # waves_carried_forward1 (list of number of beats subsetted for further analysis in the script, at iso dose level)
  # waves_carried_forward3 (list of number of beats subsetted for further analysis in the script, at saline dose level)

  # Output:
  # A histogram is automatically plotted indicating the distribution of frequency of subsetted beats across the dataset
  ################################################################################################################################
  test_vec1 <- waves_carried_forward1
  waves_carried_over <- c()                                                    # Create vector of number of subsetted beats for each
  for(i in 1:length(test_vec1)){                                               # sample at the isoprenaline dose level
    if(!is.null(test_vec1[[i]])){
      waves_carried_over[i] <- test_vec1[[i]]
    }
  }

  if(length(waves_carried_forward3) > 0){                                      # Create vector of number of subsetted beats for each
    test_vec2 <- waves_carried_forward3                                        # sample at the saline dose level
    waves_carried_over2 <- c()
    for(i in 1:length(test_vec2)){
      if(!is.null(test_vec2[[i]])){
        waves_carried_over2[i] <- test_vec2[[i]]
      }
    }
    waves_carried_over <- c(waves_carried_over, waves_carried_over2)           # Combine both inputted lists into a single vector
  }

  hist(waves_carried_over, breaks = 30, xlim = c(0, 200))                      # Create histogram from waves_carried_over vector

}
