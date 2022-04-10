createSFA <- function(histogramType, SUP = FALSE, lowerBounding = TRUE) {
  sfa <- list()
  sfa[["initialized"]] <- FALSE
  sfa[["histogramType"]] <- histogramType
  sfa[["SUP"]] <- SUP
  sfa[["lowerBounding"]] <- lowerBounding
  return(sfa)
}

# SFA

initializeSFA <- function(sfa, wordLength, symbolSet, normMean) {
  sfa[["initialized"]] <- TRUE
  sfa[["wordLength"]] <- wordLength
  sfa[["maxWordLength"]] <- wordLength
  sfa[["symbolSet"]] <- symbolSet
  sfa[["normMean"]] <- normMean
  sfa[["alphabetSize"]] <- symbolSet
  sfa[["transformation"]] <- NULL
  sfa[["orderLine"]] <- NULL
  sfa[["bins"]] <- data.frame(matrix(Inf, wordLength, sfa[["alphabetSize"]]))
  sfa[["bins"]][, 1] = -Inf
  return(sfa)
}

# fits new SFAs - creates windows and transforms them into letters
fitWindowing <- function(sfa, samples, windowSize, wordLength, symbolSet, normMean, lowerBounding) {
  sfa[["transformation"]] <- createMFT(windowSize, normMean, lowerBounding)
  sa <- list()
  sa[["tseries"]] <- list()
  sa[["labels"]] <- list()
  index = 0
  
  for (i in 1:samples[["samples"]]) {
    newList <- getDisjointSequences(samples, i, windowSize, normMean)
    for (j in 1:length(newList)) {
      index = index + 1
      sa[["tseries"]][[index]] <- newList[[j]]
      sa[["labels"]][[index]] <- newList[[j]][["label"]]
    }
  }
  sa[["samples"]] <- index
  sa[["size"]] <- windowSize
  catch <- list()
  if (sfa[["SUP"]] == TRUE) 
    catch <- fitTransformSupervised(sfa, sa, wordLength, symbolSet, normMean)
  else catch <- fitTransform(sfa, sa, wordLength, symbolSet, normMean)
  return(catch[[1]])
}

# fitTransform <- function(sfa, samples, wordLength, symbols, normMean) {
#   catch <- list()
#   catch <- fitTransformDouble(sfa, samples, wordLength, symbols, normMean)
#   transformed <- multiTransform(catch[[1]], samples, catch[[2]])
#   return(catch[[1]], transformed)
# }

# applies transformation function (singleTransform) to each independent window
multiTransform <- function(sfa, samples, approximate) {
  transformed <- list()
  for (i in 1:samples[["samples"]]) {
    transformed[[i]] <- singleTransform(sfa, samples[["tseries"]][[i]][["data"]], approximate[[i]])
  }
  return(transformed) # list of vector of numbers
}

# transforms a given independent window and assigns letters to it
singleTransform <- function(sfa, series, singleApproximate) {
  if (is.null(singleApproximate))
    singleApproximate <- transform(sfa[["transformation"]], series, sfa[["maxWordLength"]])
  if (sfa[["SUP"]] == TRUE) 
    value <- quantizationSupervised(sfa, singleApproximate)
  else value <- quantization(sfa, singleApproximate)
  return(value) # vector of numbers
}

# creates letters from a time series' sliding windows
SFAWindowing <- function(sfa, series) {
  catch <- list()
  catch <- transformWindowing(sfa[["transformation"]], series, sfa[["maxWordLength"]])
  sfa[["transformation"]] <- catch[[1]]
  transformedMFT <- catch[[2]]
  words <- vector("list", length = length(transformedMFT))
  for (i in 1:length(transformedMFT)) {
    if (sfa[["SUP"]] == TRUE) 
      words[[i]] <- quantizationSupervised(sfa, transformedMFT[[i]])
    else words[[i]] <- quantization(transformedMFT[[i]])
  }
  value <- list()
  value[[1]] <- sfa
  value[[2]] <- words
  return(value)
}

# creates words from a time series' sliding windows
SFAWindowingInt <- function(sfa, series, wordLength) {
  catch <- list()
  catch <- SFAWindowing(sfa, series)
  sfa <- catch[[1]]
  words <- catch[[2]]
  intWords <- vector(length = length(words))
  for (i in 1:length(words)) {
    bits = calcBits(sfa[["alphabetSize"]])
    intWords[i] = createUnlimitedWord(sfa, words[[i]], wordLength, bits)
  }
  value <- list()
  value[[1]] <- sfa
  value[[2]] <- intWords
  return (value)
}

# creates bins from independent windows
fitTransformDouble <- function(sfa, samples, wordLength, symbolSet, normMean) {
  if (sfa[["initialized"]] == FALSE) {
    sfa <- initializeSFA(sfa, wordLength, symbolSet, normMean)
    if (is.null(sfa[["transformation"]])) 
      sfa[["transformation"]] <- createMFT(length(samples[["tseries"]][[1]][["data"]]),
                                           normMean, sfa[["lowerBounding"]])
  }
  catch <- list()
  catch <- fillOrderLine(sfa, samples, wordLength)
  sfa <- catch[[1]]
  transformedSamples <- catch[[2]]
  if (sfa[["histogramType"]] == "EQUI_DEPTH") 
    sfa <- divideEquiDepthHistogram(sfa)
  else if (sfa[["histogramType"]] == "EQUI_FREQUENCY") 
    sfa <- divideEquiWidthHistogram(sfa)
  else if (sfa[["histogramType"]] == "INFORMATION_GAIN") 
    sfa <- divideHistogramInformationGain(sfa)
  value <- list()
  value[[1]] <- sfa
  value[[2]] <- transformedSamples
  return(value)
}

# transforms independent windows with FFT and divides them into groups 
# with respect to classes
fillOrderLine <- function(sfa, samples, wordLength) {
  sfa[["orderLine"]] <- vector("list", wordLength)
  for (i in 1:wordLength) {
    sfa[["orderLine"]][[i]] <- vector("list", samples[["samples"]])
  }
  transformedSamples <- list()
  
  for(i in 1:samples[["samples"]]) {
    transformedSamplesUnit <- transform(sfa[["transformation"]], samples[["tseries"]][[i]][["data"]], 
                                        wordLength)
    transformedSamples[[i]] <- transformedSamplesUnit
    for (j in 1:length(transformedSamplesUnit)) {
      obj <- list()
      obj[["value"]] <- round(transformedSamplesUnit[[j]], 2)
      obj[["label"]] <- samples[["labels"]][[i]]
      sfa[["orderLine"]][[j]][[i]] <- obj
    }
  }
  
  for (j in 1:length(sfa[["orderLine"]])) {
    delList <- sfa[["orderLine"]][[j]]
    newList <- list()
    index = 0
    while (length(delList) != 0) {
      index = index + 1
      currentMinValue = Inf
      currentMinLocation = -1
      currentLabel = -Inf
      for (i in 1:length(delList)) {
        if (delList[[i]][[1]] < currentMinValue | (delList[[i]][[1]] == currentMinValue & delList[[i]][[2]] < currentLabel)) {
          currentMinValue <- delList[[i]][[1]]
          currentLabel <- delList[[i]][[2]]
          currentMinLocation <- i
        }
      }
      newList[[index]] <- delList[[currentMinLocation]]
      delList <- delList[-currentMinLocation]
    }
    sfa[["orderLine"]][[j]] <- newList
  }
  value <- list()
  value[[1]] <- sfa
  value[[2]] <- transformedSamples
  return(value)
}

# quantization <- function(sfa, singleApproximate) {
#   i = 1
#   word <- rep(0, length(singleApproximate))
#   for (v in singleApproximate) {
#     c = 1
#     for (n in 1:ncol(sfa[["bins"]])) {
#       if (v < sfa[["bins"]][i, c]) break
#       else c = c + 1
#     }
#     word[i] = c - 1
#     i = i + 1
#   }
#   return(word)
# }

# assigns letters to values of independent windows for each F coefficient 
# (using bins)
quantizationSupervised <- function(sfa, singleApproximate) {
  signal <- rep(0, min(length(singleApproximate), length(sfa[["bestValues"]])))
  for (a in 1:length(signal)) {
    i = sfa[["bestValues"]][a]
    b = 1
    for (x in 1:ncol(sfa[["bins"]])) {
      if (singleApproximate[i] < sfa[["bins"]][i, b]) break
      else b = b + 1
    }
    signal[a] = b - 1
  }
  return(signal)
}

# divideEquiDepthHistogram <- function(sfa) {
#   for (i in 1:nrow(sfa[["bins"]])) {
#     tryCatch(
#       expr = {
#         depth = length(sfa[["orderLine"]][[i]]) / sfa[["alphabetSize"]]
#       },
#       error = function(e){
#         depth = 0
#       }
#     )
#     pos = 0
#     count = 0
#     tryCatch(
#       expr = {
#         for (j in 1:length(sfa[["orderLine"]][i])) {
#           count = count + 1
#           if (count > ceiling(depth * pos) & (pos == 1 | sfa[["bins"]][i, pos - 1] != sfa[["orderLine"]][[i]][[j]][[1]])) {
#             sfa[["bins"]][i, pos] = round(sfa[["orderLine"]][[i]][[j]][[1]], 2)
#             pos = pos + 1
#           }
#         }
#       }
#     )
#   }
#   sfa[["bins"]][, 1] = -Inf
#   return(sfa)
# }

# divideEquiWidthHistogram <- function(sfa) {
#   i = 0
#   for (elem in sfa[["orderLine"]]) {
#     if (length(elem) != 0) {
#       intervalWidth = (elem[[1]][[1]] - elem[[length(elem)]][[1]]) / sfa[["alphabetSize"]]
#       for (c in 1:(sfa[["alphabetSize"]] - 1)) {
#         sfa[["bins"]][i, c] = intervalWidth * (c + 1) + elem[[1]][[1]]
#       }
#     }
#     i = i + 1
#   }
#   sfa[["bins"]][, 1] = -Inf
#   return(sfa)
# }

# puts values of split points to bins
divideHistogramInformationGain <- function(sfa) {
  for (i in 1:length(sfa[["orderLine"]])) {
    sfa[["splitPoints"]] <- vector()
    sfa <- findBestSplit(sfa, sfa[["orderLine"]][[i]], 1,
                         length(sfa[["orderLine"]][[i]]), sfa[["alphabetSize"]])
    sfa[["splitPoints"]] <- sort(sfa[["splitPoints"]])
    for (j in 1:length(sfa[["splitPoints"]])) {
      sfa[["bins"]][i, j + 1] = sfa[["orderLine"]][[i]][[sfa[["splitPoints"]][j] + 1]][[1]]
    }
  }
  return(sfa)
}

# finds points giving highest information gain, divides split domain 
# to subdomians
findBestSplit <- function(sfa, element, startIndex, endIndex, remainingSymbols) {
  bestGain = -1
  bestPos = -1
  total = endIndex - startIndex + 1
  
  cOut <- list()
  for (pos in startIndex:endIndex) {
    label = element[[pos]][[2]]
    if (as.character(label) %in% names(cOut)) 
      cOut[[as.character(label)]] = cOut[[as.character(label)]] + 1
    else cOut[[as.character(label)]] <- 1
  }
  sfa[["cOut"]] <- cOut
  sfa[["cIn"]] <- list()
  
  classEntropy = entropy(sfa, sfa[["cOut"]], total)
  i = startIndex
  lastLabel = element[[i]][[2]]
  sfa <- moveElement(sfa, lastLabel)
  for (split in (startIndex + 1):endIndex) {
    i = i + 1
    label = element[[i]][[2]]
    sfa <- moveElement(sfa, label)
    if (label != lastLabel) {
      gain <- calculateInformationGain(sfa, classEntropy, i, total)
      if (gain > bestGain) {
        bestGain = gain
        bestPos = split
      }
    }
    lastLabel = label
  }
  
  if (bestPos > -1) {
    sfa[["splitPoints"]] <- c(sfa[["splitPoints"]], bestPos)
    remainingSymbols = remainingSymbols / 2
    if (remainingSymbols > 1) {
      if (bestPos - startIndex > 2 & endIndex - bestPos >= 2) {
        sfa <- findBestSplit(sfa, element, startIndex, bestPos - 1, remainingSymbols)
        sfa <- findBestSplit(sfa, element, bestPos, endIndex, remainingSymbols)
      }
      else if (endIndex - bestPos >= 4) {
        sfa <- findBestSplit(sfa, element, bestPos, (endIndex - bestPos + 1) %/% 2 - 1, 
                             remainingSymbols)
        sfa <- findBestSplit(sfa, element, (endIndex - bestPos + 1) %/% 2, endIndex, 
                             remainingSymbols)
      }
      else if (bestPos - startIndex > 4) {
        sfa <- findBestSplit(sfa, element, startIndex, (bestPos - startIndex) %/% 2 - 1, 
                             remainingSymbols)
        sfa <- findBestSplit(sfa, element, (bestPos - startIndex) %/% 2, endIndex, 
                             remainingSymbols)
      }
    }
  }
  return(sfa)
}

# uses cIn and cOut to trace the usage of independent windows by labels
moveElement <- function(sfa, label) {
  if (as.character(label) %in% names(sfa[["cIn"]])) 
    sfa[["cIn"]][[as.character(label)]] = sfa[["cIn"]][[as.character(label)]] + 1
  else sfa[["cIn"]][[as.character(label)]] <- 1
  if (as.character(label) %in% names(sfa[["cOut"]])) 
    sfa[["cOut"]][[as.character(label)]] = sfa[["cOut"]][[as.character(label)]] - 1
  else sfa[["cOut"]][[as.character(label)]] <- -1
  return(sfa)
}

entropy <- function(sfa, freq, total) {
  e = 0
  if (total != 0) {
    logarithm = 1 / log(2)
    for (k in names(freq)) {
      p = freq[[k]] / total
      if (p > 0) e = e - p * log(p) * logarithm
    }
  }
  else e = Inf
  return(e)
}

calculateInformationGain <- function(sfa, classEntropy, totalCIn, total) {
  totalCOut = total - totalCIn
  value = classEntropy
  eIn = entropy(sfa, sfa[["cIn"]], totalCIn)
  eOut = entropy(sfa, sfa[["cOut"]], totalCOut)
  if (eIn != Inf) value = value - totalCIn / total * eIn
  if (eOut != Inf) value = value - totalCOut / total * eOut
  return(value)
}

# takes all letters and joins them into a number representing word
createUnlimitedWord <- function(sfa, numbers, maxF, bits) {
  shortsPerLong = round(60 / bits)
  to = min(length(numbers), maxF)
  b = 0
  shiftOffset = 1
  for (i in 1:min(to, shortsPerLong)) {
    shift = 1
    for (j in 1:bits) {
      if (altAnd(numbers[i] - 1, shift) != 0) b = altOr(b, shiftOffset)
      shiftOffset = altShiftL(shiftOffset, 1)
      shift = altShiftL(shift, 1)
    }
  }
  return(b)
}

# assigns letters to independent windows
fitTransformSupervised <- function(sfa, samples, wordLength, symbolSet, normMean) {
  len = length(samples[["tseries"]][[1]][["data"]])
  catch <- list()
  catch <- fitTransformDouble(sfa, samples, len, symbolSet, normMean)
  sfa <- catch[[1]]
  transformedSignal <- catch[[2]]
  best <- calcBestCoefficients(sfa, samples, transformedSignal)
  sfa[["bestValues"]] <- rep(0, min(length(best), wordLength))
  sfa[["maxWordLength"]] = 0
  for (i in 1:length(sfa[["bestValues"]])) {
    sfa[["bestValues"]][[i]] <- best[[i]][[1]]
    sfa[["maxWordLength"]] <- max(best[[i]][[1]], sfa[["maxWordLength"]])
  }
  sfa[["maxWordLength"]] <- sfa[["maxWordLength"]] + sfa[["maxWordLength"]] %% 2
  value <- list()
  value[[1]] <- sfa
  value[[2]] <- multiTransform(sfa, samples, transformedSignal)
  return(value)
}

# divides normalized independent windows into classes, performs ANOVA-F test 
# and sorts calculated coefficients by their information gain
calcBestCoefficients <- function(sfa, samples, transformedSignal) {
  classes <- list()
  for (i in 1:samples[["samples"]]) {
    if (as.character(samples[["labels"]][[i]]) %in% names(classes)) {
      classes[[as.character(samples[["labels"]][[i]])]] <- append(classes[[as.character(samples[["labels"]][[i]])]], 
                                                                  transformedSignal[i])
    }
    else {
      classes[[as.character(samples[["labels"]][[i]])]] <- list()
      classes[[as.character(samples[["labels"]][[i]])]] <- append(classes[[as.character(samples[["labels"]][[i]])]], 
                                                                  transformedSignal[i])
    }
  }
  nSamples = length(transformedSignal)
  nClasses = length(classes)
  len = length(transformedSignal[[1]])
  f = getFoneway(sfa, len, classes, nSamples, nClasses)
  fSorted = sort(f, decreasing = TRUE)
  best <- list()
  i = 0
  listIndex = 0
  for (value in fSorted) {
    if (value == Inf) {
      index = match(value, f) + i
      i = i + 1
    }
    else index = match(value, f)
    listIndex = listIndex + 1
    best[[listIndex]] <- c(index, value)
  }
  return(best)
}

# calculates information gain for ANOVA-F coefficients
getFoneway <- function(sfa, len, classes, nSamples, nClasses) {
  ss <- rep(0, len)
  sumsArgs <- list()
  keysOfClasses <- names(classes)
  
  for (key in keysOfClasses) {
    allTs <- classes[[key]]
    sums <- rep(0, len)
    for (ts in allTs) {
      for (i in 1:length(ts)) {
        ss[i] = ss[i] + ts[i] * ts[i]
        sums[i] = sums[i] + ts[i]
      }
    }
    sumsArgs[[key]] <- sums
  }
  squaredData <- rep(0, len)
  squaredArgs <- list()
  for (key in keysOfClasses) {
    sums <- sumsArgs[[key]]
    for (i in 1:length(sums)) {
      squaredData[i] = squaredData[i] + sums[i]
    }
    squares <- rep(0, length(sums))
    for (i in 1:length(sums)) {
      squares[i] = squares[i] + sums[i] * sums[i]
    }
    squaredArgs[[key]] <- squares
  }
  for (i in 1:len) {
    squaredData[i] = squaredData[i] * squaredData[i]
  }
  ssTotal <- rep(0, len)
  for (i in 1:len) {
    ssTotal[i] = ss[i] - squaredData[i] / nSamples
  }
  ssBetween <- rep(0, len)
  ssWithin <- rep(0, len)
  for (key in keysOfClasses) {
    sums <- squaredArgs[[key]]
    nSamplesPerClass = length(classes[[key]])
    for (i in 1:length(sums)) {
      ssBetween[i] = ssBetween[i] + sums[i] / nSamplesPerClass
    }
  }
  for (i in 1:len) {
    ssBetween[i] = ssBetween[i] - squaredData[i] / nSamples
  }
  
  degOfFreedomBetween = nClasses - 1
  degOfFreedomWithin = nSamples - nClasses
  meanSqBetween <- rep(0, len)                                 # between classes
  meanSqWithin <- rep(0, len)                                   # within classes
  f <- rep(0, len)                                                     # f-ratio
  for (i in 1:len) {
    ssWithin[i] = ssTotal[i] - ssBetween[i]
    meanSqBetween[i] = ssBetween[i] / degOfFreedomBetween
    meanSqWithin[i] = ssWithin[i] / degOfFreedomWithin
    if (meanSqWithin[i] != 0) f[i] = meanSqBetween[i] / meanSqWithin[i]
    else f[i] = Inf
  }
  return(f)
}
