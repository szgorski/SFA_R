createMFT <- function(windowSize, normMean, lowerBounding) {
  mft <- list()
  mft[["windowSize"]] <- windowSize
  if (normMean == TRUE) mft[["startOffset"]] <- 2
  else mft[["startOffset"]] <- 0
  if (lowerBounding == TRUE) mft[["norm"]] <- 1 / sqrt(windowSize)
  else mft[["norm"]] <- 1
  return(mft)
}

# MFT

# applies FFT on an independent window, choosing first length/2 coefficients,
# dropping the first one - last two values are set to 0
transform <- function(mft, series, wordLength) {
  FFTSeries <- fft(unlist(series))
  data_new <- vector()
  windowSize = length(series)
  for (i in 1:ceiling(windowSize / 2)) {
    data_new <- c(data_new, Re(FFTSeries[i]), Im(FFTSeries[i]))
  }
  data_new[2] = 0
  data_new <- data_new[1:windowSize]
  dataLength = min(windowSize - mft[["startOffset"]], wordLength)
  data_copy = data_new[(mft[["startOffset"]] + 1):(mft[["startOffset"]] + dataLength)]
  while (length(data_copy) != wordLength) {
    data_copy <- c(data_copy, 0)
  }
  
  sign = 1
  for (i in 1:length(data_copy)) {
    data_copy[i] = data_copy[i] * mft[["norm"]] * sign
    sign = -sign
  }
  return(data_copy)
}

# applies truncated FT to sliding windows and normalizes them
transformWindowing <- function(mft, fullSeries, wordLength) {
  series <- fullSeries[["data"]]             # tseries from uv_load, row of data
  newWordLength = min(mft[["windowSize"]], wordLength + mft[["startOffset"]])
  newWordLength = newWordLength + newWordLength %% 2
  phis <- rep(0, newWordLength)
  for (i in seq.int(1, newWordLength, by = 2)) {
    half_i = - (i - 1) / 2
    phis[i] = cospi(2 * half_i / mft[["windowSize"]])
    phis[i+1] = - sinpi(2 * half_i / mft[["windowSize"]])
  }
  
  final = max(1, length(series) - mft[["windowSize"]] + 1)
  mft[["means"]] <- vector()
  mft[["STDs"]] <- vector()
  catch <- list()
  catch <- calcIncrementalMeanStddev(mft[["windowSize"]], series, mft[["means"]], 
                                     mft[["STDs"]])
  mft[["means"]] <- catch[[1]]
  mft[["STDs"]] <- catch[[2]]
  
  transformed <- list()
  mftData <- vector()
  for (t in 1:final) {
    if (t == 1) {
      fftData <- fft(unlist(series[1:mft[["windowSize"]]]))
      mftData <- rep(0, newWordLength)
      i = 1
      for (j in 1:min(mft[["windowSize"]], newWordLength)) {
        if (j %% 2 == 1) mftData[j] = Re(fftData[i])
        else {
          mftData[j] = Im(fftData[i])
          i = i + 1
        }
      }
      mftData[2] = 0
    }
    else {
      k = 1
      while (k <= newWordLength) {
        real_raw = mftData[k] + series[t + mft[["windowSize"]] - 1] - series[t - 1]
        imag_raw = mftData[k + 1]
        real = real_raw * phis[k] - imag_raw * phis[k + 1]
        imag = real_raw * phis[k + 1] + imag_raw * phis[k]
        mftData[k] = real
        mftData[k + 1] = imag
        k = k + 2
      }
    }
    
    copy <- rep(0, wordLength)
    copy_value <- mftData[(mft[["startOffset"]] + 1):min(length(mftData), (mft[["startOffset"]] + wordLength))]
    copy[1:length(copy_value)] <- copy_value
    copy <- createTimeSeries(copy, fullSeries[["label"]], fullSeries[["normCheck"]])
    copy <- normalizeFT(mft, copy, mft[["STDs"]][t])
    transformed[[t]] <- copy
  }
  value <- list()
  value[[1]] <- mft
  value[[2]] <- transformed
  return(value)
}

# normalizes transformed sliding windows
normalizeFT <- function(mft, copy, std) {
  if (copy[["normCheck"]] == TRUE & std > 0) normalizingFactor = 1 / std
  else normalizingFactor = 1
  normalizingFactor = normalizingFactor * mft[["norm"]]
  
  sign = 1
  for (i in 1:length(copy[["data"]])) {
    copy[["data"]][i] = copy[["data"]][i] * sign * normalizingFactor
    sign = -sign
  }
  return(copy[["data"]])
}
