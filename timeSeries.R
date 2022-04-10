# TIMESERIES

createTimeSeries <- function(data, label, normCheck = TRUE) {
  ts <- list()
  ts[["data"]] <- data
  ts[["label"]] <- label
  ts[["normed"]] <- FALSE
  ts[["mean"]] <- 0
  ts[["std"]] <- 1
  ts[["normCheck"]] <- normCheck
  return(ts)
}

normalize <- function(tseries, normMean) {
  tseries[["mean"]] <- mean(tseries[["data"]])
  tseries <- calcStddev(tseries)
  thisNorm = !(tseries[["normed"]])
  if (tseries[["normCheck"]] & thisNorm) tseries <- normWork(tseries, normMean)
  return(tseries)
}

# calculates standard deviation of a given time series
calcStddev <- function(tseries) {
  var = 0
  for (i in 1:length(tseries[["data"]])) {
    var = var + tseries[["data"]][i] * tseries[["data"]][i]
  }
  myNorm = 1 / length(tseries[["data"]])
  buf = myNorm * var - tseries[["mean"]] * tseries[["mean"]]
  if (buf != 0) tseries[["std"]] <- sqrt(buf)
  else tseries[["std"]] <- 0
  return(tseries)
}

# normalizes a given time series
normWork <- function(tseries, normMean) {
  if (tseries[["std"]] == 0) ISTD = 1
  else ISTD = 1 / tseries[["std"]]
  if (normMean != 0) {
    data_copy <- vector(length = length(tseries[["data"]]))
    for (i in 1:length(tseries[["data"]])) {
      data_copy[i] = (tseries[["data"]][i] - tseries[["mean"]]) * ISTD
    }
    tseries[["data"]] <- data_copy
    tseries[["mean"]] <- 0
  }
  else if (ISTD != 1) {
    data_copy <- vector(length = length(tseries[["data"]]))
    for (i in 1:length(tseries[["data"]])) {
      data_copy[i] = tseries[["data"]][i] * ISTD
    }
    tseries[["data"]] <- data_copy
  }
  tseries[["normed"]] <- TRUE
  return(tseries)
}

# NONE

# creates independent windows and normalizes them
getDisjointSequences <- function(samples, index, windowSize, normMean) {
  tseries <- samples[["tseries"]][[index]]
  amount = floor(length(tseries[["data"]]) / windowSize)
  subseqences <- list()
  for (i in 1:amount) {
    subseqences_data <- createTimeSeries(tseries[["data"]][((i - 1) * windowSize + 1):(i * windowSize)], 
                                         samples[["labels"]][[index]], TRUE)
    subseqences_data <- normalize(subseqences_data, normMean)
    subseqences[[i]] <- subseqences_data
  }
  return(subseqences)
}

# calculates means and standard deviations of all sliding windows
calcIncrementalMeanStddev <- function(windowLength, series, means, STDs) {
  mySum = 0
  sumSquared = 0
  rWindowLength = 1 / windowLength
  for (i in 1:windowLength) {
    mySum = mySum + series[[i]]
    sumSquared = sumSquared + series[[i]] * series[[i]]
  }
  means <- c(means, mySum * rWindowLength)
  buf = sumSquared * rWindowLength - means[1] * means[1]
  
  if (buf > 0) STDs <- c(STDs, sqrt(buf))
  else STDs <- c(STDs, 0)
  
  # moved by 1
  for (i in 2:(length(series) - windowLength + 1)) {
    mySum = mySum + series[[i + windowLength - 1]] - series[[i - 1]]
    means <- c(means, mySum * rWindowLength)
    sumSquared = sumSquared + series[[i + windowLength - 1]] * series[[i + windowLength - 1]] - series[[i - 1]] * series[[i - 1]]
    buf = sumSquared * rWindowLength - means[i] * means[i]
    if (buf > 0) STDs <- c(STDs, sqrt(buf))
    else STDs <- c(STDs, 0)
  }
  value <- list()
  value[[1]] <- means
  value[[2]] <- STDs
  return(value)
}

# createWord <- function(numbers, maxF, bits) {
#   shortsPerLong = round(60 / bits)
#   to = min(length(numbers), maxF)
#   b = 0
#   shiftOffset = 1
#   for (i in 1:min(to, shortsPerLong)) {
#     shift = 1
#     for (j in 1:bits) {
#       if (altAnd(numbers[i], shift) != 0) b = altOr(b, shiftOffset)
#       shiftOffset = altShiftL(shiftOffset, 1)
#       shift = altShiftL(shift, 1)
#     }
#   }
#   limit = 2147483647
#   total = 2147483647 + 2147483648
#   while (b > limit) {
#     b = b - total - 1
#   }
#   return(b)
# }

# no alternative classifiers implemented
# compareTo <- function(score, bestScore) {
#   value = 1
#   if (score[1] > bestScore[1]) value = -1
#   else if (score[1] == bestScore[1] & score[4] > bestScore[4]) value = -1
#   return(value)
# }

# UNLIMITED BITWISE CALCULATIONS
# these operations are not limited to 32-bit numbers
altAnd <- function(v1, v2) {
  i = 0
  ans = 0
  while (v1 > 0 & v2 > 0) {
    if (v1 %% 2 == v2 %% 2) ans = ans + 2 ^ i
    v1 = v1 %/% 2
    v2 = v2 %/% 2
    i = i + 1
  }
  return(ans)
}

altOr <- function(v1, v2) {
  i = 0
  ans = 0
  while (v1 > 0 | v2 > 0) {
    if (v1 %% 2 == 1 | v2 %% 2 == 1) ans = ans + 2 ^ i
    v1 = v1 %/% 2
    v2 = v2 %/% 2
    i = i + 1
  }
  return(ans)
}

altShiftL <- function(v1, v2) {
  for (i in 1:v2) {
    v1 = v1 * 2
  }
  return(v1)
}

altShiftR <- function(v1, v2) {
  for (i in 1:v2) {
    v1 = v1 %/% 2
  }
  return(v1)
}

# calculates how many bits are necessary to represent a given number as binary
# combinations - limited to 64-bit numbers
calcBits <- function(number) {
  ans = 0
  if (number >= 4294967296) {
    number = altShiftR(number, 32)
    ans = 32
  }
  if (number >= 65536) {
    number = altShiftR(number, 16)
    ans = ans + 16
  }
  if (number >= 256) {
    number = altShiftR(number, 8)
    ans = ans + 8
  }
  if (number >= 16) {
    number = altShiftR(number, 4)
    ans = ans + 4
  }
  if (number >= 4) {
    number = altShiftR(number, 2)
    ans = ans + 2
  }
  ans = ans + altShiftR(number, 1)
  return(ans)
}

# loads and normalizes chosen training and testing datasets
uvLoad <- function(datasetName) {
  name_train = paste("./datasets/univariate/", datasetName, "/", datasetName, "_TRAIN", sep = "")
  name_test = paste("./datasets/univariate/", datasetName, "/", datasetName, "_TEST", sep = "")
  train_raw <- read.csv2(name_train, header = FALSE, sep = " ", dec = ".")
  test_raw <- read.csv2(name_test, header = FALSE, sep = " ", dec = ".")
  
  train <- list()
  test <- list()
  train[["samples"]] <- nrow(train_raw) 
  train[["size"]] <- ncol(train_raw) - 1                         # no 1st column
  train[["labels"]] <- vector()
  train[["tseries"]] <- list()
  test[["samples"]] <- nrow(test_raw)
  test[["size"]] <- ncol(test_raw) - 1
  test[["labels"]] <- vector()
  test[["tseries"]] <- list()
  
  for (i in 1:train[["samples"]]) {
    label = train_raw[i, 1]
    train[["labels"]] <- c(train[["labels"]], label)
    series <- unlist(train_raw[i, -1])
    train[["tseries"]][[i]] <- createTimeSeries(series, label)
    train[["tseries"]][[i]] <- normalize(train[["tseries"]][[i]], 1)
    # train[["tseries"]][[i]][["normed"]] = TRUE
  }
  message <- paste("Loaded training data from set ", datasetName, ": ", train[["samples"]], 
                   " train samples of length ",  train[["size"]], ".", sep = "")
  print(message)
  
  
  for (i in 1:test[["samples"]]) {
    label = test_raw[i, 1]
    test[["labels"]] <- c(test[["labels"]], label)
    series <- unlist(test_raw[i, -1])
    test[["tseries"]][[i]] <- createTimeSeries(series, label)
    test[["tseries"]][[i]] <- normalize(test[["tseries"]][[i]], 1)
    # test[["tseries"]][[i]][["normed"]] = TRUE
  }
  message <- paste("Loaded testing data from set ", datasetName, ": ", test[["samples"]], 
                   " test samples of length ",  test[["size"]], ".", sep = "")
  print(message)
  
  value <- list()
  value[[1]] <- train
  value[[2]] <- test
  return(value)
}