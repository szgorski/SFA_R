maxWindowLength = 20

createWEASEL <- function(maxF, maxS, windowLength, normMean) {
  weasel <- list()
  weasel[["maxF"]] <- maxF
  weasel[["symbolSet"]] <- maxS
  weasel[["windowLengths"]] <- windowLength
  weasel[["normMean"]] <- normMean
  weasel[["signature"]] <- vector("list", length(weasel[["windowLengths"]]))
  weasel[["dict"]] <- createDictionary()
  return(weasel)
}

createDictionary <- function() {
  dictionary <- list()
  dictionary[["dict"]] <- list()
  dictionary[["dictChi"]] <- list()
  return(dictionary)
}

createBagOfBigrams <- function(label) {
  bag <- list()
  bag[["bob"]] <- list()
  bag[["label"]] <- label
  return(bag)
}

# WEASEL

# applies word creation function (createWords) for all window lengths
createFitting <- function(weasel, samples) {
  weasel[["words"]] <- vector("list", length = length(weasel[["windowLengths"]]))
  for (i in 1:length(weasel[["windowLengths"]])) {
    message <- paste("Fitting for windowLength = ", weasel[["windowLengths"]][i],
                     "...", sep = "")
    print(message)
    weasel <- createWords(weasel, samples, i, NULL)            # no progress bar
  }
  value <- list()
  value[[1]] <- weasel
  value[[2]] <- weasel[["words"]]
  return (value)
}

# creates words of a given length
createWords <- function(weasel, samples, index, bar) {
  if (is.null(weasel[["signature"]][[index]])) {
    weasel[["signature"]][[index]] <- createSFA("INFORMATION_GAIN", TRUE, FALSE)
    weasel[["signature"]][[index]] <- fitWindowing(weasel[["signature"]][[index]], 
        samples, weasel[["windowLengths"]][index], weasel[["maxF"]], 
        weasel[["symbolSet"]], weasel[["normMean"]], FALSE)
  }
  wordsList <- list()
  catch <- list()
  for (i in 1:samples[["samples"]]) {
    catch <- SFAWindowingInt(weasel[["signature"]][[index]], samples[["tseries"]][[i]], 
                             weasel[["maxF"]])
    weasel[["signature"]][[index]] <- catch[[1]]
    wordsList[[i]] <- catch[[2]]
  }
  weasel[["words"]][[index]] <- wordsList
  if (!is.null(bar)) {                                         # no progress bar
    setTxtProgressBar(bar, index)
  }
  return(weasel)
}

# creates bigrams, puts unigrams and bigrams into boxes of patterns 
# (for each sample and window length)
createBOP <- function(weasel, words, samples, f) {
  bagOfPatterns <- list()
  for (i in 1:samples[["samples"]]) {
    bagOfPatterns[[i]] <- createBagOfBigrams(samples[["labels"]][[i]])
  }
  usedBits = calcBits(weasel[["symbolSet"]])
  mask = altShiftL(1, usedBits * f) - 1
  highestBit = calcBits(maxWindowLength) + 1
  
  catch <- list()
  for (i in 1:samples[["samples"]]) {
    message <- paste("Creating BOP for sample no. ", i, "...")
    print(message)
    
    for (w in 1:length(weasel[["windowLengths"]])) {
      for (myOffset in 1:length(words[[w]][[i]])) {
        catch <- getWord(weasel[["dict"]], altOr(altShiftL(altAnd(words[[w]][[i]][myOffset], 
                                                                  mask), highestBit), w - 1))
        weasel[["dict"]] <- catch[[1]]
        word <- catch[[2]]
        
        if (as.character(word) %in% names(bagOfPatterns[[i]][["bob"]]))
          bagOfPatterns[[i]][["bob"]][[as.character(word)]] <- bagOfPatterns[[i]][["bob"]][[as.character(word)]] + 1
        else bagOfPatterns[[i]][["bob"]][[as.character(word)]] <- 1
        
        if (myOffset - weasel[["windowLengths"]][w] > 0) {
          catch <- getWord(weasel[["dict"]], altOr(altShiftL(altAnd(words[[w]][[i]][myOffset - weasel[["windowLengths"]][w]], 
                                                                    mask), highestBit), w - 1))
          weasel[["dict"]] <- catch[[1]]
          prevWord <- catch[[2]]
          
          catch <- getWord(weasel[["dict"]], altShiftL(altOr(altShiftL(prevWord, 32), word), highestBit))
          weasel[["dict"]] <- catch[[1]]
          newWord <- catch[[2]]
          
          if (as.character(newWord) %in% names(bagOfPatterns[[i]][["bob"]])) 
            bagOfPatterns[[i]][["bob"]][[as.character(newWord)]] <- bagOfPatterns[[i]][["bob"]][[as.character(newWord)]] + 1
          else bagOfPatterns[[i]][["bob"]][[as.character(newWord)]] <- 1
        }
      }
    }
  }
  value <- list()
  value[[1]] <- weasel
  value[[2]] <- bagOfPatterns
  return(value)
}

# calculates chi^2 for all patterns in all classes and removes elements below 
# the chi limit
filterChiSquared <- function(weasel, bob, chiLimit) {
  classFrequencies <- list()
  featureCount <- list()
  classProb <- list()
  observed <- list()
  chiSquare <- list()
  print("Filtering Chi^2...")
  print("    Step 1 of 5")
  for (i in bob) {
    label = i[["label"]]
    if (as.character(label) %in% names(classFrequencies))
      classFrequencies[[as.character(label)]] <- classFrequencies[[as.character(label)]] + 1
    else classFrequencies[[as.character(label)]] <- 1
  }
  print("    Step 2 of 5")
  counter = 0
  for (bop in bob) {                   # counts number of samples with this word
    counter = counter + 1
    message <- paste("        Allocating BOB no. ", counter, " of ", length(bob), "...")
    print(message)
    
    label = bop[["label"]]
    bagDict = bop[["bob"]]
    for (key in names(bagDict)) {
      if (bagDict[[key]] > 0) {
         if (key %in% names(featureCount)) 
           featureCount[[key]] <- featureCount[[key]] + 1
         else featureCount[[key]] <- 1
         
         key_new = altOr(altShiftL(label, 32), as.numeric(key))
         if (as.character(key_new) %in% names(observed)) 
           observed[[as.character(key_new)]] <- observed[[as.character(key_new)]] + 1
         else observed[[as.character(key_new)]] <- 1
      }
    }
  }
  print("    Step 3 of 5")
  for (lis in bob) {                                  # counts samples per class
    label = lis[["label"]]
    if (as.character(label) %in% names(classProb))
      classProb[[as.character(label)]] <- classProb[[as.character(label)]] + 1
    else classProb[[as.character(label)]] <- 1
  }
  print("    Step 4 of 5")
  counter = 0
  for (prob in names(classProb)) {
    counter = counter + 1
    message <- paste("        Calculating Chi^2 for class no. ", counter, " of ",
                     length(classProb), "...")
    print(message)
    
    classProb[[prob]] = classProb[[prob]] / length(bob)
    for (feature in names(featureCount)) {
      key = altOr(altShiftL(as.numeric(prob), 32), as.numeric(feature))
      expected = classProb[[prob]] * featureCount[[feature]]
      chi = getKey(observed, key) - expected
      chi_new = chi * chi / expected
      if (chi_new >= chiLimit & chi_new > getKey(chiSquare, as.numeric(feature))) 
        chiSquare[[feature]] <- chi_new
    }
  }
  print("    Step 5 of 5")
  for (i in 1:length(bob)) {                         # best elements above limit
    for (key in names(bob[[i]][["bob"]])) {
      if (getKey(chiSquare, as.numeric(key)) < chiLimit) 
        bob[[i]][["bob"]][[key]] <- 0
    }
  }
  catch <- list()
  catch <- remap(weasel[["dict"]], bob)
  weasel[["dict"]] <- catch[[1]]
  bob <- catch[[2]]
  value <- list()
  value[[1]] <- weasel
  value[[2]] <- catch[[2]]
  return(value)
}

# DICTIONARY

resetDictionary <- function(dictionary) {
  dictionary[["dict"]] <- list()
  dictionary[["dictChi"]] <- list()
  return(dictionary)
}

# returns position of a given word in a dictionary and counts its occurrences, 
# if one hasn't occurred yet, puts it at the end
getWord <- function(dictionary, word) {
  if (as.character(word) %in% names(dictionary[["dict"]])) 
    word_new <- dictionary[["dict"]][[as.character(word)]]
  else {
    word_new <- length(names(dictionary[["dict"]])) + 1
    dictionary[["dict"]][[as.character(word)]] <- word_new
  }
  value <- list()
  value[[1]] <- dictionary
  value[[2]] <- word_new
  return(value)
}

# returns position of a given feature in a dictionary and counts its
# occurrences, if one hasn't occurred yet, puts it at the end
getWordChi <- function(dictionary, word) {
  if (as.character(word) %in% names(dictionary[["dictChi"]])) 
    word_new <- dictionary[["dictChi"]][[as.character(word)]]
  else {
    word_new <- length(names(dictionary[["dictChi"]])) + 1
    dictionary[["dictChi"]][[as.character(word)]] <- word_new
  }
  value <- list()
  value[[1]] <- dictionary
  value[[2]] <- word_new
  return(value)
}

# returns the size of a dictionary
getSize <- function(dictionary) {
  if (length(dictionary[["dictChi"]]) != 0) 
    value <- length(dictionary[["dictChi"]]) + 1
  else value <- length(dictionary[["dict"]])
  return(value)
}

# removes entries filtered by chi limit
remap <- function(dictionary, bagOfPatterns) {
  print("Remapping...")
  catch <- list()
  for (i in 1:length(bagOfPatterns)) {
    oldMap <- bagOfPatterns[[i]][["bob"]]
    bagOfPatterns[[i]][["bob"]] <- list()
    for (key in names(oldMap)) {
      if (oldMap[[key]] > 0) {
        catch <- getWordChi(dictionary, as.numeric(key))
        dictionary <- catch[[1]]
        wordChi <- catch[[2]]
        bagOfPatterns[[i]][["bob"]][[as.character(wordChi)]] <- oldMap[[key]]
      }
    }
  }
  value <- list()
  value[[1]] <- dictionary
  value[[2]] <- bagOfPatterns
  return(value)
}

# searches for a given entry in a dictionary
getKey <- function(dictionary, key) {
  if (as.character(key) %in% names(dictionary)) 
    value = dictionary[[as.character(key)]]
  else value = 0
  return(value)
}
