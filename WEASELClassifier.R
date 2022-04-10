source("linearModel.R")
source("timeSeries.R")
source("MFT.R")
source("SFA.R")
source("WEASEL.R")

createClassifier <- function(){
  classifier <- list()
  classifier[["maxF"]] <- 6
  classifier[["minF"]] <- 4
  classifier[["maxS"]] <- 4
  classifier[["chi"]] <- 2
  classifier[["bias"]] <- 1
  classifier[["p"]] <- 0.1
  classifier[["iter"]] <- 5000
  classifier[["c"]] <- 1
  classifier[["maxWindowLength"]] <- 20
  classifier[["wordModel"]] <- NULL
  classifier[["linearModel"]] <- NULL
  return(classifier)
}

createWEASELModel <- function(myNorm, f, correct, size, nFeatures) {
  model <- list()
  model[["norm"]] <- myNorm
  model[["f"]] <- f
  model[["correct"]] <- correct
  model[["size"]] <- size
  model[["nFeatures"]] <- nFeatures
  return(model)
}

# CLASSIFIER

# runs the WEASEL algorithm
evaluate <- function(classifier, train, test) {
  catch <- list()
  catch <- fitWEASEL(classifier, train)
  classifier <- catch[[1]]
  scores <- catch[[2]]
  catch <- prediction(classifier, scores, test)
  classifier <- catch[[1]]
  accuracy <- catch[[2]]
  labels <- catch[[3]]
  
  print("---- RESULTS ----")
  info <- paste("Trained correctly: ", round(100 * scores[["correct"]] / train[["samples"]], 2), 
                "%, classified correctly: ", round(100 * accuracy, 2), "%")
  print(info)
  sink()
  print("---- WORK FINISHED ----")
  return(classifier)
}

# builds and trains a classifier
fitWEASEL <- function(classifier, train) {
  maxCorrect = -1
  bestF = -1
  bestNorm = FALSE
  classifier[["minWindowLength"]] <- 4
  maxWindowLength = classifier[["maxWindowLength"]]
  for (i in 1:train[["samples"]]) {
    maxWindowLength = min(length(train[["tseries"]][[i]][["data"]]), maxWindowLength)
  }
  classifier[["windows"]] <- classifier[["minWindowLength"]]:(maxWindowLength - 1)
  
  keepGoing = TRUE
  catch <- list()
  for (normMean in c(TRUE, FALSE)) {
    print("")
    if (normMean == TRUE) print("---- FITTING FOR NORM_MEAN = TRUE ----")
    else print("---- FITTING FOR NORM_MEAN = FALSE ----")
    
    if (keepGoing) {
      weasel <- createWEASEL(classifier[["maxF"]], classifier[["maxS"]], 
                             classifier[["windows"]], normMean)
      catch <- createFitting(weasel, train)
      weasel <- catch[[1]]
      words <- catch[[2]]
      f = classifier[["minF"]]
      
      while (f <= classifier[["maxF"]] & keepGoing) {
        weasel[["dict"]] <- resetDictionary(weasel[["dict"]])
        catch <- createBOP(weasel, words, train, f)
        weasel <- catch[[1]]
        bop <- catch[[2]]
        catch <- filterChiSquared(weasel, bop, classifier["chi"])
        weasel <- catch[[1]]
        bop <- catch[[2]]
        problem <- initializeLibLinearProblem(classifier, bop, weasel[["dict"]],
                                              classifier[["bias"]])
        correct <- trainLibLinear(classifier, problem, 10)
        if (correct > maxCorrect) {
          maxCorrect = correct
          bestF = f
          bestNorm = normMean
        }
        if (correct == train[["samples"]]) keepGoing = FALSE
        f = f + 2
      }
    }
  }
  
  print("")
  print("---- BUILDING A CLASSIFIER ----")
  catch <- list()
  classifier[["wordModel"]] <- createWEASEL(classifier[["maxF"]], classifier[["maxS"]], 
                                            classifier[["windows"]], bestNorm)
  catch <- createFitting(classifier[["wordModel"]], train)
  classifier[["wordModel"]] <- catch[[1]]
  words <- catch[[2]]
  
  catch <- createBOP(classifier[["wordModel"]], words, train, bestF)
  classifier[["wordModel"]] <- catch[[1]]
  bop <- catch[[2]]
  
  catch <- filterChiSquared(classifier[["wordModel"]], bop, classifier[["chi"]])
  classifier[["wordModel"]] <- catch[[1]]
  bop <- catch[[2]]
  
  problem <- initializeLibLinearProblem(classifier, bop, classifier[["wordModel"]][["dict"]], 
                                        classifier[["bias"]])
  param <- createParameter(classifier[["c"]], classifier[["iter"]], 
                           classifier[["p"]])

  classifier[["linearModel"]] <- trainModel(problem, param)
  myModel <- createWEASELModel(bestNorm, bestF, maxCorrect, train[["samples"]], 
                               problem[["n"]])
  
  value <- list()
  value[[1]] <- classifier
  value[[2]] <- myModel
  return(value)
}

# uses a classifier to test data
prediction <- function(classifier, scores, test) {
  print("")
  print("---- TESTING A CLASSIFIER ----")
  catch <- list()
  catch <- createFitting(classifier[["wordModel"]], test)
  classifier[["wordModel"]] <- catch[[1]]
  words <- catch[[2]]
  
  catch  <- createBOP(classifier[["wordModel"]], words, test, scores[["f"]])
  classifier[["wordModel"]] <- catch[[1]]
  bop <- catch[[2]]
  
  catch <- remap(classifier[["wordModel"]][["dict"]], bop)
  classifier[["wordModel"]][["dict"]] <- catch[[1]]
  bop <- catch[[2]]
  classifier[["features"]] <- initializeLibLinear(bop, scores[["nFeatures"]])
  
  predLabels <- vector()
  for (f in classifier[["features"]]) {
    predLabels <- append(predLabels, predictModel(classifier[["linearModel"]], f))
  }
  acc = 0
  for (i in 1:test[["samples"]]) {
    acc = acc + as.numeric(predLabels[[i]] == test[["labels"]][[i]])
  }
  acc = acc / test[["samples"]]
  value <- list()
  
  print("Printing output...")
  sink("output.txt")
  print("---- TEST OUTPUT ----")
  print(classifier)
  print(predLabels)
  print(acc)
  value[[1]] <- classifier
  value[[2]] <- acc
  value[[3]] <- predLabels
  return(value)
}

# LIBLINEAR

initializeLibLinearProblem <- function(classifier, bop, dict, bias) {
  n <- getSize(dict)
  y <- vector("list", length = length(bop))
  for (i in 1:length(bop)) {
    y[[i]] <- bop[[i]][["label"]]
  }
  features <- initializeLibLinear(bop, n)
  prob <- createProblem(length(features), n, y, features, bias)
  return(prob)
}

# NONE
# sorts all features from boxes of patterns
initializeLibLinear <- function(bop, maxFeature) {
  featuresTrain <- vector("list", length = length(bop))
  for (b in 1:length(bop)) {
    features <- list()
    bob <- bop[[b]]
    index = 0
    for (key in names(bob[["bob"]])) {
      if (bob[["bob"]][[key]] > 0  & as.numeric(key) <= maxFeature) {
        index = index + 1
        features[[index]] <- createFeatureNode(as.numeric(key), bob[["bob"]][[key]])
      }
    }
    myList <- vector("list", length = 2)
    if (length(features) > 0) {
      for (i in 1:length(features)) {
        myList[[1]] <- append(myList[[1]], features[[i]][["index"]])
        myList[[2]] <- append(myList[[2]], features[[i]][["value"]])
      }
    }
    myFrame <- data.frame(myList)
    if (length(features) > 0) {
      myFrame <- myFrame[order(myFrame[[1]]),]
    }
    
    newFeature <- list()
    if (length(features) > 0) {
      for (i in 1:nrow(myFrame)) {
        newFeature[[i]] <- createFeatureNode(myFrame[i, 1], myFrame[i, 2])
      }
    }
    featuresTrain[[b]] <- newFeature
  }
  return(featuresTrain)
}

# randomzes a training sample, creates 10 folds (subproblems) in which 
# predicts labels
trainLibLinear <- function(classifier, prob, nFolds = 10) {
  param <- createParameter(classifier[["c"]], classifier[["iter"]], classifier[["p"]])
  if (nFolds > prob[["l"]]) nFolds = prob[["l"]]
  startFold <- vector()
  startFold[[1]] <- 0
  perm <- vector(length = prob[["l"]])
  for (i in 1:prob[["l"]]) {
    perm[i] = i
  }
  perm <- sample(perm)
  for (i in 2:nFolds) {
    startFold[[i]] <- floor((i - 1) * prob[["l"]] / nFolds)
  }
  startFold <- append(startFold, prob[["l"]])
  correct = 0
  
  # cross validation of training set
  for (i in 1:nFolds) {
    b = startFold[i]
    e = startFold[i + 1]
    subProb <- createProblem(prob[["l"]] - (e - b), prob[["n"]], NULL, NULL, 
                             prob[["bias"]])
    rows <- vector()
    subProb[["y"]] <- list()
    if (b > 0) {
      for (j in 1:b) {
        rows <- append(rows, perm[j])
        subProb[["y"]] <- append(subProb[["y"]], prob[["y"]][[perm[j]]])
      }
    }
    if (prob[["l"]] > e)
    {
      for (j in (e + 1):prob[["l"]]) {
        rows <- append(rows, perm[j])
        subProb[["y"]] <- append(subProb[["y"]], prob[["y"]][[perm[j]]])
      }
    }
    subProb[["x"]] <- list()
    index = 0
    for (j in rows) {
      index = index + 1
      subProb[["x"]][[index]] <- prob[["x"]][[j]]
    }
    
    foldModel <- trainModel(subProb, param)
    foldX <- vector("list", length = e - b)
    foldY <- vector()
    index = 0
    for (j in (b + 1):e) {
      index = index + 1
      foldX[[index]] <- prob[["x"]][[perm[j]]]
      foldY <- append(foldY, prob[["y"]][[perm[j]]])
    }
    foldLabels <- vector(length = length(foldY))
    for (j in 1:length(foldY)) {
      foldLabels[j] <- predictModel(foldModel, foldX[[j]])
      message <- paste("Predicted model no. ", b + j, " of ", length(perm))
      print(message)
    }
    for (j in 1:length(foldY)) {
      if (foldY[j] == foldLabels[j]) correct = correct + 1
    }
  }
  return(correct)
}
