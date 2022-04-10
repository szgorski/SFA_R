createProblem <- function(l, n, y, x, bias = 0) {
  problem <- list()
  problem[["l"]] <- l
  problem[["n"]] <- n
  problem[["y"]] <- y
  problem[["x"]] <- x
  problem[["bias"]] <- bias
  return(problem)
}

createModel <- function(bias, label, nr_class, nr_feature, w) {
  myModel <- list()
  myModel[["bias"]] <- bias
  myModel[["label"]] <- label
  myModel[["nr_class"]] <- nr_class
  myModel[["nr_feature"]] <- nr_feature
  myModel[["w"]] <- w
  return(myModel)
  
}

createParameter <- function(c, iter, p) {
  param <- list()
  param[["c"]] <- c
  param[["iter"]] <- iter
  param[["p"]] <- p
  return(param)
}

# FEATURENODE

createFeatureNode <- function(index, value) {
  fn <- list()
  fn[["index"]] <- index
  fn[["value"]] <- value
  return(fn)
}

# NONE

# divides chosen time series to groups, divides a problem to subproblems
# with respect to groups in order to train a model on each class
trainModel <- function(problem, param) {
  print("Training a model...")
  n = length(problem[["x"]])
  indexBefore = 0
  for (i in 1:n) {
    nodes <- problem[["x"]][[i]]
    nr_class = length(nodes)
    if (nr_class > 0) {
      for (j in 1:nr_class) {
        node <- nodes[[j]]
        # if (node[["index"]] <= indexBefore) {}      # originally a return info
        # indexBefore = node[["index"]]
      }
    }
  }
  
  if (problem[["bias"]] >= 0) nr_feature = problem[["n"]] - 1
  else nr_feature = problem[["n"]]
  myModel <- createModel(problem[["bias"]], NULL, 0, nr_feature, NULL)
  
  perm <- rep(0, problem[["l"]])
  catch <- list()
  catch <- groupClasses(problem, perm)
  perm <- catch[[1]]
  nr_class <- catch[[2]]
  label <- catch[[3]]
  start <- catch[[4]]
  count <- catch[[5]]
  myModel[["nr_class"]] <- nr_class
  myModel[["label"]] <- label
  weightedC <- rep(param[["c"]], nr_class)
  
  x <- vector("list", length(perm))
  for (j in 1:length(perm)) {
    x[[j]] <- problem[["x"]][[perm[j]]]
  }
  
  sub_prob <- createProblem(problem[["l"]], problem[["n"]], NULL, NULL, 0)
  sub_prob[["x"]] <- vector("list", length = sub_prob[["l"]])
  for(i in 1:sub_prob[["l"]]) {
    sub_prob[["x"]][[i]] <- x[[i]]
  }
  sub_prob[["y"]] <- rep(0, sub_prob[["l"]])

  if (nr_class == 2) {
    myModel[["w"]] <- rep(0, problem[["n"]])
    i = start[1] + count[1]
    for (j in 1:i) {
      sub_prob[["y"]][j] = 1
    }
    i = i + 1
    while (i <= sub_prob[["l"]]) {
      sub_prob[["y"]][i] = -1
      i = i + 1
    }
    w = trainSingle(sub_prob, param, myModel[["w"]], weightedC[1], weightedC[2])
  } 
  else {
    myModel[["w"]] <- rep(0, problem[["n"]] * nr_class)
    w <- rep(0, length(myModel[["n"]]))
    for (i in 1:nr_class) {
      message <- paste("    Training class no. ", i, " of ", nr_class, "...")
      print(message)
      
      si = start[i]
      ei = si + count[i] - 1
      K = 1
      if (si > 1) {
        for (j in 1:(si - 1)) {
          sub_prob[["y"]][K] = -1
          K = K + 1
        }
      }
      while (K <= ei) {
        sub_prob[["y"]][K] = 1
        K = K + 1          
      }
      while (K <= sub_prob[["l"]]) {
        sub_prob[["y"]][K] = -1
        K = K + 1
      }
      w = trainSingle(sub_prob, param, w, weightedC[i], param[["c"]])
      for (j in 1:problem[["n"]]) {
        myModel[["w"]][(j - 1) * nr_class + i] = w[j]
      }
    }
  }
  return(myModel)
}

# calculates weights of each feature
trainSingle <- function(prob, param, w, Cp, Cn) {
  eps = param[["p"]]
  iter = 0
  xTx <- rep(0, prob[["l"]])
  maxIter = 1000
  index <- rep(0, prob[["l"]])
  alpha <- rep(0, 2 * prob[["l"]])
  y <- rep(0, prob[["l"]])
  maxInnerIter = 100
  innerEps = 0.01
  minInnerEps = min(10^(-8), eps)
  upperBound <- c(Cn, 0, Cp)
  
  for (i in 1:prob[["l"]]) {
    if (prob[["y"]][i] > 0) y[i] = 1
    else y[i] = -1
  }
  for (i in 1:prob[["l"]]) {
    alpha[2 * i - 1] = min(0.001 * upperBound[y[i] + 2], 10^(-8))
    alpha[2 * i] = upperBound[y[i] + 2] - alpha[2 * i]
  }
  for (i in 1:prob[["n"]]) {
    w[i] = 0
  }
  
  for (i in 1:prob[["l"]]) {
    if (length(prob[["x"]][[i]]) > 0) {
      for (j in 1:length(prob[["x"]][[i]])) {
        xi = prob[["x"]][[i]][[j]]
        C = xi[["value"]]
        xTx[i] = xTx[i] + C * C
        varIndex = xi[["index"]] - 1
        w[varIndex] = w[varIndex] + y[i] * alpha[2 * i - 1] * C
      }
  }
    index[i] = index[i] + i
  }
  
  while(iter < maxIter) {
    for (i in 1:prob[["l"]]) {
      newtonIter = i + sample(0:(prob[["l"]] - i - 1), 1)
      index <- swap(index, i, newtonIter)
    }
    newtonIter = 0
    maxG = 0
    for (s in 1:prob[["l"]]) {
      i = index[s]
      C = upperBound[y[i] + 2]
      ywTx = 0
      if (length(prob[["x"]][[i]]) > 0) {
        for (j in 1:length(prob[["x"]][[i]])) {
          xi = prob[["x"]][[i]][[j]]
          ywTx = ywTx + w[xi[["index"]]] * xi[["value"]]
        }
      }
      ywTx = ywTx * y[i]
      
      # a = xisq = xTx[i]
      # b = ywTx
      index1 = 2 * i - 1
      index2 = 2 * i
      sign = 1
      cond = 0.5 * xTx[i] * (alpha[index2] - alpha[index1]) * ywTx
      if (cond < 0) {
        index1 = 2 * i
        index2 = 2 * i - 1
        sign = -1
      }
      alphaOld = alpha[index1]
      z = alphaOld
      if (C - alphaOld < 0.5 * C) z = 0.1 * alphaOld
      gp = xTx[i] * (z - alphaOld) + sign * ywTx + log(z / (C - z))
      maxG = max(maxG, abs(gp))
      #eta = 0.1
      
      innerIter = 0
      while (innerIter <= maxInnerIter & abs(gp) >= innerEps) {
        gpp = xTx[i] + C / (C - z) / z
        z_temp = z - gp / gpp
        if (z_temp < 0) z = z * 0.1
        else z = z_temp
        gp = xTx[i] * (z - alphaOld) + sign * ywTx + log(z / (C - z))
        newtonIter = newtonIter + 1
        innerIter = innerIter + 1
      }
      if (innerIter > 0) {
        alpha[index1] = z
        alpha[index2] = C - z
        if (length(prob[["x"]][[i]]) > 0) {
          for (j in 1:length(prob[["x"]][[i]])) {
            xi = prob[["x"]][[i]][[j]]
            w[xi[["index"]]] = w[xi[["index"]]] + sign * (z - alphaOld) * y[i] * xi[["value"]]
          }
        }
      }
    }
    iter = iter + 1
    if (maxG < eps) break
    if (newtonIter <= prob[["l"]] / 10) innerEps = max(minInnerEps, 0.1 * innerEps)
  }
  
  # ordering of W is different than in JAVA
  # v = 0
  # for (i in prob[["n"]]) {
  #   v = v + w[i] * w[i]
  # }
  # v = 0.5 * v
  # for (i in 1:prob[["l"]]) {
  #   v = v + alpha[2 * i] * log(alpha[2 * i]) + alpha[2 * i + 1] * log(alpha[2 * i + 1]) - upperBound[y[i] + 2] * log(upperBound[y[i] + 2])
  # }
  return(w)
}

# gruops shuffled time series so that labels appear in order one after another
groupClasses <- function(prob, perm) {
  label <- vector()
  for (i in prob[["y"]]) {
    if (!(i %in% label)) label <- c(label, i)
  }
  count <- rep(0, length(label))
  dataLabel <- vector(length = length(prob[["y"]]))
  for (i in 1:length(prob[["y"]])) dataLabel[i] <- match(prob[["y"]][[i]], label)
  nr_class = length(label)
  
  for (i in 1:length(prob[["y"]])) 
    count[match(prob[["y"]][[i]], label)] = count[match(prob[["y"]][[i]], label)] + 1
  if (nr_class == 2 & label[1] == -1 & label[2] == 1) {
    label <- swap(label, 1, 2)
    count <- swap(count, 1, 2)
    for (i in 1:prob[["l"]]) {
      if (dataLabel[i] == 0) dataLabel[i] = 1
      else dataLabel[i] = 0
    }
  }
  
  start <- rep(0, nr_class)
  start[1] = 1
  for (i in 2:nr_class) start[i] = start[i - 1] + count[i - 1]
  for (i in 1:prob[["l"]]) {
    perm[start[dataLabel[i]]] = i
    start[dataLabel[i]] = start[dataLabel[i]] + 1
  }
  start[1] = 1
  for (i in 2:nr_class) start[i] = start[i - 1] + count[i - 1]
  
  value <- list()
  value[[1]] <- perm
  value[[2]] <- nr_class
  value[[3]] <- label
  value[[4]] <- start
  value[[5]] <- count
  return(value)
}

swap <- function(arr, a, b) {
  temp <- arr[[a]]
  arr[[a]] <- arr[[b]]
  arr[[b]] <- temp
  return(arr)
}

# copyOf <- function(original, newLength) {
#   copy <- vector(length = newLength)
#   len = min(length(original), newLength)
#   copy[1:len] <- original[1:len]
#   return(copy)
# }

# MODEL

# gives features, returns best-fitting labels
predictModel <- function(model, x) {
  values <- rep(0, model[["nr_class"]])
  if (model[["bias"]] >= 0) n = model[["nr_feature"]] + 1
  else n = model[["nr_feature"]]
  nr_w = model[["nr_class"]]
  
  if (length(x) > 0) {
    for (i in 1:length(x)) {
      index = x[[i]][["index"]]
      if (index < n) {
        for (j in 1:nr_w) {
          if ((x[[i]][["index"]] - 1) * nr_w + j <= length(model[["w"]]))
            values[j] = values[j] + model[["w"]][[(x[[i]][["index"]] - 1) * nr_w + j]] * x[[i]][["value"]]
        }
      }
    }
    if (model[["nr_class"]] == 2) {
      if (values[1] > 0) value = model[["label"]][[1]]
      else value = model[["label"]][[2]]
    }
    else {
      maxIdx = 1
      for (i in 2:model[["nr_class"]]) {
        if (values[i] > values[maxIdx]) maxIdx = i
      }
      value = model[["label"]][[maxIdx]]
    }
  }
  
  else {
    value = model[["label"]][[1]]
  }
  return(value)
}
