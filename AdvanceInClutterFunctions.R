nncleanEngine <- function (X, k, d, ..., tol = 0.001, maxit = 50, plothist = FALSE, 
                           lineargs = list(), verbose = TRUE, Xname = "X") 
{
  kthNND <- nndist(X, k = k)
  n <- length(kthNND)
  if (k >= n) {
    if (verbose) 
      cat(paste("Cannot compute neighbours of order k =", 
                k, "for a pattern of", n, "data points;", "treating all points as noise"), 
          call. = FALSE)
    return(list(z = rep(0, n), probs = rep(0, n), lambda1 = NA, 
                lambda2 = NA, p = 0, kthNND = kthNND, d = d, n = n, 
                k = k, niter = 0, maxit = maxit, converged = TRUE, 
                hist = NULL))
  }
  d <- ensure2vector(d)
  alpha.d <- (2 * pi^(d/2))/(d * gamma(d/2))
  kNNDpowd1 <- kthNND^(d[1])
  kNNDpowd2 <- kthNND^(d[2])
  probs <- numeric(n)
  thresh <- (min(kthNND) + diff(range(kthNND))/3)
  high <- (kthNND > thresh)
  delta <- as.integer(high)
  p <- 0.5
  lambda1 <- k/(alpha.d[1] * mean(kNNDpowd1[!high]))
  lambda2 <- k/(alpha.d[2] * mean(kNNDpowd2[high]))
  loglik.old <- 0
  loglik.new <- 1
  Z <- !kthNND
  niter <- 0
  while (abs(loglik.new - loglik.old)/(1 + abs(loglik.new)) > 
         tol) {
    if (niter >= maxit) {
      warning(paste("E-M algorithm failed to converge in", 
                    maxit, ngettext(maxit, "iteration", "iterations")), 
              call. = FALSE)
      break
    }
    niter <- niter + 1
    f1 <- dknn(kthNND[!Z], lambda = lambda1, k = k, d = d[1])
    f2 <- dknn(kthNND[!Z], lambda = lambda2, k = k, d = d[2])
    delta[!Z] <- (p * f1)/(p * f1 + (1 - p) * f2)
    delta[Z] <- 0
    sumdelta <- sum(delta)
    negdelta <- 1 - delta
    p <- sumdelta/n
    lambda1 <- (k * sumdelta)/(alpha.d[1] * sum(kNNDpowd1 * 
                                                  delta))
    lambda2 <- (k * (n - sumdelta))/(alpha.d[2] * sum(kNNDpowd2 * 
                                                        negdelta))
    loglik.old <- loglik.new
    loglik.new <- sum(-p * lambda1 * alpha.d[1] * (kNNDpowd1 * 
                                                     delta) - (1 - p) * lambda2 * alpha.d[2] * (kNNDpowd2 * 
                                                                                                  negdelta) + delta * k * log(lambda1 * alpha.d[1]) + 
                        negdelta * k * log(lambda2 * alpha.d[2]))
    if (verbose) 
      cat(paste("Iteration", niter, "\tlogLik =", loglik.new, 
                "\tp =", signif(p, 4), "\n"))
  }
  if (plothist) {
    dotargs <- list(...)
    if (spatstat.options("monochrome")) 
      dotargs <- col.args.to.grey(dotargs)
    xlim <- c(0, max(kthNND))
    H <- do.call(hist, resolve.defaults(list(quote(kthNND), 
                                             plot = FALSE, warn.unused = FALSE), dotargs, list(nclass = 40)))
    barheights <- H$density
    support <- seq(from = xlim[1], to = xlim[2], length.out = 200)
    fittedy <- p * dknn(support, lambda = lambda1, k = k, 
                        d = d[1]) + (1 - p) * dknn(support, lambda = lambda2, 
                                                   k = k, d = d[2])
    ylim <- range(c(0, barheights, fittedy))
    xlab <- paste("Distance to", ordinal(k), "nearest neighbour")
    reallyplot <- resolve.1.default("plot", list(...), list(plot = TRUE))
    H <- do.call(hist, resolve.defaults(list(quote(kthNND), 
                                             probability = TRUE), dotargs, list(plot = TRUE, warn.unused = reallyplot, 
                                                                                nclass = 40, xlim = xlim, ylim = ylim, xlab = xlab, 
                                                                                ylab = "Probability density", axes = TRUE, main = "")))
    H$xname <- xlab
    if (reallyplot) {
      box()
      lineargs <- resolve.defaults(lineargs, list(col = "green", 
                                                  lwd = 2))
      if (spatstat.options("monochrome")) 
        lineargs <- col.args.to.grey(lineargs)
      do.call(lines, append(list(x = support, y = fittedy), 
                            lineargs))
    }
  }
  delta1 <- dknn(kthNND[!Z], lambda = lambda1, k = k, d = d[1])
  delta2 <- dknn(kthNND[!Z], lambda = lambda2, k = k, d = d[2])
  probs[!Z] <- delta1/(delta1 + delta2)
  probs[Z] <- 1
  if (verbose) {
    cat("Estimated parameters:\n")
    cat(paste("p [cluster] =", signif(p, 5), "\n"))
    cat(paste("lambda [cluster] =", signif(lambda1, 5), "\n"))
    cat(paste("lambda [noise]   =", signif(lambda2, 5), "\n"))
  }
  return(list(delta1 = delta1, delta2 = delta2, z = round(probs), probs = probs, lambda1 = lambda1, 
              lambda2 = lambda2, p = p, kthNND = kthNND, d = d, n = n, 
              k = k, niter = niter, maxit = maxit, converged = (niter >= 
                                                                  maxit), hist = if (plothist) H else NULL))
}


class_net_Engine <- function(X, K){
  
  em_n <- nncleanEngine(X = X, k = K, d = 2, verbose = F) 
  pp_n <- em_n$probs 
  pp_n[which(pp_n == 0)] <- 0.000000001
  entr <-  - sum(pp_n * log2(pp_n))
  zz_n <- em_n$z
  zz_n <- factor(zz_n, levels = c(0, 1))
  levels(zz_n) <- c("noise", "feature")
  clas <- as.ppp(unmark(X))
  clas$marks <- zz_n
  
  if(is.multitype(X)){
    out <- list(rates = calculate_statistics(table(clas$marks, X$marks)),
                entropy = entr, X = X, clas = clas, K = K)
  } else {
    out <- list(entropy = entr, X = X, clas = clas, K = K)
  }
  
  class(out) <- "localdetectE"
  return(out)
  
}


calculate_statistics <- function(tbl) {
  P <- sum(tbl[,2])
  N <- sum(tbl[,1])
  TP <- tbl[2,2]
  TN <- tbl[1,1]
  FP <- tbl[2,1]
  FN <- tbl[1,2]
  TPR <- TP / P
  FPR <- FP / N
  ACC <- (TP + TN) / (P + N)
  data.frame(TPR = TPR, FPR = FPR, ACC = ACC)
}



class_net <- function(X, K = NULL, n_feat = 1, n_entr = 5:35, verbose = T){
  
  if(n_feat > 4) stop("n_feat > 4 not allowed")
  
  X0 <- X
  K0 <- K
  
  Z <-  if(is.multitype(X)) unmark(X) else X
  if(is.multitype(X)) Z$marks <- factor(1:npoints(X))
  
  class_track <- list()
  K_track <- vector(l = n_feat)
  
  deltas_track <- list()
  
  for(n in 1:n_feat){
    
    if(verbose) cat("Classifying feature at iteration", n, "\n", "\n")
    
    
    deltas <- vector()
    Z0 <-  if(is.multitype(Z)) unmark(Z) else Z
    for(i in n_entr){
      
      deltas[i] <- class_net_Engine(Z0, K = i)$entropy
    }
    
    deltas <- deltas[complete.cases(deltas)]
    if(is.null(K0)){
      o_track <- list()
      if(verbose)  cat("Selecting K ... ")
      
      
      x0 <- n_entr
      z0 <- - x0
      out.lm <- lm(deltas ~ 1)
      o <- try(segmented(out.lm, ~ z0), silent = T) 
      K <- - round(o$psi[2])
      
      o_track[[n]] <- o
      
      
      if(verbose) cat(paste("Estimated K =", K, "\n", "\n"))
    } 
    deltas_track[[n]] <- deltas
    
    kthNND_n <- nndist(Z, k = K)

    em_n <- nncleanEngine(X = Z, k = K, d = 2, verbose = F) 
    pp_n <- em_n$probs 
    entr <-  - sum(pp_n * log2(pp_n))
    zz_n <- em_n$z
    zz_n <- factor(zz_n, levels = c(0, 1))
    levels(zz_n) <- c("noise", "feature")
    
    clas <- if(is.multitype(X)){
      ppp(Z$x, Z$y, owin(range(Z$x), range(Z$y)), marks = data.frame(zz_n, Z$marks))
    } else {
      ppp(Z$x, Z$y, owin(range(Z$x), range(Z$y)), marks = data.frame(zz_n))
    }
    
    
    class_track[[n]] <- clas
    K_track[n] <- K
    
    Z <- if(is.multitype(X)){
      clas[clas$marks[, 1] == "feature"]
    } else {
      clas[clas$marks == "feature"]
    }
    
    if(n == n_feat){
      id <- if(is.multitype(X)){
        as.numeric(Z$marks$Z.marks)
      } else {
        as.numeric(Z$marks)
      }
      v <- vector(l = npoints(X0))
      v[id] <- "feature"
      v[ - id] <- "clutter"
    } 
    
  }
  
  if(verbose) cat(paste("Procedure ended. \n"))
  
  if(is.null(K0)){
    out <- list(n_feat = n_feat, K_track = K_track, K0 = K0,
                X = X0, id = v, o_track = o_track,
                class_track = class_track, n_entr = n_entr,
                deltas_track = deltas_track)
    
  } else{
    out <- list(n_feat = n_feat, K_track = K_track, K0 = K0, 
                X = X0, id = v,
                class_track = class_track, n_entr = n_entr,
                deltas_track = deltas_track)
  }
  class(out) <- "localdetect"
  return(out)
  
}


rates.localdetect <- function(x){
  if(is.multitype(x$X)){
    calculate_statistics(table(x$id, x$X$marks))
  } else {
    stop("No rates to compute, as no true classification is known")
  }
}


