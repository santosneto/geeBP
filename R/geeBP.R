print.geeBP <- function(x, digits = NULL, quote = FALSE, prefix = "", ...){
  if(is.null(digits)) digits <- options()$digits else options(digits =
                                                                digits)
  cat("\n", x$title)
  cat("n"," Model:\n")
  cat(" Link:                     ", x$model$link, "\n")
  cat(" Variance to Mean Relation:", x$model$varfun, "\n")
  cat(" Correlation Structure:    ", x$model$corstr, "\n")
  cat("\nCall:\n")
  dput(x$call)                        #       cat("\nTerms:\n")
  ###        ys <- matrix(rep(as.matrix(x$id, ncol = 1), 5), ncol = 5)
  ys <- matrix(rep(matrix(x$id, ncol = 1), 5), ncol = 5)
  ys[, 2] <- x$y
  ys[, 3] <- x$linear.predictors
  ys[, 4] <- x$fitted.values
  ys[, 5] <- x$residuals
  dimnames(ys) <- list(1:length(x$y), c("ID", "Y", "LP", "fitted",
                                        "Residual")) #       cat("\nFitted Values:\n")
  cat("\nNumber of observations : ", x$nobs, "\n")
  cat("\nMaximum cluster size   : ", x$max.id, "\n")
  cat("\n\nCoefficients:\n")
  print(t(x$coefficients), digits = digits)
  cat("\nEstimated Scale Parameter: ", format(round(x$scale, digits)))
  cat("\nNumber of Iterations: ", x$iterations)
  cat("\n\nEstimated Correlation Parameters:\n")
  print(x$alpha, digits = digits)
  invisible(x)
}

summary.geeBP <- function(object,...){
  value = list()
  class(value) = "summary.geeBP"
  value$call = object$call
  value$scale = object$scale
  value$corstr = object$model$corstr
  value$link = object$model$link
  value$alpha = object$alpha
  value$nclusters = object$nclusters
  value$max.id = object$max.id
  valorp = 0
  wald = 0
  for(i in 1:length(object$coefficients)){
    wald[i] = wald.test(object$robust.variance, object$coefficients, Terms = i)$result$chi2[1]
    valorp[i] = wald.test(object$robust.variance, object$coefficients, Terms = i)$result$chi2[3]
  }
  value$coefficients = data.frame(object$coefficients, object$robust.se, wald, valorp)
  colnames(value$coefficients) = c("Estimate","Robust.Std.err", "Wald", "Pr(>|W|)")
  return(value)
}

print.summary.geeBP <- function(x, digits = NULL, quote = FALSE, prefix = "", ...){
  if(is.null(digits)) digits <- options()$digits else options(digits = digits)
  cat("\nCall:\n")
  dput(x$call)                        #       cat("\nTerms:\n")
  cat("\nCoefficients:\n")
  printCoefmat(as.matrix(x$coefficients), digits = digits)
  cat("\nNumber of clusters : ", x$nclusters, "\n")
  cat("\nMaximum cluster size   : ", x$max.id, "\n")
  cat("\nEstimated Scale Parameter: ", format(round(x$scale, digits)))
  cat("\n\nEstimated Correlation Parameters:\n")
  print(x$alpha, digits = digits)
  invisible(x)
}

residuals.geeBP <- function (object, type = c("standard", "pearson", "response"), ...)
{
  type <- match.arg(type)
  y   <- object$y
  r   <- object$residuals
  mu  <- object$fitted.values
  X = object$comp$X
  W = object$comp$W
  u = object$comp$u
  sqrtW = sqrtm(W)$B
  slam = solve(object$comp$Lambda)
  H = sqrtW%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%sqrtW
  hij = diag(H)
  rij = 0
  for(i in 1:nrow(X)){
    mij = rep(0,nrow(X))
    mij[i] = 1
    rij[i] = (t(mij)%*%sqrtW%*%slam%*%u)/sqrt(1-hij[i])
  }
  res <- switch(type,
                pearson = r,
                response = y - mu,
                standard = rij)
  if (!is.null(object$na.action))
    res <- naresid(object$na.action, res)
  #    if (type == "partial")
  #      res <- res + predict(object, type = "terms")
  res
}

plot.geeBP <- function(x,...){
  xx <- fitted(x)
  rp <- residuals(x,"standard")
  qplot(xx,rp) + ylab("Standard residuals") + xlab("Fitted values")
}

geeBP = function(formula, data, id, tol = 0.001, maxiter = 25,
                 corstr = "independence", linkmu = "log"){

  namescor = c("independence", "unstructured", "exchangeable", "AR-1", "one-dependent",
               "one-dependent-stat","two-dependent","two-dependent-stat")
  if(all(namescor != corstr)){
    stop("the correlation structure is not defined")
  }
  nameslink = c("log", "identity")
  if(all(nameslink != linkmu)){
    stop("the link function is not defined")
  }
  formula = as.formula(formula)
  nformula = all.vars(formula)
  fnames = 0
  jaux = 1
  listaux = list(NULL)
  for(i in 1:length(nformula)){
    if(is.factor(data[,nformula[i]])){
      fnames[jaux] = nformula[i]
      listaux[[jaux]] = "contr.sum"
      jaux = jaux+1
    }
  }
  call <- match.call()
  if(jaux>1){
    names(listaux) = fnames
    X = as.matrix(model.matrix(formula, data = data, contrasts = listaux))
  }
  else{
    X = as.matrix(model.matrix(formula, data = data))
  }
  p = ncol(X)
  y = model.frame(formula, data = data)[,1]
  t = as.vector(table(id))
  n = length(table(id))
  N = nrow(X)
  if(linkmu == "log"){
    if(jaux>1){
      mod0 = gamlss(formula,family = BP(mu.link = "log"), trace = FALSE, data = data, contrasts = listaux)
    }
    else{
      mod0 = gamlss(formula,family = BP(mu.link = "log"), trace = FALSE, data = data)
    }
  }
  if(linkmu == "identity"){
    if(jaux>1){
      mod0 = gamlss(formula,family = BP(mu.link = "log"), trace = FALSE, data = data, contrasts = listaux)
    }
    else{
      mod0 = gamlss(formula,family = BP(mu.link = "log"), trace = FALSE, data = data)
    }
  }
  if(corstr == "independence"){
    cat("a gamlss object was returned")
    return(mod0)
  }
  beta = mod0$mu.coefficients
  phi = mod0$sigma.coefficients
  if(phi<0){
    phi = 1
  }

  cont = 1
  repeat{
    eta = X%*%beta
    if(linkmu == "log"){
      mu = as.vector(exp(eta))
    }
    if(linkmu == "identity"){
      mu = as.vector(eta)
    }
    vmu = mu^2 + mu

    # y estrela
    ys = log(y/(y+1))

    # mi estrela
    mus = digamma(mu*(phi+1))-digamma(mu*(phi+1)+phi+2)

    vmus = trigamma(mu*(phi+1))-trigamma(mu*(phi+1)+phi+2)

    # Vetor b_i
    u = ys-mus

    if(linkmu == "log"){
      G = diag(as.vector(mu))
    }
    if(linkmu == "identity"){
      G = diag(1,N,N)
    }
    A = diag(as.vector(vmus))
    Lambda = (phi+1)*G%*%A

    uc = split(u,id)
    if(corstr == "unstructured"){
      Rg = matrix(0,max(t),max(t))
      cnum = den1 = den2 = 0
      for(j in 1:(max(t))){
        for(k in j:(max(t))){
          for(i in 1:n){
            if(is.na(uc[[i]][j])||is.na(uc[[i]][k])){
              cnum = cnum
              den1 = den1
              den2 = den2
            }
            else{
              cnum = cnum + (uc[[i]][j])*(uc[[i]][k])
              den1 = den1 + (uc[[i]][j])^2
              den2 = den2 + (uc[[i]][k])^2
            }
          }
          Rg[j,k] = cnum/(sqrt(den1)*sqrt(den2))
          Rg[k,j] = Rg[j,k]
        }
      }
      diag(Rg) = 1
      R = list(NULL)
      for(i in 1:n){
        R[[i]] = Rg[1:t[i],1:t[i]]
      }
      Rm = bdiag(R)
    }
    if(corstr == "AR-1"){
      cnum = cden1 = cden2 = 0
      for(i in 1:n){
        for(j in 1:(t[i]-1)){
          cnum = cnum + uc[[i]][j]*uc[[i]][j+1]
          cden1 = cden1 + (uc[[i]][j]^2)
          cden2 = cden2 + (uc[[i]][j+1]^2)
        }
      }
      alpha = cnum/sqrt(cden1*cden2)
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      R = list(NULL)
      for(i in 1:n){
        R[[i]] = matrix(0,t[i],t[i])
        for(j in 1:t[i]){
          for(l in 1:t[i]){
            R[[i]][j,l] = alpha^(abs(j-l))
          }
        }
      }
      Rm = as.matrix(bdiag(R))
      R=R[[1]]
    }
    if(corstr == "exchangeable"){
      cnum = cden = 0
      for(i in 1:n){
        aux = uc[[i]]%*%t(uc[[i]])
        cnum = cnum + sum(aux[upper.tri(aux)])*(2/(t[i]-1))
        for(j in 1:(t[i])){
          cden = cden + (uc[[i]][j]^2)
        }
      }
      alpha = (cnum/cden)
      Rm = matrix(0,N,N)
      R = list(NULL)
      for(i in 1:n){
        R[[i]] = matrix(alpha,t[i],t[i])
        diag(R[[i]]) = 1
      }
      Rm = as.matrix(bdiag(R))
      R=R[[1]]
    }
    if(corstr == "one-dependent"){
      alpha = 0
      den = 0
      for(i in 1:n){
        for(j in 1:t[1]){
          den = den + (scomb[i,j]^2)
        }
      }
      for(j in 1:(t[1]-1)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+1])
        }
        alpha[j] = num/den
      }
      alpha = (N/n)*alpha
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      R = matrix(0,t[1],t[1])
      for(i in 1:t[1]){
        for(j in 1:t[1]){
          if(j==(i+1)){
            R[i,j] = alpha[i]
            R[j,i] = R[i,j]
          }
        }
      }
      diag(R) = 1
      Rm = kronecker(diag(n),R)
    }
    if(corstr == "two-dependent"){
      alpha1 = 0
      den = 0
      for(i in 1:n){
        for(j in 1:t[1]){
          den = den + (scomb[i,j]^2)
        }
      }
      for(j in 1:(t[1]-1)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+1])
        }
        alpha1[j] = num/den
      }
      alpha1 = (N/n)*alpha1
      alpha2 = 0
      for(j in 1:(t[1]-2)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+2])
        }
        alpha2[j] = num/den
      }
      alpha2 = (N/n)*alpha2
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      R = matrix(0,t[1],t[1])
      for(i in 1:t[1]){
        for(j in 1:t[1]){
          if(j==(i+1)){
            R[i,j] = alpha1[i]
            R[j,i] = R[i,j]
          }
          if(j == (i+2)){
            R[i,j] = alpha2[i]
            R[j,i] = R[i,j]
          }
        }
      }
      diag(R) = 1
      Rm = kronecker(diag(n),R)
    }
    if(corstr == "one-dependent-stat"){
      alpha = 0
      den = 0
      for(i in 1:n){
        for(j in 1:t[1]){
          den = den + (scomb[i,j]^2)
        }
      }
      for(j in 1:(t[1]-1)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+1])
        }
        alpha[j] = num/den
      }
      alpha = (N/n)*alpha
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      R = matrix(0,t[1],t[1])
      for(i in 1:t[1]){
        for(j in 1:t[1]){
          if(j==(i+1)){
            R[i,j] = sum(alpha)/(t[1]-1)
            R[j,i] = R[i,j]
          }
        }
      }
      diag(R) = 1
      Rm = kronecker(diag(n),R)
    }
    if(corstr == "two-dependent-stat"){
      alpha1 = 0
      den = 0
      for(i in 1:n){
        for(j in 1:t[1]){
          den = den + (scomb[i,j]^2)
        }
      }
      for(j in 1:(t[1]-1)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+1])
        }
        alpha1[j] = num/den
      }
      alpha1 = (N/n)*alpha1
      alpha2 = 0
      for(j in 1:(t[1]-2)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+2])
        }
        alpha2[j] = num/den
      }
      alpha2 = (N/n)*alpha2
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      R = matrix(0,t[1],t[1])
      for(i in 1:t[1]){
        for(j in 1:t[1]){
          if(j==(i+1)){
            R[i,j] = sum(alpha1)/(t[1]-1)
            R[j,i] = R[i,j]
          }
          if(j == (i+2)){
            R[i,j] = sum(alpha2)/(t[1]-1)
            R[j,i] = R[i,j]
          }
        }
      }
      diag(R) = 1
      Rm = kronecker(diag(n),R)
    }
    Omega = as.matrix(sqrt(A)%*%Rm%*%sqrt(A))
    W = Lambda%*%solve(Omega)%*%t(Lambda)
    z = eta + solve(Lambda)%*%u

    #Novo valor de beta
    beta1 = solve(t(X)%*%W%*%X)%*%(t(X)%*%W%*%z)
    dif = abs(beta1-beta)
    if(sum(dif)<=(tol*p)){
      beta = beta1
      cat("The algorithm converged")
      converg = 1
      break
    }

    if(cont == maxiter){
      cat("Maximum number of iterations reached")
      converg = 0
      break
    }
    beta = beta1

    r = (y-mu)*(1/sqrt(vmu))

    phi = 1/(sum(r^2)/(N-p))
    cont = cont + 1
  }

    # Matriz de sensibilidade
    S = -t(X)%*%W%*%X
    invOmega = solve(Omega)

    VarBeta = solve(S)%*%t(X)%*%Lambda%*%invOmega%*%u%*%t(u)%*%invOmega%*%Lambda%*%X%*%solve(S)

  SEbeta = sqrt(diag(VarBeta))

  fit <- list()
  attr(fit, "class") <- c("geeBP")
  fit$title <- "geeBP:  BETA PRIME GENERALIZED ESTIMATING EQUATIONS"
  fit$model <- list()
  fit$model$link <- linkmu
  fit$model$varfun <- "mu(1+mu)"
  fit$model$corstr <- corstr
  fit$call <- call
  fit$formula <- formula
  fit$nclusters = n
  fit$clusters = t
  fit$nobs <- N
  fit$iterations <- cont
  fit$coefficients <- beta
  eta <- as.vector(X %*% fit$coefficients)
  fit$linear.predictors <- eta
  mu <- as.vector(mu)
  fit$fitted.values <- mu
  fit$residuals <- r
  fit$family <- "Beta prime"
  fit$y <- as.vector(y)
  fit$id <- as.vector(id)
  fit$max.id <- max(t)
  fit$working.correlation <- R
  fit$scale <- phi
  fit$robust.variance <- VarBeta
  fit$robust.se = SEbeta
  if(corstr == "unstructured"){
    fit$alpha = fit$working.correlation[upper.tri(fit$working.correlation)]
  }
  if(corstr == "AR-1"||corstr == "exchangeable"||corstr == "one-dependent"){
    fit$alpha = alpha
  }
  if(corstr == "two-dependent"){
    fit$alpha = list(diag1 = alpha1, diag2 = alpha2)
  }
  if(corstr == "one-dependent-stat"){
    fit$alpha = R[1,2]
  }
  if(corstr == "two-dependent-stat"){
    fit$alpha = c(R[1,2],R[1,3])
  }
  fit$comp$X = X
  fit$comp$W = W
  fit$comp$u = u
  fit$comp$Lambda = Lambda
  fit$comp$vmu = vmu
  fit$comp$A = A
  fit$comp$G = G
  fit$comp$Rm = Rm
  fit$comp$Omega = Omega
  # QIC
  Q = phi*(y*(log(mu)-log(y))-(y+1)*(log(1+mu)-log(1+y)))
  psiq = sum(Q)
  VI = A
  oi = t(X)%*%Lambda%*%solve(VI)%*%Lambda%*%X
  QIC = -2*psiq + 2*sum(diag(oi%*%VarBeta))
  CIC = sum(diag(oi%*%VarBeta))
  fit$QIC = QIC
  fit$CIC = CIC
  devi = -2*(1/phi)*Q
  D = sum(devi)
  phia = D/n
  fit$EQIC = (phi)*D + sum(log(2*pi*(1/phi)*(diag(A)+1/6))) + 2*CIC
  return(fit)
}
