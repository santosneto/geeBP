diag.geeBP = function(model, type = c("cooks distance", "hat")){
  X = model$comp$X
  W = model$comp$W
  u = model$comp$u
  p = length(model$coefficients)
  sqrtW = sqrtm(W)$B
  slam = solve(model$comp$Lambda)
  H = sqrtW%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%sqrtW
  if(type == "hat"){
    return(H)
  }
  if(type == "cooks distance"){
    hij = diag(H)
    rij = residuals(model, "standard")
    CD = (rij^2)*(hij/(p*(1-hij)))
    return(CD)
  }
}

LI.geeBP = function(model, pert = c("case-weight", "response", "precision", "correlation",
                                    "covariate"), covariate = NULL,
                    type = c("normal", "standard")){
  X = model$comp$X
  W = model$comp$W
  u = model$comp$u
  sqrtW = sqrtm(W)$B
  slam = solve(model$comp$Lambda)

  if(pert == "case-weight"){
    Maux = diag(u)%*%slam%*%t(W)%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%slam%*%diag(u)
    resul = eigen(Maux, symmetric=T)$vectors[,1]
  }

  if(pert == "response"){
    phi = model$scale
    vmu = model$comp$vmu
    sij = sqrt(vmu/phi)
    y = model$y
    vy = as.vector(y*(1+y))
    B = diag(sij)%*%solve(diag(vy))
    Maux = B%*%slam%*%t(W)%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%slam%*%B
    resul = eigen(Maux, symmetric=T)$vectors[,1]
  }

  if(pert == "precision"){
    mu = model$fitted.values
    A = model$comp$A
    Rm = model$comp$Rm
    G = model$comp$G
    phi = model$scale
    Sigma = (sqrtm(A)$B)%*%Rm%*%(sqrtm(A)$B)
    invSigma = solve(Sigma)
    ws = diag(mu*psigamma(mu*(1+phi),2)-(1+mu)*psigamma((1+phi)*mu + phi + 2,2))
    wss = diag(mu*psigamma(mu*(1+phi),1)-(1+mu)*psigamma((1+phi)*mu + phi + 2,1))
    N = model$nobs
    Phi = diag(1+phi,N)
    Phi1 = Phi - diag(N)
    Aw = (sqrtm(A)$B)%*%((Phi1)^2)%*%wss
    dLam = -(1/phi)*G%*%((Phi1)^2)%*%(Phi%*%ws+A)
    dSigma = -0.5*(1/phi)*((sqrtm(A)$B)%*%Rm%*%Aw+Aw%*%Rm%*%(sqrtm(A)$B))
    dinvSigma = -invSigma%*%dSigma%*%invSigma
    du = (1/phi)*(Phi-1)^2%*%wss
    Delta = t(X)%*%(model$comp$Lambda)%*%(invSigma%*%du+dinvSigma%*%diag(u,N))+t(X)%*%dLam%*%invSigma%*%diag(u,N)
    S = -t(X)%*%W%*%X
    Maux = -t(Delta)%*%solve(S)%*%Delta
    resul = eigen(Maux, symmetric=T)$vectors[,1]
  }

  if(pert == "correlation"){
    R = bdiag(model$working.correlation)
    n = model$nclusters
    t = model$max.id
    dR = -R + diag(t)
    mu = model$fitted.values
    A = model$comp$A
    phi = model$scale
    daux10 = list(NULL)
    for(i in 1:n){
      daux10[[i]] = dR
    }
    dSl = sqrt(A)%*%bdiag(daux10)%*%sqrt(A)
    dSigma = bdiag(dSl)
    A = model$comp$A
    Rm = model$comp$Rm
    Sigma = (sqrtm(A)$B)%*%Rm%*%(sqrtm(A)$B)
    invSigma = solve(Sigma)
    dinvSigma <- -invSigma%*%dSigma%*%invSigma
    Delta = as.matrix(t(X)%*%(model$comp$Lambda)%*%dinvSigma%*%diag(u))
    S = -t(X)%*%W%*%X
    Maux = -t(Delta)%*%solve(S)%*%Delta
    resul = eigen(Maux, symmetric=T)$vectors[,1]
  }

  if(pert == "covariate"){
    if(is.null(covariate)){
      stop("Indicate the position of the covariate of interest")
    }
    l = covariate
    sl <-  sd(X[,l])
    eta = X%*%model$coefficients
    N = model$nobs
    dG <- diag(N)
    phi = model$scale
    mu = model$fitted.values
    Ww = diag((1+phi)*(psigamma(mu*(1+phi),3)-psigamma(mu*(1+phi)+phi+2, 3)))
    A = model$comp$A
    G = model$comp$G
    AA = solve(sqrt(A))%*%G%*%Ww
    p = length(model$coefficients)
    dX = matrix(0,p,N)
    dX[l,] <- sl
    beta = model$coefficients
    dLam = (1+phi)*beta[l]*sl*(Ww%*%G^2+dG^2%*%A)
    dSigma = 0.5*beta[l]*sl*(sqrt(A)%*%(model$comp$Rm)%*%AA+AA%*%(model$comp$Rm)%*%sqrt(A))
    dinvSigma = -solve(model$comp$Omega)%*%dSigma%*%solve(model$comp$Omega)
    du = -(1+phi)*beta[l]*sl*(G%*%A)
    Lambda = model$comp$Lambda
    Rm = model$comp$Rm
    Sigma = (sqrtm(A)$B)%*%Rm%*%(sqrtm(A)$B)
    invSigma = solve(Sigma)
    Delta = t(X)%*%Lambda%*%(invSigma%*%du+dinvSigma%*%diag(u))+(t(X)%*%dLam+dX%*%Lambda)%*%invSigma%*%diag(u)
    Maux <- -t(Delta)%*%solve(-t(X)%*%W%*%X)%*%Delta
    resul <- eigen(Maux, symmetric=T)$vectors[,1]
  }

  if(type == "normal"){
    return(resul)
  }else{
    dmax <- abs(resul)
    dn = dmax/sum(dmax)
    return(dn)
  }
}
