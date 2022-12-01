

#' @name sparseMLRD
#'
#' @title two-stage estimation of LRD parameter and sparse inverse G using penalization methods
#' @description This function estimates sparse long-run variance (complex) matrix of multivariate LRD model
#' using local Whittle graphical lasso, spice (studentized version of graphical lasso) and CLIME.
#' @param data matrix of dim*Tt
#' @param m number of frequencies used in estimation.  If missed, m=[T^.8].
#' @param lambda "ebic" uses extended BIC criteria to select penalty parameter in graphical lasso.
#' User also can provide numerical value.
#' @param type Types of thresholding in ADMM algorithm. "hard", "soft" and "adaptive" threshold functions are possible.
#' @param method Three methods are possible. If "glasso", LW grapchial lasso is used. If "spice", long-run correlation matrix is used instead of covariance matrix.
#' If "clime" is used, LW CLIME estimator is provided. Default is "glasso".
#' @param approxTF If TRUE, univariate LRD parameter is used in the estimation. Otherwise, multivariate LRD parameter estimator is used.
#' @param gridTF If TRUE, penalty parameter is searched over interval provided on bound argument. Otherwise, optim function searches optimal
#' lambda.
#' @param gg The tuning parameter in the extended BIC criteria. Default value is 1. If gg=0, it is a usual BIC.
#' @param bound Bound of grid search in extended BIC. Default value is (.05, 1)
#' @param debiasTF If TRUE, debiased by applying constrained MLE introduced in the paper.
#' @keywords sparse Multivariate Local Whittle estimation
#' @export sparseMLRD
#' @examples
#' sparseMLRD(data, m)
#'
#'
#'

sparseMLRD = function(data, m, lambda = "ebic", type="soft", method="glasso", approxTF=TRUE, gridTF=TRUE, gg=1, bound=c(0.05, 1), debiasTF=TRUE){
  # data is dim*length matrix
  Tt = ncol(data);
  p = nrow(data);

  if(missing(m)){ m = floor(Tt^.8); }

  #  ### If it takes too much to estimate D use univariate case
  if(approxTF){
    dinit = numeric(p);
    for(i in 1:p){
      dinit[i] = lwe(data[i,], m, lowerd=0, upperd=1)[2];
    }
    ret = list(); ret$dhat = dinit;
    I = Periodogram(data);
    Ghat = LWGhat(I, dinit, m=m, Tt=Tt); ret$Ghat = Ghat;
    What = sqrt(1/Re(diag(Ghat))); What = diag(What);
    Gam = What%*%Ghat%*%What; ret$Gam = Gam;

  } else{
    ### First estimate D  ##############
    ret = lwMLRD(data, m);
    id = which(ret$dhat == 0 || ret$dhat == .5); ## If it hits boundary, use univariate one.
    ret$dhat[id] = ret$dinit[id];
    Ghat = ret$Ghat;
    What = sqrt(1/Re(diag(Ghat))); What = diag(What);
    Gam = What%*%Ghat%*%What; ret$Gam = Gam;
  }

  if(method == "spice"){
    out = glasso.complex(Gam, Tt=Tt, lambda= lambda, type=type, gridTF=gridTF, gg=gg, bound=bound, debiasTF=debiasTF);
#   out$iGam = out$iG;
    out$R = 1*(out$iG != 0);
    out$iG = What%*%out$iG%*%What;
#   out$iG = out$iG0 = What%*%out$iG%*%What;

    if(debiasTF){
      fit = glassoD2(S=Ghat, R = out$R);
      out$iG = fit$Theta;
    } }

  if(method == "glasso"){
    out = glasso.complex(Ghat, Tt=Tt, lambda= lambda, type=type, gridTF=gridTF, gg=gg, bound=bound, debiasTF=debiasTF);
  }

  if(method == "clime"){
    out = clime.complex(Ghat, Tt=Tt, lambda= lambda, type=type, gridTF=gridTF, gg=gg, bound=bound, debiasTF=debiasTF);
  }

  fit = append(ret, out)
  fit$R = 1*( out$iG != 0);
  return(fit)
}



#' glasso.complex
#'
#' @title Sparse estimation of inverse G using graphical lasso.
#' @description This function estimates sparse long-run variance (complex) matrix using graphical Lasso.
#' @param Ghat Estimated (nonsparse) long-run variance matrix to be sparsely estimated using glasso
#' @param Tt Data length
#' @param lambda "ebic" uses extended BIC criteria to select penalty parameter in graphical lasso.
#' User also can provide numerical value.
#' @param type Types of thresholding in ADMM algorithm. "hard", "soft" and "adaptive" threshold functions are possible.
#' @param approxTF If TRUE, univariate LRD parameter is used in the estimation. Otherwise, multivariate LRD parameter estimator is used.
#' @param gridTF If TRUE, penalty parameter is searched over interval provided on bound argument. Otherwise, optim function searches optimal
#' lambda.
#' @param gg The tuning parameter in the extended BIC criteria. Default value is 1. If gg=0, it is a usual BIC.
#' @param bound Bound of grid search in extended BIC. Default value is (.05, 1)
#' @param debiasTF If TRUE, debiased by applying constrained MLE introduced in the paper.
#' @keywords Local Whittle estimation
#' @export
#' @examples
#' glasso.complex(Ghat, Tt, lambda="ebic")
#'
#'
#'
#'

glasso.complex = function(Ghat, Tt, lambda="ebic", type="soft", gridTF=TRUE, gg=1, bound=c(0.05, 1), debiasTF=FALSE){

  ## Complex Glasso, also works for Real as well.
  glasso.admm = function(S, lambda=.2, gg=1, Tt=1, type, debiasTF=TRUE){

    ## here mu is 1/mu of the previous program
    dim = ncol(S)
    X = solve(S); Y= diag(diag(X));  Lam = 0*X;
    fc=10e20; tol.Xrel = tol.Yrel = tol.frel = 1e-9;
    maxK = 1000;  muf = 1e-3; # final mu
    rmu = 1/4; # ratio of decreasing mu

    mu=.01;
    st=0; iter = 1;
    ## Update X
    while(st < 1 && iter < maxK){
      Xp= X; Yp=Y;
      W = Y + mu*(Lam - S);
      #    W = (W + t(Conj(W)))/2; ## Make sure it is a symmetric matrix
      out = svd(W);
      gamma = (out$d + sqrt((out$d)^2 + 4*mu))/2;
      X = out$u%*%diag(gamma)%*%t(Conj(out$v));
      gradfX = S - out$u%*%diag(1/gamma)%*%t(Conj(out$v));

      ## Update Y
      if(type == "adaptive"){
        Y = adaptivesoft.thre(X - mu*gradfX, mu*lambda, diagTF=FALSE);}
      if(type == "soft"){
        Y = soft.thre(X - mu*gradfX, mu*lambda, diagTF=FALSE);}
      if(type == "hard"){
        Y = hard.thre(X - mu*gradfX, mu*lambda, diagTF=FALSE);}
      # Lam = -soft.thre(Y, lambda);
      Lam = gradfX + (X - Y)/mu; # This works for the Complex number

      ## Determine stopping rule
      fp = fc;
      fc = -sum(log(gamma)) + Re(sum(diag(S%*%X))) + lambda*sum(sum(abs(X)));
      frel = abs(fp - fc)/max(abs(fp), abs(fc), 1);

      Xrel = norm(abs(X-Xp), type="F")/max(1, norm(abs(X), "F"), norm(abs(Xp), "F"))
      Yrel = norm(abs(Y-Yp), type="F")/max(1, norm(abs(Y), "F"), norm(abs(Yp), "F"))
      st = ( frel <  tol.frel )*( Xrel <  tol.Xrel )*( Yrel < tol.Yrel );
      iter = iter+1;

      #  ## Another stopping rule is to check (Sigma - Lambda) is positive definite
      if( iter %% 10 == 0 ){
        #  chk = svd(S - Lam)
        #  gap = fc -sum(log(chk$d)) - dim;
        #  st = st +  (gap < .01);
        mu = max(mu*rmu, muf);
      }
    }

    R = 1*(Y != 0); iG = Y;
    ## Debiasing using structure.
    if(debiasTF){
      fit = glassoD2(S=S, R = R);
      iG = iG.db = fit$Theta;
    }
    num.para = sum(R) - dim(R)[1] + sum(diag(R) != 0);
    n2loglik = Tt*(Re(sum(diag(S%*%iG))) - log(Det.complex(iG)));
    bic = n2loglik  + num.para*(log(Tt) + 4*gg*log(dim));
    # bic = Tt*(log(Det.complex(fit$W)))  + num.para*(log(Tt) + 4*gg*log(dim));
    return(list(iG=iG, n2loglik = Re(n2loglik), bic=bic, iter=iter))
  }

  dim = nrow(Ghat);
  if( lambda == "ebic") {
    bic.rho = function(rho, Ghat, Tt, gg, type, debiasTF){
      out = glasso.admm(Ghat, lambda = rho, gg=gg, type=type, Tt=Tt, debiasTF=debiasTF);
      return(out$bic)
    }
    if(gridTF){
      ## Grid search?? Most stable..
      #    sqrho = seq(0.05, 1, by=.05);
      sqrho = seq(bound[1],  bound[2], length=20);
      est = lapply(sqrho, bic.rho, Ghat=Ghat, Tt=Tt, gg=gg, type=type, debiasTF=debiasTF);
      id = which.min(unlist(est)); rho = sqrho[id];
    } else{
      #   est = optim(par=.1, fn=bic.rho, Ghat=Ghat, Tt=Tt, gg=gg, debiasTF=debiasTF, method = "L-BFGS-B", lower=bound[1], upper=bound[2], control = list(maxit = 400, factr=1e2, lmm=4));
      est = optim(par=.1, fn=bic.rho, Ghat=Ghat, Tt=Tt, gg=gg, type=type, debiasTF=TRUE, method = "Brent", lower=bound[1], upper=bound[2], control = list(maxit = 200));
      rho = est$par;
      ## If it gives the boundary values / initial values, then do grid search again
      if(rho == bound[1] || rho == bound[2] || rho == .1){
        sqrho = seq(bound[1],  bound[2], length=20);
        est = lapply(sqrho, bic.rho, Ghat=Ghat, Tt=Tt, gg=gg, debiasTF=TRUE);
        id = which.min(unlist(est)); rho = sqrho[id];
      }
    }
  } else { rho = lambda;  } ## Fixed constant

  out = glasso.admm(Ghat, lambda = rho, Tt=Tt, type=type, gg=gg, debiasTF=debiasTF);

  result = out;
  result$lambda = rho;
  return(result)
}


#' glassoD2
#'
#' @title Constrained MLE for long-run variance matrix
#' @description This function provides constrained MLE for long-run variance matrix. It is used to debias the estimated
#' sparse long-run variance matrix using graphical Lasso.
#'
#' @param S input complex long-run variance matrix
#' @param R constraint matrix consisting of 1 (if nonzero) and 0 (if zero)
#' @param maxIt Maximum iteration of algorithm. Default is 100
#' @param tol Tolerance to stop iteratrion. Default is 1e-6
#'
#' @keywords Constrained estimation of long-run variance matrix.
#' @export
#' @examples
#' glassoD2(S, R)
#'
#'
#'


############################################
# Constrained problem. use this for debiasing
#############################################

glassoD2 = function(S, R, maxIt=100, tol = 1e-6){

  p = nrow(S);

  # Initialization
  W = S;   # diagonal of W remains unchanged
  W_old = W; Th = matrix(0, p, p);
  i =st = 0;

  # Graphical Lasso loop
  while( st < 1 && i < maxIt){
    i = i+1; sq = seq(from=p, to =1, by=-1);
    for( k in sq){
      j = sq[k]
      jminus = setdiff(1:p,j);
      W11 = W[jminus,jminus];
      s12 = S[jminus,j];    # W_11^(-1/2) * s_12
      r12 = R[jminus, j];
      b = solve(W11, s12); b1 = b*r12;
      w12  = W11 %*% b;
      w12 = w12*r12;
      ## Padding zero with constraints
      W[jminus,j] = w12;
      W[j,jminus] = t(Conj(W[jminus,j]));
      W = (W + t(Conj(W)))/2

      w22 = W[j,j];
      th22 = 1/(w22 - t(Conj(b1))%*%W[jminus, jminus]%*%b1);
      th22 = c(Re(th22));
      if(th22 < 0){ th22 = Re(1/W[j,j]); }
      Th[jminus, j] = -1*c(b1/th22);
      Th[j, jminus] =  Conj(Th[jminus, j]);
      Th[j,j] = th22;
    }
    # Stop criterion
    if( sum(Mod(W-W_old)) < tol ){ st = 1 }
    W_old = W;
  }

  # Theta = solve(W);
  return(list(Theta=solve(W), W=W))
}




#' clime.complex
#'
#' @title Sparse estimation of inverse G (precison matrix) using LW CLIME method.
#' @description This function estimates sparse long-run variance (complex) matrix using LW-CLIME method.
#' @param Ghat Estimated (nonsparse) long-run variance matrix to be sparsely estimated using LW-CLIME
#' @param Tt Data length, used in calculating BIC
#' @param lambda "ebic" uses extended BIC criteria to select penalty parameter in LW-CLIME.
#' User also can provide numerical value.
#' @param type Types of thresholding in ADMM algorithm. "hard", "soft" and "adaptive" threshold functions are possible.
#' @param approxTF If TRUE, univariate LRD parameter is used in the estimation. Otherwise, multivariate LRD parameter estimator is used.
#' @param gridTF If TRUE, penalty parameter is searched over interval provided on bound argument. Otherwise, optim function searches optimal
#' lambda.
#' @param gg The tuning parameter in the extended BIC criteria. Default value is 1. If gg=0, it is a usual BIC.
#' @param bound Bound of grid search in extended BIC. Default value is (.05, 1)
#' @param debiasTF If TRUE, debiased by applying constrained MLE introduced in the paper.
#' @keywords Local Whittle CLIME estimation
#' @export
#' @examples
#' clime.complex(Ghat, Tt, lambda="ebic")
#'
#'
#'
#'
#'
clime.complex = function(Ghat, Tt, lambda="ebic", type="soft", gridTF=TRUE, gg=1, bound=c(0.05, 1), debiasTF=FALSE){

  ## ADMM algorithm for CLIME estimator
  clime.admm = function(S, lambda=.2, gg=1, Tt=1000, type="soft", debiasTF=TRUE){

    dim = ncol(S);
    E = diag(dim);
    # Pertubation work well in the Real covariance
    X = solve(S + .1*E);
    #  X = solve(S);
    X= diag(diag(X));
    V = Y = 0*X; U = Z = E;

    tol.Vrel = tol.Yrel = tol.frel = 1e-5;
    maxK = 1000;  muf = 1e-3; # final mu
    rmu = 1/4; # ratio of decreasing mu

    mu=.01; rr=.2;
    st=0; iter = 1;
    ## Update X
    while(st < 1 && iter < maxK){

      Xp= X; Yp=Y; Vp=V; Up=U; Zp=Z;

      ## Update X
      if(type == "adaptive"){
        X = adaptivesoft.thre(Xp - Vp, mu, diagTF=FALSE);}
      if(type == "soft"){
        X = soft.thre(Xp - Vp, mu, diagTF=FALSE);}
      if(type == "hard"){
        X = hard.thre(Xp - Vp, mu, diagTF=FALSE);}

      U1 = S%*%X + Yp;
      Z = U1 + soft.thre(E - U1, lambda);
      Y = U1 - Z;
      V = S%*%(2*Y - Yp)*mu*rr;
      ## Make V as Hermitian
      V = (V + t(Conj(V)))/2

      ## Determine stopping rule
      Vrel = norm(abs(V-Vp), type="F")/max(1, norm(abs(V), "F"), norm(abs(Vp), "F"))
      #    Yrel = norm(abs(Y-Yp), type="F")/max(1, norm(abs(Y), "F"), norm(abs(Yp), "F"))
      #    st = 1*( Yrel < tol.Yrel )*(Vrel < tol.Vrel);
      #    print(Vrel)
      st = 1*(Vrel < tol.Vrel);
      iter= iter+1;
      mu = max(mu*rmu, muf);
      #    if( iter %% 10 == 0 ){
      #        mu = max(mu*rmu, muf);
      #    }

    }

    R = 1*(X != 0); iG = X;
    ## Debiasing using structure.
    if(debiasTF){
      fit = glassoD2(S=S, R = R);
      iG = iG.db = hard.thre(fit$Theta, 1e-6);
    } else{
      iG = hard.thre((X + t(Conj(X)))/2, 1e-6);
    }
    R = 1*(iG != 0);
    num.para = sum(R) - dim(R)[1] + sum(diag(R) != 0);
    n2loglik = Tt*(Re(sum(diag(S%*%iG))) - log(Det.complex(iG)));
    bic = n2loglik  + num.para*(log(Tt) + 4*gg*log(dim));
    return(list(iG=iG, n2loglik = Re(n2loglik), bic=bic, iter=iter, R=R))
  }


  dim = nrow(Ghat);
  if( lambda == "ebic") {
    bic.rho = function(rho, Ghat, Tt, gg, type, debiasTF){
      out = clime.admm(Ghat, lambda = rho, gg=gg, type=type, Tt=Tt, debiasTF=FALSE);
      return(out$bic)
    }
    if(gridTF){
      ## Grid search?? Most stable..
      sqrho = seq(bound[1],  bound[2], length=50);
      est = lapply(sqrho, bic.rho, Ghat=Ghat, Tt=Tt, gg=gg, type=type, debiasTF=debiasTF);
      id = which.min(unlist(est)); rho = sqrho[id];
    } else{
      #   est = optim(par=.1, fn=bic.rho, Ghat=Ghat, Tt=Tt, gg=gg, debiasTF=debiasTF, method = "L-BFGS-B", lower=bound[1], upper=bound[2], control = list(maxit = 400, factr=1e2, lmm=4));
      est = optim(par=.1, fn=bic.rho, Ghat=Ghat, Tt=Tt, gg=gg, type=type, debiasTF=TRUE, method = "Brent", lower=bound[1], upper=bound[2], control = list(maxit = 200));
      rho = est$par;
      ## If it gives the boundary values / initial values, then do grid search again
      if(rho == bound[1] || rho == bound[2] || rho == .1){
        sqrho = seq(bound[1],  bound[2], length=20);
        est = lapply(sqrho, bic.rho, Ghat=Ghat, Tt=Tt, gg=gg, debiasTF=TRUE);
        id = which.min(unlist(est)); rho = sqrho[id];
      }
    }
  } else { rho = lambda;  } ## Fixed constant

  out = clime.admm(Ghat, lambda = rho, Tt=Tt, type=type, gg=gg, debiasTF=debiasTF);

  result = out;
  result$lambda = rho;
#  result$est = unlist(est);
  return(result)
}



#' threGhat.cv
#'
#'
#' @title Sparse estimation using thresholding method
#' @description This function sparsely estimates long-run variance matrix of multivariate LRD series using threshoulding method.
#' @param data input data
#' @param m number of frequencies used in estimation. If missed, $m=T^.8$.
#' @param bound Bound of threshold values in grid search. Default value is (.01, 1)
#' @param eta The parameter in the adaptive threshould function. Default is 1.
#' @param gg The tuning parameter in calculating extended BIC. Default value is 1. If gg=0, it is a usual BIC.
#' @param debiasTF If TRUE, debiased by applying constrained MLE introduced in the paper.
#' @keywords Thresholding LW long-run variance matrix
#' @export
#' @examples
#' threGhat.cv(data)
#'
#'
#'
#'
# Randomly sample periodogram to calculate Frobenius norm
threGhat.cv = function(data, m, bound=c(0.01, 1), eta=1, gg=1, debiasTF=TRUE){
  Tt = dim(data)[2]; p = dim = dim(data)[1];

  dinit = numeric(dim);
  for(i in 1:p){
    dinit[i] = lwe(data[i,], m, 0)[2];
  }
  I1 = Periodogram(data);
  ret = lwMLRD(data, m);
  if(missing(bound)){
    ss = summary(Mod(as.vector(ret$Ghat)));
    bound = c(max(ss[1], sqrt(log(dim)/Tt)), ss[6])
  }

  LWGhat.index = function(I, d, index, Tt){
    p = length(d);
    G=matrix(0, p,p);
    len = length(index)
    ## Phase parameter estimate
    for(j in 1:len){
      k = index[j];
      ell = 2*pi*k/Tt;
      Psi = diag(ell^d);
      G=G+Psi%*%I[,,k]%*%Conj(Psi);
    }
    Ghat=G/len;
    return(Ghat)
  }

  len=50;
  rhoset = seq(from=bound[1], to = bound[2], length=len);
  E1 = matrix(0, 30, len)
  for(r in 1:30){
    index = sample(1:m, floor(m/2), replace=FALSE);
    Ghat1 = LWGhat.index(I=I1, d=dinit, index = index, Tt=Tt)

    index2 = setdiff(1:m, index)
    Ghat2 = LWGhat.index(I=I1, d=dinit, index = index2, Tt=Tt)


    out1 = lapply(rhoset, adaptivesoft.thre, G=Ghat1, eta=eta);
    err1 = lapply(out1, dist.Frobenius, H=Ghat2);
    E1[r,] = simplify2array(err1);
  }
  E = colMeans(E1);
  id = which.min(E);

  G = adaptivesoft.thre(ret$Ghat, rhoset[id], eta=eta);
  ff = list();
  R = 1*(G != 0);
  ff$R = R;

  if(debiasTF){
    G = glassoD2(ret$Ghat, R = R)$W;
    num.para = sum(R) - dim(R)[1] + sum(diag(R) != 0);
    n2loglik = Tt*(Re(sum(diag(ret$Ghat%*%solve(G)))) + log(Det.complex(G)));
    bic.db = n2loglik  + num.para*(log(Tt) + 4*gg*log(dim));
    ff$bic.db = bic.db;
  }

  ff$Ghat = ret$Ghat;
  ff$G = G;
  ff$bound = bound;
  ff$rhoset = rhoset;
  ff$dhat = dinit;
  ff$err = E;
  ff$lambda = rhoset[id];
  ff$debiasTF= debiasTF;
  return(ff)
}




