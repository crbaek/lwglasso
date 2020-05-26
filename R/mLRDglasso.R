

#' MLRD.glasso() Function
#'
#' @title two-stage estimatin of LRD parameter and sparse G using glasso
#' @param data, matrix of dim*Tt
#' @param m number of frequencies used in estimation
#' @keywords Local Whittle estimation
#' @export
#' @examples
#' MLRD.glasso(data, m)
#'
#'
#'


MLRD.glasso = function(data, m, lambda = "ebic", type="soft", approxTF=TRUE, gridTF=TRUE, gg=1, bound=c(0.05, 1), debiasTF=FALSE){
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
  } else{
    ### First estimate D  ##############
    ret = lwMLRD(data, m);
    id = which(ret$dhat == 0 || ret$dhat == .5); ## If it hits boundary, use univariate one.
    ret$dhat[id] = ret$dinit[id];
    Ghat = ret$Ghat;
  }

  out = glasso.complex(Ghat, Tt=Tt, lambda= lambda, type=type, gridTF=gridTF, gg=gg, bound=bound, debiasTF=debiasTF);

  # Avagan et al. (2017) JCGS suggest to do glasso on half to reduce bias...
  #  Ghalf = Matpower(Ghat, 1/2);
  #  out1 = glasso.complex(Ghalf, Tt=Tt, lambda= lambda, gridTF=gridTF, gg=gg, bound=bound, debiasTF=debiasTF);
  #  xx = out1$iG%*%out1$iG

  fit = append(ret, out)
  fit$R = 1*( out$iG != 0);
  return(fit)
}


#' glasso.complex() Function
#'
#' @title Sparse estimation of inverse G using graphical lasso.
#' @param Ghat, Estimated Ghat to be sparsely estimated using glasso
#' @param m number of frequencies used in estimation
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
        Y = adaptivesoft.thre(X - mu*gradfX, mu*lambda);}
      if(type == "soft"){
        Y = soft.thre(X - mu*gradfX, mu*lambda);}
      if(type == "hard"){
        Y = hard.thre(X - mu*gradfX, mu*lambda);}
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


#' glassoD2() Function
#'
#' @title Sparse constrained estimation in complex Graphical model
#' @param S, input complex matrix
#' @param R, constraint matrix
#' @keywords constrained estimation of complex matrix
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



#' threGhat.bic() Function
#'
#' @title Sparse estimation using thresholding method
#' @param data, input data
#' @param m, number of frequencies used in LW estimation
#' @param lambda, Selection of lambda, ebic is default.
#' @keywords Thresholding LW long-run variance matrix
#' @export
#' @examples
#' threGhat.bic(data)
#'
#'
#'
#'
#'#####################################################
## Optimal lambda from BIC in the Thresholding method

threGhat.bic = function(data, m, lambda = "ebic", bound, approxTF=TRUE, adaptiveTF=TRUE, eta=1, gg=1, debiasTF=FALSE){
  Tt = dim(data)[2]; p = dim = dim(data)[1];

  if(missing(m)){ m = floor(Tt^.8); }

  #  ### If it takes too much to estimate D use univariate case
  if(approxTF){
    dinit = numeric(dim);
    for(i in 1:p){
      dinit[i] = lwe(data[i,], m, 0)[2];
    }
    ret = list(); ret$dhat = dinit;
    I = Periodogram(data);
    Ghat = LWGhat(I, dinit, m=m, Tt=Tt); ret$Ghat = Ghat;
  } else{
    ### First estimate D  ##############
    ret = lwMLRD(data, m);
    id = which(ret$dhat == 0 || ret$dhat == .5); ## If it hits boundary, use univariate one.
    ret$dhat[id] = ret$dinit[id];
    Ghat = ret$Ghat;
  }

  if(missing(bound)){
    ss = summary(Mod(as.vector(Ghat)));
    bound = c(max(ss[1], sqrt(log(dim)/Tt)), ss[6])
  }

  rhoset = seq(from=bound[1], to = bound[2], length=20);
  if(adaptiveTF){
    rr = lapply(rhoset, adaptivesoft.thre, G=Ghat, eta=eta, diagTF=FALSE)
  }else{
    rr = lapply(rhoset, soft.thre, G=Ghat, diagTF=FALSE)
  }

  bic = numeric(20);
  for(j in 1:20){
    Y = rr[[j]];
    if(debiasTF){
      Y = glassoD2(Ghat, R = 1*(Y!= 0) )$W; rr[[j]] = Y; }
    n2loglik = Tt*(log(Det.complex(Y)) + Re(sum(diag(Ghat%*%solve(Y)))));
    ell = sum( Y != 0) ;
    bic[j] = n2loglik + ell*(log(Tt) + 4*gg*log(dim));
  }

  id = which.min(bic);
  G = G.db = rr[[id]]
  bic.db = min(bic);

  out = list();
  out$rhoset = rhoset; out$bic = bic;
  out$lambda = rhoset[id];
  out$Ghat = Ghat; # Full Ghat estimation
  out$dhat  = ret$dhat;
  out$G0 = G; # only thresholding
  out$G = G.db; # reestimated to reduce bias
  out$R = 1*(G.db != 0);
  out$bic.db = bic.db;
  out$debiasTF = debiasTF;
  out$adaptiveTF = adaptiveTF;
  return(out)
}





