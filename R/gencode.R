

#' GenMultivarGauss () Function
#'
#' @title Generating Multivariate Gaussian process using circular embedding
#' @param R, given autocovariance matrix
#' @param N, sample size
#' @keywords Gaussian Process generation
#' @export
#' @examples
#' GenMultivarGauss(R,N, isAsymm=TRUE)
#'
#'
#'

###############################
# Generation Code
###############################
GenMultivarGauss = function(R,N, isAsymm=TRUE){

  ### Internal Fucntios
  circembed = function(Rtemp, isAsymm=TRUE){

    nd = dim(Rtemp)[1];   N = dim(Rtemp)[3];
    Rcirc = array(0, c(nd, nd, 2*N-2))
    for(k in 1:nd){
      Rcirc[k,k,] = c(Rtemp[k,k,], Rtemp[k,k,seq(from=N-1, to=2, by=-1)]);
    }

    for(k in 1:(nd-1)){
      for(j in (k+1):nd){
        if(isAsymm){
          Rcirc[k,j,] = 2*c( Rtemp[k,j,-N], Rtemp[j,k, seq(from=N, to=2, by=-1)]);
        } else{
          Rcirc[k,j,] = 2*c( Rtemp[k,j,-N], Rtemp[k,j, seq(from=N, to=2, by=-1)]);
        }
      }
    }

    return(Rcirc)
  }


  InitMultivarGauss = function(Rcirc){

    P = dim(Rcirc)[1];
    M = dim(Rcirc)[3];

    Lambda = 0*Rcirc
    for( k in 1:P){
      for(j in 1:P){
        Lambda[k,j,] = fft(Rcirc[k,j,]);
      }
    }

    sqrtLambda = matrix(0, P, M);
    for(k in 1:P){
      sqrtLambda[k,] = sqrt(Re(Lambda[k,k,])); }

    crossSpec = 0*Rcirc;
    for(k in 1:(P-1)){
      for(j in (k+1):P){
        Lxy = Lambda[k,j,];
        Sxy = numeric(M);
        tmp = sqrtLambda[k,]*sqrtLambda[j,];
        nonzeroInd = which(tmp !=0);
        Sxy[nonzeroInd] = Lxy[nonzeroInd]/tmp[nonzeroInd];
        crossSpec[k,j,] = Sxy;
      }
    }
    return(list(sqrtLambda=sqrtLambda, crossSpec=crossSpec))
  }

  SynthStepMultivarGauss = function(sqrtLambda,crossSpec,N){

    P = nrow(sqrtLambda);
    M = ncol(sqrtLambda);
    Z = matrix(0, P, M);
    for(m in 1:M){
      Cupper = crossSpec[,,m];
      Cm = Cupper + t(Cupper);
      diag(Cm) = 2;
#      Am = t(chol(Re(Cm)));
      aa = svd(Cm);
      Am = aa$u%*%diag(sqrt(aa$d))%*%t(Conj(aa$v)) ## Use SVD to get half
      W = complex(real=rnorm(P), imaginary=rnorm(P))/sqrt(2);
      Z[,m] = Am%*%W;
    }

    tmp = sqrtLambda*Z;
    Xtilde = 1/sqrt(M)*t(apply(tmp, 1, fft)); # Note fft should be applied to rowwise, apply gives column matrix
    X = Re(Xtilde[,1:N]);
#   Xiid = Im(Xtilde[,1:N]);
    return(X)
  }

  ################################################
  Rcirc = circembed(R, isAsymm);
  init = InitMultivarGauss(Rcirc);
  X = SynthStepMultivarGauss(init$sqrtLambda, init$crossSpec, N);
  return(X)
}



#' CovarFARIMA () Function
#'
#' @title Covariance function in mlw paper (21)
#' @param N, sample size
#' @param d, LRD parameter
#' @param delta delta parameter
#' @keywords Gaussian Process generation
#' @export
#' @examples
#' CovarFARIMA(N,d,delta,Sigmae)
#'
#'
#'
#'


CovarFARIMA = function(N,d,delta,Sigmae){

# ensure the elements in d satisfy -1/2<d(k)<1/2

if( sum((d<=-1/2 | d>=1/2)) > 0 ){
  message("d does not satisfy -1/2<d(k)<1/2, for all k")
};

if( sum((delta<=-1/2 | delta>=1/2)) > 0 ){
  message("delta does not satisfy -1/2<delta(k)<1/2, for all k")
};

p = length(d);

# Calculate R
# Variable storing covariance, RR(:,:,n+1)=EX_n X_0*

Gamma0.R = function(n, d, delta, Sigmae){
  p = length(d);
  Gamma0 = matrix(0, p, p);

  for(j in 1:p){
    for(k in 1:p){

  if(j == k ){
  num2 = (1/pi)*gamma(1-d[j]-d[j])*sin(d[j]*pi);
  Gamma0[j,j] = Sigmae[j,j]*num2*exp(lgamma(n+d[j])-lgamma(n+1-d[j]));
  } else{
  num1 = (1/pi)*gamma(1-delta[j]-delta[k])*sin(delta[k]*pi);
  num3 = (1/pi)*gamma(1-delta[k]-delta[j])*sin(delta[j]*pi);
  Gamma0[j,k] = Sigmae[j,k]*num1*exp(lgamma(n+delta[k])-lgamma(n+1-delta[j]));
  Gamma0[k,j] = Sigmae[k,j]*num3*exp(lgamma(n+delta[j])-lgamma(n+1-delta[k]));
  }

    }
  }
return(Gamma0)
}

  RR  = lapply(0:(N-1), Gamma0.R, d=d, delta=delta, Sigmae=Sigmae)
  RRneg = lapply(RR, t);
## Return object is a list
return(list(R=simplify2array(RR), RRneg=simplify2array(RRneg)))
}



#' CovtwoVARFIMA() Function
#'
#' @title two-sided VARFIMA(0, D, 0) model in SPIE paper equation (18)
#' @param N, sample size
#' @param d, LRD parameter, ensure the elements in d satisfy -1/2<d(k)<1/2
#' @param R, R is a real upper-triangular matrix representing sparsity
#' @param X, X is also complex valued upper-triangular matrix such that invG = (RW)^H%*%(RH)
#' @keywords Gaussian Process generation
#' @export
#' @examples
#' CovtwoVARFIMA(N, d, R, X)

CovtwoVARFIMA = function(N, d, R, X){
# R is a real upper-triangular matrix representing sparsity
# X is also complex valued upper-triangular matrix such that invG = (RW)^H%*%(RH)
# d is a LRD parameters
  # ensure the elements in d satisfy -1/2<d(k)<1/2

  if( sum((d<=-1/2 | d>=1/2)) > 0 ){
    message("d does not satisfy -1/2<d(k)<1/2, for all k")
  };

  p = length(d);
  rX = R*X; # Rademacher product
  iG = t(Conj(rX))%*%rX;
  Z = solve(rX);
  eD = diag(exp(1i*(pi/2)*d));
  eDm = diag(exp(-1i*(pi/2)*d));
  Z = Z%*%eDm;
  G = Z%*%(t(Conj(Z)))

  Qp = Qn = matrix(0, p, p);
  diag(Qp) = sqrt(2*pi)*eD%*%diag(Z);

  for(j in 2:p){
    for(k in 1:j){
      c1 = Re(Z[j,k])/cos(d[j]*pi/2);
      s1 = Im(Z[j,k])/sin(d[j]*pi/2);
      Qp[j,k] = sqrt(2*pi)*(c1-s1)/2;
      Qn[j,k] = sqrt(2*pi)*(c1+s1)/2;
    }
  }
  Qp = Re(Qp); Qn = Re(Qn);
  ## Check whether this gives the same G, equation (19) of SPIE paper
  #W = eDm%*%Qp + eD%*%Qn
  #W = W/sqrt(2*pi);
  #G1 = W%*%(t(Conj(W)))
  # Seems to be correct

  ## Apply Prposition 5.1 to get the covariance function
  b1 = Qn%*%(t(Qn));
  b2 = Qn%*%(t(Qp));
  b3 = Qp%*%(t(Qp));
  b4 = Qp%*%(t(Qn));

  Gamma0.R = function(n, b1, b2, b3, b4){

  Gamma0 = matrix(0, p, p);
  if(n > 0){
  for(j in 1:p){
    for(k in 1:p){
      g1 = 2*gamma(1-d[j]-d[k])*sin(d[k]*pi)*exp(lgamma(n+d[k])-lgamma(n+1-d[j]));
      g3 =  2*gamma(1-d[k]-d[j])*sin(d[j]*pi)*exp(lgamma(n+d[j])-lgamma(n+1-d[k]));
      g4 = 2*pi*exp(lgamma(n+d[j]+d[k])-lgamma(d[j]+d[k]) -lgamma(1+n));
      Gamma0[j,k] = g1*b1[j,k]+g3*b3[j,k]+g4*b4[j,k];
    }}
  } else{
    for(j in 1:p){
      for(k in 1:p){
        g1 = 2*gamma(1-d[j]-d[k])*sin(d[k]*pi)*exp(lgamma(n+d[k])-lgamma(n+1-d[j]));
        g3 =  2*gamma(1-d[k]-d[j])*sin(d[j]*pi)*exp(lgamma(n+d[j])-lgamma(n+1-d[k]));
        g4 = 2*pi*exp(lgamma(n+d[j]+d[k])-lgamma(d[j]+d[k]) -lgamma(1+n));
        Gamma0[j,k] = g1*b1[j,k]+g3*b2[j,k]+g3*b3[j,k]+g4*b4[j,k];
      }}
  }
  return(Gamma0/(2*pi))
  }

  RR  = lapply(0:(N-1), Gamma0.R, b1=b1, b2=b2, b3=b3, b4=b4);
  RRneg = lapply(RR, t);
  ## Return object is a list
  return(list(R=simplify2array(RR), RRneg=simplify2array(RRneg), G=G, iG=iG))
}

#' CovtwoVARFIMA2() Function
#'
#' @title two-sided VARFIMA(0, D, 0) model in SPIE paper equation (18)
#' @param N, sample size
#' @param d, LRD parameter, ensure the elements in d satisfy -1/2<d(k)<1/2
#' @param G, G is long-run variance matrix such that G = Z*t(Con(Z));
#' @keywords Gaussian Process generation
#' @export
#' @examples
#' CovtwoVARFIMA2(N, d, G)

CovtwoVARFIMA2= function(N, d, G){
  # G is long-run variance matrix such that G = Z*t(Con(Z));
  # Cholesky decomposition of complex-valued matrix G
  # d is a LRD parameters
  # ensure the elements in d satisfy -1/2<d(k)<1/2

  if( sum((d<=-1/2 | d>=1/2)) > 0 ){
    message("d does not satisfy -1/2<d(k)<1/2, for all k")
  };

  p = length(d);
  eD = diag(exp(1i*(pi/2)*d));
  eDm = diag(exp(-1i*(pi/2)*d));
  Qp = Qn = matrix(0, p, p);
  Z = chol.complex(G);
  Z = Z%*%eDm;
  G = Z%*%(t(Conj(Z)))
  iG = solve(G);

  Qp = Qn = matrix(0, p, p);
  diag(Qp) = sqrt(2*pi)*eD%*%diag(Z);

  for(j in 2:p){
    for(k in 1:j){
      c1 = Re(Z[j,k])/cos(d[j]*pi/2);
      s1 = Im(Z[j,k])/sin(d[j]*pi/2);
      Qp[j,k] = sqrt(2*pi)*(c1-s1)/2;
      Qn[j,k] = sqrt(2*pi)*(c1+s1)/2;
    }
  }
  Qp = Re(Qp); Qn = Re(Qn);
  ## Check whether this gives the same G, equation (19) of SPIE paper
  #W = eDm%*%Qp + eD%*%Qn
  #W = W/sqrt(2*pi);
  #G1 = W%*%(t(Conj(W)))
  # Seems to be correct

  ## Apply Prposition 5.1 to get the covariance function
  b1 = Qn%*%(t(Qn));
  b2 = Qn%*%(t(Qp));
  b3 = Qp%*%(t(Qp));
  b4 = Qp%*%(t(Qn));

  Gamma0.R = function(n, b1, b2, b3, b4){

    Gamma0 = matrix(0, p, p);
    if(n > 0){
      for(j in 1:p){
        for(k in 1:p){
          g1 = 2*gamma(1-d[j]-d[k])*sin(d[k]*pi)*exp(lgamma(n+d[k])-lgamma(n+1-d[j]));
          g3 =  2*gamma(1-d[k]-d[j])*sin(d[j]*pi)*exp(lgamma(n+d[j])-lgamma(n+1-d[k]));
          g4 = 2*pi*exp(lgamma(n+d[j]+d[k])-lgamma(d[j]+d[k]) -lgamma(1+n));
          Gamma0[j,k] = g1*b1[j,k]+g3*b3[j,k]+g4*b4[j,k];
        }}
    } else{
      for(j in 1:p){
        for(k in 1:p){
          g1 = 2*gamma(1-d[j]-d[k])*sin(d[k]*pi)*exp(lgamma(n+d[k])-lgamma(n+1-d[j]));
          g3 =  2*gamma(1-d[k]-d[j])*sin(d[j]*pi)*exp(lgamma(n+d[j])-lgamma(n+1-d[k]));
          g4 = 2*pi*exp(lgamma(n+d[j]+d[k])-lgamma(d[j]+d[k]) -lgamma(1+n));
          Gamma0[j,k] = g1*b1[j,k]+g3*b2[j,k]+g3*b3[j,k]+g4*b4[j,k];
        }}
    }
    return(Gamma0/(2*pi))
  }

  RR  = lapply(0:(N-1), Gamma0.R, b1=b1, b2=b2, b3=b3, b4=b4);
  RRneg = lapply(RR, t);
  ## Return object is a list
  return(list(R=simplify2array(RR), RRneg=simplify2array(RRneg), G=G, iG=iG))
}


#' CovoneVARFIMA() Function
#'
#' @title One-sided VARFIMA(0, D, 0) model in SPIE paper equation (18)
#' @param N, sample size
#' @param d, LRD parameter, ensure the elements in d satisfy -1/2<d(k)<1/2
#' @param invS is a sparse inverse matrix of innovations Sigma matrix
#' @keywords Gaussian Process generation
#' @export
#' @examples
#' CovoneVARFIMA(N, d, invS)
#'

CovoneVARFIMA = function(N, d, invS){
  # invS is a sparse inverse matrix of innovations Sigma matrix
  # d is a LRD parameters
  # ensure the elements in d satisfy -1/2<d(k)<1/2

  if( sum((d<=-1/2 | d>=1/2)) > 0 ){
    message("d does not satisfy -1/2<d(k)<1/2, for all k")
  };

  p = length(d);
  S = solve(invS);
  Qp = t(chol(S));

  ei = exp(-1i*pi/2*d)
  trueG = diag(ei)%*%S%*%diag(Conj(ei))/(2*pi);
  trueiG = solve(trueG);
  b3 = S;

  Gamma0.R = function(n, b3){

    Gamma0 = matrix(0, p, p);
      for(j in 1:p){
        for(k in 1:p){
          g3 =  2*gamma(1-d[k]-d[j])*sin(d[j]*pi)*exp(lgamma(n+d[j])-lgamma(n+1-d[k]));
          Gamma0[j,k] = g3*b3[j,k];
        }}
    return(Gamma0/(2*pi))
  }

  RR  = lapply(0:(N-1), Gamma0.R, b3=b3);
  RRneg = lapply(RR, t);
  ## Return object is a list
  return(list(R=simplify2array(RR), RRneg=simplify2array(RRneg), G = trueG, iG=trueiG))
}

