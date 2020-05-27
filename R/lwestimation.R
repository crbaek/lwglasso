
#' lwMLRD
#'
#' @title Local Whittle estimation of multivariate time series
#' @description This function provides local Whittle estimation of multivariate long range dependent time series
#'
#' @param data Input data (dim*length)
#' @param m The number of frequencies used in estimation.  If missed, $m=[T^.8]$.
#' @keywords Local Whittle estimation
#' @export
#' @examples
#' lwMLRD(data, m)
#'
#'
#'


lwMLRD = function(data, m){
  # Input data is dim*length

  Tt=ncol(data); p = nrow(data);
  N=floor(Tt/2);
  omega = 2*pi*(1:N)/Tt;

  I = Periodogram(data);

  ## Initial estimators
  dinit = numeric(p);
  for(i in 1:p){
    dinit[i] = lwe(data[i,], m, 0)[2];
  }

  mLWobj = function(d, I, Tt, m){

    Ghat = LWGhat(I, d, m, Tt);
    detG = Det.complex(Ghat);
    if(detG > 1e-6){
      lam =  2*pi*(1:m)/Tt;
      z1 =log(detG) - 2*sum(d)*(1/m)*sum(log(lam));
    } else{
      z1 = 10^20; }
    return(z1)
  }

  lowerd =rep(0, p); upperd =rep(.5, p);
  est = optim(par=dinit, fn=mLWobj, I=I, Tt=Tt, m=m, method = "L-BFGS-B", lower=lowerd, upper=upperd, control = list(maxit = 400, factr=1e2, lmm=4));

  # est =  nlminb(dinit, mLWobj, gradient = NULL, hessian = NULL,I=I, Tt=Tt, m=m,
  #         scale = 1, control = list(iter.max=600), lower = lowerd, upper = upperd)

  dhat = est$par;

  out = list();
  out$I = I; # Periodogram
  out$dinit=dinit; # initial value
  out$dhat = dhat;  # Estimated value
  out$Ghat = LWGhat(I, dhat, m, Tt);
  return(out)
}



#' LWGhat
#'
#' @title Long-run variance estimation based on local Whittle estimation of multivariate time series
#' @describeIn This function provides long-run variance of multivariate LRD series for given LRD parameter d and periodograms.
#' @param I Periodogram
#' @param d LRD parameters
#' @param m number of frequencies used in LW estimation of LRD parameters
#' @param Tt The number of sample size
#' @keywords Long-run variance estimation
#' @export
#' @examples
#' LWGhat(I, d, m, Tt)
#'
#'
#'
#'

LWGhat = function(I, d, m, Tt){
  p = length(d);
  G=matrix(0, p,p);
  ## Phase parameter estimate
  for(j in 1:m){
    ell = 2*pi*j/Tt;
    Psi = diag(ell^d);
    G=G+Psi%*%I[,,j]%*%Conj(Psi);
  }
  Ghat=G/m;
  return(Ghat)
}



#' lwe
#'
#' @title Local Whittle estimation of LRD parameter for univariate time series
#' @description This function provides the local Whittle estimation of LRD parameter.
#' @param data Input data (univariate time series)
#' @param m The number of frequencies used in estimation. If missing, then $T^{1/2}$ is used
#' @param lowerd the lower bound of LW object function, default = 0
#' @param upperd the upper bound of LW object function, default = 1
#' @keywords Local Whittle LRD parameter estimation
#' @export
#' @examples
#' lwe(data,m, lowerd=0, upperd=1)
#'
#'
#'
#'
#####################################################################
## Local whittle estimator for fixed frequencies #############
#####################################################################
lwe <-function(data,m, lowerd, upperd){

  if(missing(lowerd)){ lowerd = 0};
  if(missing(upperd)){ upperd = 1};


  T = length(data);
  #if(missing(m)){ m = T^(2/3)};
  if(missing(m)){ m = T^(1/2)};
  m = floor(m);


  ######################################################################
  # LW whittle estimator object function
  ######################################################################
  lw_R <- function(d, omega, power, m){
    lambda = omega[1:m];
    Ix = power[1:m];
    lwd = log(mean((lambda^(2*d))*Ix)) - 2*d*mean(log(lambda));
    return(lwd);
  }

  # First calculate periodogram
  M=floor(T/2);
  dft = fft(data);
  dft = dft[2:M];
  omega = 2*pi*seq(from=1, to =M, by=1)/T;
  power   = (abs(dft)^2)/(2*pi*T);
  y = log(power); y = y[1:M];

  nf = m;
  lwd = chat = 0*seq(1:length(nf));
  for(i in 1:length(lwd)){
    m1 = nf[i];
    y1 = y[1:m1];
    x = log(2*sin(omega[1:m1]/2))
    out = coef(lm(y1~1+x))
    GPH = -out[2]/2;
    init = GPH*(GPH > lowerd)*(GPH < upperd) + .1*(1-(GPH > lowerd)*(GPH < upperd));
    est = optim(par=init, lw_R, omega = omega, power=power, m= m1, method = "L-BFGS-B", lower=lowerd, upper=upperd);
    lwd[i] = est$par;
    chat[i] = mean(power[1:m]*(omega[1:m])^(2*lwd[i]));
  }
  if(length(m) < 2){
    out = c(m, lwd, chat);
    names(out) = c("nf", "LW-dhat", "chat"); }
  if(length(m) > 1 ){
    out = list();
    out$nf = m;
    out$lwd = lwd;
    out$chat = chat;}
  return(out)
}


#' lweplot
#'
#' @title Plot of local Whittle estimator over the number of frequencies used
#' @description lweplot depics local Whittle estimators over the number of frequencies used
#' @param data Input data (univariate time series)
#' @param m The maximum number of frequencies used in the plot. If missing, then $T^{.8}$ is used
#' @param lowerd the lower bound of LW object function, default = -.5
#' @param upperd the upper bound of LW object function, default = 1
#' @keywords Plot of local Whittle LRD parameter estimators
#' @export lweplot
#' @examples
#' lweplot(data)
#'
#'
#'
#'
lweplot<-function(data, m, lowerd, upperd){
  T = length(data);
  if(missing(m)){ m = T^.8 };
  if(missing(lowerd)){ lowerd = -.5};
  if(missing(upperd)){ upperd = 1};

  m = floor(m);

  ######################################################################
  # LW whittle estimator object function
  ######################################################################
  lw_R <- function(d, omega, power, m){
    lambda = omega[1:m];
    Ix = power[1:m];
    lwd = log(mean((lambda^(2*d))*Ix)) - 2*d*mean(log(lambda));
    return(lwd);
  }

  # First calculate periodogram
  M=floor(T/2);
  dft = fft(data);
  dft = dft[2:M];
  omega = 2*pi*(1:M)/T;
  power   = (abs(dft)^2)/(2*pi*T);
  nf = seq(10, m, 5);
  lwd = 0*seq(1:length(nf));
  for(i in 1:length(lwd)){
    est = optim(.1, lw_R, omega = omega, power=power, m= nf[i], method = "L-BFGS-B", lower=lowerd, upper=upperd);
    lwd[i] = est$par;
  }

  ub = lwd + qnorm(.975)/(2*sqrt(nf));
  lb = lwd - qnorm(.975)/(2*sqrt(nf));

  plot( nf, lwd, type="l", lty=1, xlab="Frequency", ylab="d", ylim=c(min(lb)-.1, max(ub)+.1), lwd=2);
  points(nf, ub, type="l", lty=2, col="red")
  points(nf, lb, type="l", lty=2, col="red")
  title('Local Whittle Estimator');

  out = list();
  out$nf = nf;
  out$lwd = lwd;

  return(out)
}




