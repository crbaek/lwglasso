
#' lwMLRD() Function
#'
#' Local Whittle estimation of multivariate time series
#' @param data, m
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



#' LWGhat() Function
#'
#' Long-run variance estimation based on local Whittle estimation of multivariate time series
#' @param I, d, m, Tt
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



#' lwe() Function
#'
#' Local Whittle estimation of LRD parameter for univariate time series
#' @param data,m, lowerd, upperd
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




