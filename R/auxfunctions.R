#' soft.thre
#'
#' @title softthresholding complex matrix G
#' @description Soft-thresholding of comlex matrix G with threshould phi.
#' @param G Input complex matrix
#' @param phi Threshold level parameter
#' @param diagTF If TRUE, apply thresholding on the diagonal entries of G
#' @keywords soft thresholding function
#' @export
#' @examples
#' soft.thre(G)
#'
#'

soft.thre = function(G, phi, diagTF=TRUE){
if(is.matrix(G)){ dim = ncol(G)
} else{ dim=1}

g = as.vector(G); g2 = numeric(length(g));
id1 = which(abs(g) > phi);
g2[id1] = (1 - phi/abs(g[id1]))*g[id1];
out = matrix(g2, ncol=dim);
if(!diagTF){
  diag(out) = diag(G);
}
return(out)
}


#' hard.thre() Function
#'
#' @title Hard thresholding complex matrix G
#' @description Hard-thresholding of comlex matrix G with threshould phi.
#' @param G Input complex matrix
#' @param phi Threshold level parameter
#' @param diagTF If TRUE, apply thresholding on the diagonal entries of G
#' @keywords hard.thre
#' @export
#' @examples
#' hard.thre(G)
#'
#'
#'

hard.thre = function(G, phi, diagTF=TRUE){
  if(is.matrix(G)){ dim = ncol(G)
  } else{ dim=1}


  g = as.vector(G); g2 = numeric(length(g));
  id1 = which(abs(g) > phi);
  g2[id1] = g[id1];
  out = matrix(g2, ncol=dim);
  if(!diagTF){
    diag(out) = diag(G);
  }
  return(out)
}


#' adaptivesoft.thre
#'
#' @title Adaptive soft-thresholding complex matrix G
#' @description Hard-thresholding of comlex matrix G with threshould phi.
#' @param G Input complex matrix
#' @param phi Phi parameter in the adaptive threshold function.
#' @param eta Eta parameter in the adaptive threshould function.
#' @param diagTF If TRUE, apply thresholding on the diagonal entries of G
#' @keywords adaptivesoft.thre
#' @export
#' @examples
#' adaptivesoft.thre(G)
#'
#'
#'
#'

adaptivesoft.thre = function(G, phi, eta=2, diagTF=TRUE){
  if(is.matrix(G)){ dim = ncol(G)
  } else{ dim=1}

  g = as.vector(G); g2 = numeric(length(g));
  tmp = abs(g) - (phi^(eta+1))/(abs(g))^eta
  id1 = which(tmp > 0);
  g2[id1] = g[id1]/abs(g[id1])*tmp[id1];
  out = matrix(g2, ncol=dim);
  if(!diagTF){
    diag(out) = diag(G);
  }
  return(out)
}


#' dist.Frobenius() Function
#'
#' @title Frobenious distance between two complex matrices
#' @description dist.Frobenius(G, H) calculates Frobenius distance between two complex matrices G and H.
#' @keywords Frobenious norm
#' @export
#' @examples
#' dist.Frobenius(G, H)
#'
#'
#'
#'

dist.Frobenius = function(G, H){
  aa = (G-H); bb = aa%*%t(Conj(aa));
  sqrt(Re(sum(diag(bb)))) ; }

#' dist.spectral() Function
#'
#' @title  Spectral distance between two complex matrices
#' @description dist.spectral(G, H) calculates Spectral distance between two complex matrices G and H.
#' @keywords Spectral norm
#' @export
#' @examples
#' dist.spectral(G, H)
#'
#'
#'
dist.spectral = function(G, H){
  aa = (G-H); bb = aa%*%t(Conj(aa));
  c = svd(bb)$d;
  return(sqrt(max(c)))
}

#' Det.complex() Function
#'
#' Determinant of complex square matrix
#' @param G Input complex matrix
#' @keywords Determinant
#' @export
#' @examples
#' Det.complex(G)
#'
#'
#'


Det.complex = function(G){
  #  prod(eigen(G, only.values=TRUE)$values) # Not stable for complex valued
  prod(svd(G)$d)
}


#' Periodogram Function
#'
#' @title Periodogram of multivariate time series
#' @param data Input data of dim*length
#' @keywords Periodogram
#' @export
#' @examples
#' Periodogram(data)
#'
#'
#'


Periodogram = function(data){
  Tt=ncol(data); p = nrow(data);
  N=floor(Tt/2);
  omega = 2*pi*(1:N)/Tt;

  ## Periodogram for each dimension
  I = array(0, c(p, p, N-1));
  for(i in 1:p){
    for(j in 1:p){
      dft1 = fft(data[i,]);
      dft1 = dft1[2:N];

      dft2 = fft(data[j,]);
      dft2 = dft2[2:N];
      I[i,j,]   = dft1*Conj(dft2)/(2*pi*Tt);
    }
  }
  return(I)
}


#' chol.complex
#' @title Cholesky decomposition of complex maxtrix G
#' @usage chol.complex(G)
#' @param G Input complex matrix
#' @keywords (complex) Cholesky decomposition
#' @export chol.complex
#' @examples
#' chol.complex(G)
#'
#'
#'

chol.complex= function(G){
  dim = nrow(G);
  L = matrix(0, dim, dim);

  intchol = function(A, L){
    p = ncol(A); p1 = p-1;
    A11 = A[1:p1, 1:p1];
    A12 = A[1:p1, p];
    A22 = A[p, p];
    ell = solve(L, A12); ell = Conj(ell);
    l2 = sqrt(A22 - Re(t(Conj(ell))%*%(ell)));
    return(c(ell, l2))
  }

  L[1,1] = sqrt(G[1,1]);
  L[2,1:2] = intchol(G[1:2, 1:2], L[1,1]);

  for(k in 3:dim){
    k1 = k-1;
    L[k, 1:k] = intchol(G[1:k, 1:k], L[1:k1, 1:k1])
  }

  #  G1 = L%*%(t(Conj(L)))
  #  if(sum(abs(G-G1)) > 10e-3){
  #    warning("Error")
  #  }
  return(L)
}


