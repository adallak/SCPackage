
#' Computes smoothed estimate for a fixed tuning parameter with different penalty options.
#'
#' Solves the optimization problem in Dallakyan and Pourahmadi (2019) using
#' block-coordinate algorithm.
#'
#' @param S        Sample Covariance Matrix
#' @param lambda1  lambda value to control the sparsity
#' @param lambda2  lambda value to control smoothness
#' @param max_iter Number of maximum iteration
#' @param init.x   initial value for vectorized Cholesky factor L.
#' @param band     Positive number to select the band. If specified, algorithm forces the entries outside of band equal to zero and iterates only over specifed subdiagonals inside the band.
#' @param ABSTOL   Tolerance for algorithm convergence.
#'
#' @return Returns the estimated smoothed Cholesky factor \code{L}
#' @export 
#' 
#' smoothchol <-
#' function(S, lambda1 = 0, lambda2 , max_iter = 50, init.x = NULL,
#' band = NULL, ABSTOL = 1e-3)
#'
#' @examples 
#' set.seed(12)
#' p <- 50
#' band <- 5
#' L_true <- generateL(p = p, band = band)$L
#' library(varband)
#' X <- sample_gen(L = L_true, n = n)
#' sample covariance matrix
#' S <- crossprod(scale(X, center = TRUE, scale = FALSE)) / n
#' L_fused = smoothchol(S, lambda1 = 0, lambda2 = 0.2, type = "fused")$L
#' L_trend = smoothchol(S, lambda1 = 0, lambda2 = 0.2, type = "l1trend")$L
#' L_HP = smoothchol(S, lambda1 = 0, lambda2 = 0.2, type = "HP")$L
#' 
#' #' @seealso \code{\link{smoothcholCV}}
#' 
smoothchol<-function(S, lambda1 = 0, lambda2, max_iter = 70, init.x = NULL, type= c("fused", "l1trend", "HP"), 
                     band=NULL , ABSTOL   = 1e-3 )
{
  type = match.arg(type)
  lambda2 <- lambda2
  lambda1 <- lambda1
  max_iter <- max_iter
  p <- nrow(S)
  mat <- matGenerate(p)
  A <- mat$A
  ND <- mat$ND
  if(is.null(band))
  {
    band = p
  }
  
  if(is.null(init.x))
  {
    x <- matrix(0, p * (p+1) / 2, 1);
    x[A[[1]]] = 1 / sqrt(diag(S))
  }else{
    x <- init.x
  }
  if (type == "fused"){
    fused = iter_fused(x, S, lambda1, lambda2, band, A, max_iter, ABSTOL )
    x = fused$x
    history = fused$history
  }
  if(type == "l1trend"){
    l1trend = iter_trend(x, S, lambda1, lambda2, band,A, max_iter, ABSTOL )
    x = l1trend$x
    history = l1trend$history
  }
  if(type == "HP")
  {
    HP = iter_hp(x, S, lambda1, lambda2, band,A, max_iter, ABSTOL )
    x = HP$x
    history = HP$history
  }
  
  ind <- append(A[[1]] , ND)
  x <- x[ind]
  L <- Lfromx(x, p)
  return = list(x = x, history = history, L = L)
  return(return)
  
}

##################################################################
### This function converts vector x to lower triangular matrix L

Lfromx <- function(x, p)
{
  L <- matrix(0, p, p)
  diag(L) = x[1 : p]
  L[p,1] = x[ length(x) ]
  k = 0
  m <- 0
  for (i in 1 : (p - 2))
  {
    l = length(diag(L[-( 1 : i), -(p :(p - i + 1))]))
    k = k + l
    #    print(k)
    #    cat("L length",length(diag(L[-(1:i),-(p:(p-i+1))])),"\n")
    diag(L[-( 1 : i), - (p : (p - i + 1))]) = x[(p + m + 1) : (p + k)]
    m <- k
  }
  return(L)
}

###################################################
matGenerate<-function(p)
{
  ## This funciton generates index set for each subdiagonal
  p <- p
  A <- array(list(),dim=c(p,1))
  ND<-c()
  start = 1
  end = p
  for (i in 1 : (p))
  {
    ind <- p-i
    A[i]<-list(start : end)
    start = start + ind + 1
    end   = end + ind
    if( i > 1)
    {
      ND <- append(ND , A[[i]])
    }
  }
  output=list(A=A,ND=ND)
  return(output)
}
