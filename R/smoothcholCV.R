sc_seq<-function(X, lambda_seq , init.x = NULL, 
                 lambda.type = c("lambda1", "lambda2"), stand = FALSE,
                 lambda2 = 0, lambda1 = 0, max_iter=50, 
                 pen.type=c("HP","fused","l1trend"), band, ABSTOL = 1e-3,
                 Stest)
{
  n <- dim(X)[1]
  p <- dim(X)[2]
  n_lambda = length(lambda_seq)
  penalty <- match.arg(pen.type)
  lambda.type  =  match.arg(lambda.type)
  tmp = c()
  ### Standardize the input
  if (isTRUE(stand))
  {
    S =  crossprod(scale(X, center = TRUE, scale = TRUE)) / (n - 1)
  }else{
    S = cov(X)
  }
  x_mat = matrix(0, nrow = p * (p+1) / 2, ncol = n_lambda ) 
  for (i in 1 : n_lambda){
    if (lambda.type == "lambda2")
    {
      x_mat[, i] <- smoothchol(S, lambda1 = lambda1, lambda2 = lambda_seq[i],
                               max_iter=max_iter, init.x = init.x, band=band, 
                               type=penalty, ABSTOL = ABSTOL)$x
    }else{
      x_mat[, i] <- smoothchol(S, lambda1 = lambda_seq[i], lambda2 = lambda2,
                               max_iter=max_iter, init.x = init.x, band=band, 
                               type=penalty, ABSTOL = ABSTOL)$x
    }
    init.x = x_mat[, i]   ## Updates the vector x
    omega = crossprod(Lfromx(x_mat[, i], p))
    tmp[i] = likelihood(omega, Stest)
  }
  
  
  return(list(lambda_seq = lambda_seq, x_mat = x_mat, tmp = tmp))
}



#' Performs kfolds-cross validation
#'
#'Selects tuning parameters by cross validation according to the likelihood on testing data.
#'
#' @param k  Folds used in cross-validation. The default is $k = 5$
#' @param X   A n-by-p sample matrix, each row is an observation of th p-dim random vector.
#' @param both.lambda Logical. If TRUE the cross-validation implemented for both lambdas.
#' @param lambda1_seq A vector of non-negative tuning parameters for lambda1 to control sparsity.
#' @param lambda2_seq A vector of non-negative tuning parameters for lambda1 to control smoothness
#' @param max_iter    Maximum number of iterations
#' @param band        Positive number of subdiagonal to be estimated. If not provided, the algorithm iterates over all subdiagonals.
#' @param n_lambda    If lambda1_seq and lambda2_seq is not provided, create a vector of lambdas with length n_lambda. Default is 60.
#' @param pen.type    Selects penalty for smoothness. 
#' @param ABSTOL      The tolerence for convergence
#' @param stand       Logical, if TRUE the data will be standardized.
#'
#' @return
#' A list object containing
#' \itemize{
#'  \item{lambda1_min: }{Selected value of lambda1 based on cross validation.}
#'  \item{lambda2_min: }{Selected value of lambda1 based on cross validation.}
#'  \item{L_fit: }{Estimate of L corresponding to the best fit.}
#'  \item{lambda1_seq: }{lambda1 grid used in cross validation.}
#'  \item{lambda2_seq: }{lambda2 grid used in cross validation.}
#' }
#' 
#' @export
#' 
#'
#' @examples
#' set.seed(11)
#' require(varband)
#' n = 100
#' p = 50
#' L_true = generateL(p = p, case = "c")$L
#' X = sample_gen(L = L_true, n = n)
#' L_cv = smoothcholCV(k = 5,X =  X, both.lambda = FALSE, n_lambda = 30, pen.type = "fused")
#' 
#' @seealso \code{\link{smoothchol}}
#' 
#' 
smoothcholCV <- function(k = 5, X, both.lambda = FALSE, lambda1_seq = NULL, 
                         lambda2_seq = NULL, max_iter = 50
                         , init.x = NULL, band = NULL, n_lambda = 60, 
                         pen.type=c("HP","fused","l1trend"), 
                         ABSTOL   = 1e-3, stand = TRUE )
{
  n = dim(X)[1]
  p = dim(X)[2]
  penalty <- match.arg(pen.type)
  mat <- matGenerate(p)
  A <- mat$A
  if (isTRUE(stand))
  {
    S =  crossprod(scale(X, center = TRUE, scale = TRUE)) / (n - 1)
  }else{
    S = var(X)
  }
  if(is.null(band))
  {
    band = p - 1
  }
  
  if(is.null(init.x))
  {
    init.x <- matrix(0, p *(p+1) / 2, 1);
    init.x[A[[1]]] = 1 / sqrt(diag(S))
  }else{
    init.x <- init.x
  }
  
  if (!is.null(lambda2_seq)){
    lambda2_seq = sort(lambda2_seq[lambda2_seq >= 0], decreasing = TRUE)
    if (length(lambda2_seq) == 0)
    {
      warning("NO positive lambda were supplied. Ignoring given lambdas.")
      lambda2_seq = exp(seq(log(1), log(0.01), length = n_lambda))
      #print(lambda_seq)
    }
  }else{
    # If lambda_seq is not supplied, calculate lambda_max
    #(the minimal value of lambda that gives zero solution), 
    #and create a sequence of length n_lambda as
    lambda2_seq = exp(seq(log(1), log(0.01), length = n_lambda))
    
  }
  ##### Train first for lambda2 and lambda1 = 0 then 
  #for lambda 1 and fit cv result for lambda2
  n_lambda = length(lambda2_seq)
  if (length(lambda1_seq) == 1)
  {
    lambda1 = lambda1_seq
  }else{
    lambda1 = 0
  }
  fold_ids = sample(1 : n) %% k + 1
  cvm =  rep(0, n_lambda) # want to have CV(lambda)
  cv_tmp = matrix(NA, n_lambda, k) 
  for (fold in 1:k){
    index = which(fold_ids == fold)  ## Stores indexes for CV 
    xtrain = X[-index,]        ## train data
    meanx = colMeans(xtrain)
    n_train = dim(xtrain)[1]
    #Create testing data xtest and ytest, everything in fold
    xtest = X[index,]          ## test data
    xtest = scale(xtest, center = meanx, scale = FALSE)
    n_test = dim(xtest)[1]
    if (isTRUE(stand))
    {
      Stest =  crossprod(scale(xtest, center = FALSE, scale = TRUE)) / (n_test - 1)
    }else{
      Stest = cov(xtest)
    }
    
    sc_fit = sc_seq(X = xtrain, lambda_seq = lambda2_seq, 
                    init.x = init.x, lambda.type = c("lambda2"), 
                    max_iter=max_iter,pen.type= penalty,
                    band = band, ABSTOL = ABSTOL, 
                    stand = stand, Stest = Stest)$tmp
#   for ( i in 1 : n_lambda)
#    {
#      omega = crossprod(Lfromx(sc_fit[, i], p))
#      cv_tmp[i, fold] = likelihood(omega, Stest)
#    }
    cv_tmp[, fold] = sc_fit
  }
  cvm = rowMeans(cv_tmp)
  se_cvm = apply(cv_tmp,1 ,sd) / k
  ibest_cvm = which.min(cvm)
  ibest_1se = min(which(cvm < cvm[ibest_cvm] + se_cvm[ibest_cvm]))
  lambda1_min = lambda1
  lambda1_1se = 0
  lambda2_min = lambda2_seq[ibest_cvm] 
  lambda2_1se = lambda2_seq[ibest_1se]
  #############Running loop for lambda1
  if (isTRUE(both.lambda))
  {
    if (!is.null(lambda1_seq)){
      lambda1_seq = sort(lambda1_seq[lambda1_seq >= 0], decreasing = TRUE)
      if (length(lambda1_seq) == 0)
      {
        warning("NO positive lambda were supplied. Ignoring given lambdas.")
        lambda1_seq = exp(seq(log(1), log(0.01), length = n_lambda))
        #print(lambda_seq)
      }
    }else{
      # If lambda_seq is not supplied, calculate lambda_max 
      #(the minimal value of lambda that gives zero solution), 
      #and create a sequence of length n_lambda as
      lambda1_seq = exp(seq(log(1), log(0.01), length = n_lambda))
      
    }
    
    n_lambda1 = length(lambda1_seq)
    cvm_lambda1 = rep(0, n_lambda1)
    cv_tmp_lambda1 = matrix(NA, n_lambda1, k)  # want to have CV(lambda)
    for (fold in 1:k){
      index = which(fold_ids == fold)  ## Stores indexes for CV 
      xtrain = X[-index,]        ## train data
      meanx = colMeans(xtrain)
      n_train = dim(xtrain)[1]
      ##Create testing data xtest and ytest, everything in fold
      xtest = X[index,]          ## test data
      xtest = scale(xtest, center = meanx, scale = FALSE)
      n_test = dim(xtest)[1]
      if (isTRUE(stand))
      {
        Stest =  crossprod(scale(xtest, center = FALSE, scale = TRUE)) / (n_test - 1)
      }else{
        onevec = matrix(1, n_test , 1)
        Stest = var(X)
      }
      
      sc_fit_lambda1 = sc_seq(X = xtrain, lambda_seq = lambda1_seq, 
                              lambda2 = lambda2_min, init.x = init.x,
                              lambda.type = c("lambda1"), max_iter=max_iter,
                              pen.type= penalty, band = band, ABSTOL = ABSTOL,
                              stand = stand, Stest = Stest )$tmp
      cv_tmp_lambda1[, fold] = sc_fit_lambda1
      
      # for ( i in 1 : n_lambda1)
      # {
      #   omega = crossprod(Lfromx(sc_fit_lambda1[, i], p))
      #   cv_tmp_lambda1[i] = cv_tmp_lambda1[i] + likelihood(omega, Stest)
      # }
    }
    cvm_lambda1 = rowMeans(cv_tmp_lambda1)
    se_cvm_l1 = apply(cv_tmp_lambda1,1 ,sd) / k
    ibest_cvm_l1 = which.min(cvm_lambda1)
    ibest_1se_l1 = min(which(cvm_lambda1 < cvm_lambda1[ibest_cvm_l1] + se_cvm_l1[ibest_cvm_l1]))
    lambda1_1se = lambda1_seq[ibest_1se_l1]
    lambda1_min = lambda1_seq[ibest_cvm_l1] 
  }

  sc_cv_fit = smoothchol(S, lambda1 = lambda1_min, lambda2 = lambda2_min,
                         max_iter = max_iter, init.x = init.x, band = band, 
                         type = penalty, ABSTOL = ABSTOL)
  return(list(lambda1_min = lambda1_min, lambda2_min = lambda2_min, 
              lambda1_1se = lambda1_1se, lambda2_1se = lambda2_1se,
              L_fit = sc_cv_fit$L, history = sc_cv_fit$history, 
              cvm = cvm, lambda2_seq = lambda2_seq))
}    

########################################################
likelihood <- function (Omega, S){
  n = dim(Omega)[1]
  # Calculate the negative log-Gaussian likelihood with
  # precision matrix Omega and sample covariance S
  return(-n * determinant(Omega, logarithm = TRUE)$modulus[1] + n * sum(S*Omega))
}


standardizeX <- function(X){
  n = dim(X)[1]
  # [ToDo] Center Y
  Xmeans = colMeans(X)      ## Column means of X
  # [ToDo] Center and scale X
  Xcenter = X - matrix(Xmeans, nrow(X), ncol(X), byrow =T) ## Centering X
  normsX  = colSums(Xcenter^2)/n                            ## Find norm of scale X for scaling
  weights = sqrt(normsX)
  Xtilde =  Xcenter * matrix(1/weights, nrow(X), ncol(X), byrow =T)
  # Return:
  # Xtilde - centered and appropriately scaled X
  # Xmeans - means of columns of X (vector)
  # weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return(list(Xtilde = Xtilde, Xmeans = Xmeans, weights = weights))
}


### smoothcholCV <- function(k = 5, X, both.lambda = FALSE, lambda1_seq = NULL, lambda2_seq = NULL, max_iter = 50 , band = NULL, n_lambda = 60, pen.type=c("HP","fused","l1trend"), ABSTOL   = 1e-3, stand = FALSE )
