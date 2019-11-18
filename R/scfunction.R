# require(l1tf)
# require(flsa)
# require(Rcpp)
# sourceCpp("SC_C_integrate.cpp")
#source("CK_matrix_Generate.R")
#require(genlasso)
#require(Matrix)

sc<-function(S, lambda1 = 0, lambda2, max_iter=50, init.x=NULL, type= c("fused","lasso","l1trend","HP"), 
             band=NULL,  ABSTOL   = 1e-4 )
{
  history = c()
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
  for (iter in 1 : max_iter)
  {
    #    progress(iter)
    oldL <- x
    ### Update Subdiagonals
    for (i in 2 :(band))
    {
      li = length(A[[i]])
      sum = as.matrix(offsubsum(A, x, i, S))
      if (type == "lasso")
      {
        
        if(i == p)
        {
          x[A[[i]]] = 0
        }else{
          B_ii_inv <- 1 / (diag(S)[1 : li])
          y<- - 2 *(sum) #+generateB(i,1,S)%*%x[A[[1]]]))
          x[A[[i]]] = 1/2 * B_ii_inv * (sign(y) * pmax( abs(y) - lambda2 ,0 ))
        }
      }
      if(type == "fused")
      {
        Bii <- diag(S)[1 : li]
        sqrt_Bii <- sqrt(Bii)
        sqrt_Bii_inv <-  1 / sqrt_Bii
        
        if (i > (p-1))
        {
          x[A[[i]]] = 0 
        }else{
          y<- - sqrt_Bii_inv * sum
          coef <- flsa(y = c(y), lambda2 = lambda2)
          coef = sqrt_Bii_inv * coef
          x[A[[i]]] = coef
        }
        if(lambda1 > 0)
        {
          x[A[[i]]] = sign(x[A[[i]]]) * pmax(abs(x[A[[i]]]) - (1 / (2 * Bii)) * lambda1, 0)
        }
      }
      if(type == "l1trend")
      {
        
        if(i>=p-2)
        {
          x_i = 0
        }else{
          Bii <- diag(S)[1 : li]
          sqrt_Bii <- sqrt(Bii)
          sqrt_Bii_inv <- 1/ sqrt_Bii
          y<- - sqrt_Bii_inv * sum
          coef <- l1tf(x = y, lambda = lambda2)
          coef = sqrt_Bii_inv * coef
          x[A[[i]]] = coef
        }
        if(lambda1>0)
        {
          x[A[[i]]] = sign(x_i) * pmax(abs(x_i)-(1/(2 * Bii)) * lambda1,0)
        }
      }
      if (type == "HP")
      {
        Bii <- diag( diag(S)[1 : li])
        y <- -2 * sum#
        if (i == p)
        {
          Bii = diag(S)[1 : li]
          D <- length(A[[i]])
          z = (1/(2 * Bii + 2 * lambda2 * D^2)) * y
          x[A[[i]]] = sign(z) * pmax( abs(y)-lambda1, 0)
        }else if(i==(p-1))
        {
          D <- getD1dSparse(li)
          DtD = crossprod(D)
          z = solve(2 * Bii+ 2 * lambda2 * DtD, y)
          x[A[[i]]] = sign(z) * pmax(abs(z) - lambda1, 0)
        }else{
          D <- getDtfSparse(li, 1)
          DtD = crossprod(D)
          z = solve(2 * Bii+ 2 * lambda2 * DtD, y)
          x[A[[i]]] = sign(z) * pmax(abs(z) - lambda1, 0)
        }
      }
      
    }
    ### Update Diagonal 
    sum_1 = offsubsum(A,x,1,S)
    diag_S = diag(S)
    x[A[[1]]] = ( - (sum_1) + sqrt((sum_1) ^ 2 + 4 * diag_S )) / (2 * diag_S)
    #### Update the new vector
    vecL = x
    ### History based on infinity norm
    history[iter]  =  max(abs(vecL - oldL));
    if (history[iter] < ABSTOL)
    {
      break
    }
  }
  if(iter==max_iter)
  {
    cat("Algorithm does not converge","\n")
  }
  
  ind <- append(A[[1]] , ND)
  x <- x[ind]
  L <- Lfromx(x, p)
  return = list(x = x, history = history,L = L)
  return(return)
}

###################################################
matGenerate<-function(p)
{
  ## This funciton generates index set A and other selection matrices
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

### This function converts x to Matrix L

Lfromx<-function(x,p)
{
  L <- matrix(0, p, p)
  diag(L) = x[1 : p]
  L[p,1] = x[ length(x) ]
  k = 0
  m <- 0
  for (i in 1:(p-2))
  {
    l= length(diag(L[-( 1 : i),-(p :(p - i + 1))]))
    k = k + l
    #    print(k)
    #    cat("L length",length(diag(L[-(1:i),-(p:(p-i+1))])),"\n")
    diag(L[-( 1 : i), - (p : (p - i + 1))]) = x[(p + m + 1) : (p + k)]
    m <- k
  }
  return(L)
}



#####################################################
#####################################################
####################################
##################################
#  if (i==p)
#  {
#  x[A[[i]]] =shrinkage(-offsubsum(A,x,i,S),lambda1)/(2*generateB(i,i,S))
#  }else{
# sqrt_Bi<-diag(1/((diag(Bi))))
# y<-as.matrix(-sqrt_Bi%*%(offsubsum(A,x,i,S)+generateB(i,1,S)%*%x[A[[1]]])/2)
# if(abs(y[1]/(diag(Bi)[1])-y[2]/(diag(Bi)[2])) <= lambda2*(1/(2*(diag(Bi)[1])^2)+1/(2*(diag(Bi)[2])^2)))
# {
#  x[A[[i]]]=c(-(2*Bi[1,1]*y[1]+2*Bi[2,2]*y[2])/(2*(Bi[1,1])^2+2*(Bi[2,2])^2),-(2*Bi[1,1]*y[1]+2*Bi[2,2]*y[2])/(2*(Bi[1,1])^2+2*(Bi[2,2])^2))
# }else if (y[1]/(diag(sqrt_Bi)[1])>y[2]/(diag(sqrt_Bi)[2]) +(1/(2*(sqrt_Bi[1,1])^2)+1/(2*(sqrt_Bi[2,2])^2)) *lambda2){
# x[A[[i]]]=c((2*(diag(Bi)[1])*y[1]-lambda2)/(2*(diag(Bi)[1])^2),(2*(diag(Bi)[2])*y[2]+lambda2)/(2*(diag(Bi)[2])^2))
# }else{
#   x[A[[i]]]=c((2*(diag(Bi)[1])*y[1]+lambda2)/(2*(diag(Bi)[1])^2),(2*(diag(Bi)[2])*y[2]-lambda2)/(2*(diag(Bi)[2])^2))
#  }
#     x[A[[i]]] = sign(x[A[[i]]])*pmax(abs(x[A[[i]]])-(1/2*diag(sqrt_Bi))*lambda1,0)
# }

# if(abs(y[1]/Bi[1,1]-2*y[2]/Bi[2,2]+y[3]/Bi[3,3])<=lambda2*(1/(2*(Bi[1,1])^2)+1/((Bi[2,2])^2)+1/(2*(Bi[3,3])^2)))
# {
#   x[A[[i]]]=c((2*Bi[1,1]*y[1]+2*Bi[2,2]*y[2]+2*Bi[3,3]*y[3])/(2*(Bi[1,1])^2+(Bi[2,2])^2+2*(Bi[3,3])^2),
#           1/2*(2*Bi[1,1]*y[1]+2*Bi[2,2]*y[2]+2*Bi[3,3]*y[3])/(2*(Bi[1,1])^2+(Bi[2,2])^2+2*(Bi[3,3])^2),
#           (2*Bi[1,1]*y[1]+2*Bi[2,2]*y[2]+2*Bi[3,3]*y[3])/(2*(Bi[1,1])^2+(Bi[2,2])^2+2*(Bi[3,3])^2))
# }else if(y[1]/Bi[1,1]+y[3]/Bi[3,3]>2*y[2]/Bi[2,2]+lambda2*(1/(2*(Bi[1,1])^2)+1/((Bi[2,2])^2)+1/(2*(Bi[3,3])^2)))
# {
#   x[A[[i]]]=c(y[1]/Bi[1,1]- lambda2/(2*(Bi[1,1])^2),
#     y[2]/Bi[2,2]+lambda2/(Bi[2,2])^2,
#     y[3]/Bi[3,3]- lambda2/(2*(Bi[3,3])^2))
# }else{
#   x[A[[i]]]=c(y[1]/Bi[1,1]+lambda2/(2*(Bi[1,1])^2),
#     y[2]/Bi[2,2]-lambda2/(Bi[2,2])^2,
#     y[3]/Bi[3,3]+ lambda2/(2*(Bi[3,3])^2))
# }

