#' Structuring neighbourhood array
#'
#' @param m Spatial dimension vector of length 1 or 2 (integer).
#' @param sp Max spatial order (integer).
#' @param type neighbourhood type (rook or queen)
#' @param torus Logical, indicator for circular space.
#' @export
create.neighbourhood.array <- function(m = c(5,5), sp = 2, type = "queen", torus = TRUE){
  if(length(m)==1) m <- c(m,1)
  if(!(type %in% c("rook","queen"))) stop("Type must be either rook or queen")
  W <- spdep::cell2nb(m[1],m[2], type = type, torus = torus)
  Warr <- array(0, dim = c(prod(m),prod(m), sp+1))
  Warr[,,1] <- diag(prod(m))
  if(sp == 0) return(Warr)
  W <- spdep::nblag(W, max(sp,2))
  for(i in 1:(sp))
    Warr[,,i+1] <- spdep::nb2mat(W[[i]])
  return(Warr)
}


#' Simulation of STARMAGARCH
#'
#'
#' @param parameters List of parameters.
#' @param n Temporal dimension (integer).
#' @param m Spatial dimension vector of length 1 or 2 (integer).
#' @param W Neighbourhood array (optional)
#' @param burnin Temporal burnin.
#' @param type neighbourhood type (rook or queen)
#' @param torus Logical, indicator for circular space.
#' @export

#'
simSTARMAGARCH <- function(parameters, n=100, m= c(5,5), W = NULL, burnin = 300,
                           type = "rook", torus = TRUE){
  # Controlling that input makes sense:
  if(length(m)==1) m <- c(m,1)
  if(!(type %in% c("rook","queen"))) stop("Type must be either rook or queen")

  # Setting up neighbourhood array:
  if(is.null(W)){
    sp <- max(sapply(c(2,3,5,6), function(i) nrow(parameters[[i]])))
    W <- create.neighbourhood.array(m = m, sp = sp, type=type, torus=torus)
  }

  x <-  y <- matrix(0, ncol = n+burnin, nrow = prod(m[1:2]))
  sig <- x + parameters$omega
  for(t in 2:ncol(y)){
    sig[,t]<-parameters$omega +
      rowSums(sapply(1:min(ncol(parameters$alpha),t-1), function(j)
        apply(W[,,1:nrow(parameters$alpha)], 3, function(w) w %*% x[,t-j]^2)%*% parameters$alpha[,j]))+
      rowSums(sapply(1:min(ncol(parameters$beta),t-1), function(j)
       apply(W[,,1:nrow(parameters$beta)], 3, function(w) w %*% sig[,t-j])%*% parameters$beta[,j]))
    x[,t]<-sqrt(sig[,t])*rnorm(prod(m))
  }
  for(t in 2:ncol(y)){
    y[,t]<-parameters$mu +
      rowSums(sapply(1:min(ncol(parameters$phi),t-1), function(j)
        apply(W[,,1:nrow(parameters$phi)], 3, function(w) w %*% (y[,t-j]-parameters$mu))%*% parameters$phi[,j]))+
      rowSums(sapply(1:min(ncol(parameters$theta),t-1), function(j)
        apply(W[,,1:nrow(parameters$theta)], 3, function(w) w %*% x[,t-j])%*% parameters$theta[,j]))+x[,t]
  }
  return(y[,-(1:burnin)])
  #f <- MakeADFun(data = list(y=y, W=Warr,init = init),
  #               parameters = parameters, DLL = "STARMAGARCH")
}
