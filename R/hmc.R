#' @title Hamiltonian MC sampler
#' @param U potential energy function.
#' @param dUdq derivative of potential energy function with respect to position.
#' @param M the mass matrix (i.e.,~the covariance matrix in the kinetic energy function)
#' @param eps_gen a function which generates the leapfrog step-size \eqn{\epsilon}
#' @param L number of leapfrog steps per proposal
#' @param lower vector of lower constraints on the variables
#' @param upper vector of upper constraints on the variables
#' @examples
#'
#' ### Model setup -- 2d correlated Gaussian
#' S <- matrix(c(1,-0.98,-0.98,1),2,2)
#' n <- nrow(S)
#' Q <- chol2inv(chol(S))
#' cholQ <- chol(Q)
#' U <- function(q) 0.5 * crossprod(cholQ %*% q)
#' dUdq <- function(q) Q %*% q
#' M <- diag(1,n)
#' 
#' ### Sampler parameters -- set eps and L according to eigenvalues of covariance matrix
#' E <- eigen(S)
#' eps_gen <- function() round(min(E$values),2)
#' L = as.integer(max(E$values)/eps_gen())
#' print(paste0("eps = ",eps_gen(),". L = ",L))
#' sampler <- hmc_sampler(U = U,dUdq = dUdq, M = M,eps_gen = eps_gen,L = L)
#' 
#' ### Now sample
#' N <- 1000
#' q <- matrix(0,n,N)
#' for(i in 2:N) q[,i] <- sampler(q = q[,(i-1)])
#' plot(t(q),ylim = c(-4,4),xlim=c(-4,4))
hmc_sampler <- function(U,dUdq,M,eps_gen,L,lower=NULL,upper=NULL) {
  .check_args(U=U, dUdq = dUdq, M=M, eps_gen = eps_gen, L=L,lower=lower,upper=upper)
  cholM <- chol(M)
  Minv <- chol2inv(cholM)
  function(q) {
    stopifnot(length(q) == nrow(M))
    .hmc_sample (q, U=U, dUdq = dUdq, Minv = Minv, cholM = cholM, eps_gen=eps_gen, L=L,lower=lower,upper=upper)
  }
}


#' @title Leapfrog steps
#' @param qp data frame with columns \code{q} and \code{p} denoting the position and momentum of the particles respectively
#' @param dUdq derivative of potential energy
#' @param Minv inverse of the mass matrix
#' @param eps step-size
#' @param L number of steps to simulate
#' @examples
#' nt <- 30
#' qp <- data.frame(q=0,p=1)
#' eps_gen <- function() 0.3
#' Minv = matrix(1,1,1)
#' 
#' dUdq <- function(q) q
#' 
#' for(i in 1:(nt-1)) 
#'   qp[i+1,] <- leapfrog(qp = qp[1,], dUdq = dUdq, Minv = Minv,eps = eps_gen(),L=i)
#' 
#' plot(qp$q,qp$p)
leapfrog <- function(qp,dUdq,Minv,eps,L=1L,lower=NULL,upper=NULL) {
  
  q <- qp$q
  p <- qp$p
  
  p <- p - eps/2 * dUdq(q)
  for(i in 1:L) {
    q <- q + eps*Minv %*% p
    if(!(is.null(lower))) {
        q_id <- which(q < lower)
        q[q_id] <- lower[q_id] + (lower[q_id] - q[q_id])
        p[q_id] <- -p[q_id]
      }
        
    if(!(is.null(upper))) {
      q_id <- which(q > upper)
      q[q_id] <- upper[q_id] - (q[q_id] - upper[q_id])
      p[q_id] <- -p[q_id]
    }
    
    if(!(i==L)) p <- p - eps*dUdq(q) 
  }
  p <- p - eps * dUdq(q)/2
  
  data.frame(q=q,p=p)
}


.is.diag <- function(X) {
  all(X[lower.tri(X)] == 0, X[upper.tri(X)] == 0)
}
.sample_momentum <- function(cholM) {
  t(cholM) %*% rnorm(n = nrow(cholM))
}

.propose_state <- function(qp,U,dUdq,Minv,eps,L,lower,upper) {
  
  #Run the dynamics
  qp <- leapfrog(qp = qp ,dUdq = dUdq,Minv = Minv,eps = eps,L=L,lower=lower,upper=upper)
  # Negate momentum for summetric proposal
  qp$p <- -qp$p
  qp
}

.hmc_sample <- function(q,U,dUdq,Minv,cholM,eps_gen,L,lower,upper) {
  
  current_q <- q
  current_p <- .sample_momentum(cholM = cholM)
  current_U <- U(current_q)
  current_K <- 0.5 * t(current_p) %*% Minv %*% current_p 
  qp <- data.frame(q=current_q, p = current_p)
  
  proposed_qp <- .propose_state(qp,U = U,dUdq = dUdq,Minv = Minv,eps = eps_gen(),L = L,lower=lower,upper=upper)
  proposed_U <- U(proposed_qp$q)
  proposed_K <- 0.5 * t(proposed_qp$p) %*% Minv %*% proposed_qp$p
  
  if(log(runif(1)) < current_U - proposed_U + current_K - proposed_K) {
    return(proposed_qp$q)
  } else {
    return(current_q)
  }
}

.check_args <- function(qp=data.frame(q=c(0,0),p=c(1,0)),
                        U = function(){},
                        dUdq = function(){},
                        M = matrix(c(1,0,0,1),2,2),
                        eps_gen=function() 0.1,
                        L=10L,
                        lower=NULL,
                        upper=NULL) {
   
  stopifnot(is.data.frame(qp))
  stopifnot(names(qp) == c("q","p"))
  stopifnot(is.function(U))
  stopifnot(is.function(dUdq))
  stopifnot(is.matrix(M))
  stopifnot(is.function(eps_gen))
  stopifnot(is.integer(L))
  stopifnot(L > 0)
  stopifnot(is.null(lower) | is.numeric(lower))
  stopifnot(is.null(upper) | is.numeric(upper))
  if(!(is.null(lower)) & !.is.diag(M)) stop("Constraints not implemented for non-diagonal M")
  if(!(is.null(upper)) & !.is.diag(M)) stop("Constraints not implemented for non-diagonal M")
    
}