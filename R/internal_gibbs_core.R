#' Generate a cubic B-spline density basis
#' @description This function generates a cubic B-spline density basis.
#' @details \link{splineDesign} is used to generate a cubic B-spline basis.  Each B-spline is then normalised to become a B-spline density using analytical integration.  Note that the two end knots are each coincident four times.
#' @importFrom splines splineDesign
#' @importFrom Rcpp evalCpp
#' @useDynLib bsplinePsd, .registration = TRUE
#' @export
#' @param x numeric vector for which the B-spline densities are to be generated
#' @param knots knots used to generate the cubic B-spline densities
#' @return matrix of the cubic B-spline density basis
#' @seealso \link{splineDesign}
#' @examples 
#' \dontrun{
#' 
#' # Generate basis functions
#' x = seq(0, 1, length = 256)
#' knots = sort(c(0, runif(10), 1))
#' basis = dbspline(x, knots)
#' 
#' # Plot basis functions
#' plot(x, basis[1, ], type = "l", ylim = c(min(basis), max(basis)))
#' for (i in 2:nrow(basis)) lines(x, basis[i, ], col = i)
#' }
dbspline <- function(x, knots) {
  
  degree = 3  # Hard coded degree for cubic B-splines
  
  # Create coincident knots at end points
  knots.mult <- c(rep(knots[1], degree), knots, rep(knots[length(knots)], degree))
  
  # Find B-spline basis functions
  B <- splines::splineDesign(knots.mult, x, ord = degree + 1, outer.ok = TRUE)
  
  # Find normalising constant
  bs_int <- rep(NA, length = length(knots.mult) - 4)
  for (ii in 1:(length(knots.mult) - 4)) {
    part <- knots.mult[ii:(ii + 4)]
    bs_int[ii] <- AnIn1(part) + AnIn2(part) + AnIn3(part) + AnIn4(part) +
      AnIn5(part) + AnIn6(part) + AnIn7(part) + AnIn8(part)  # 8 components to analytical integral
    if (bs_int[ii] == 0) {  
      bs_int[ii] <- Inf  # Makes B.norm = 0 rather than NaN (usually caused by 0/0)
    }
  }
  
  # Convert to B-spline densities
  B.norm <- t(B) / bs_int  
  
  return(B.norm)
  
}

#' Compute unnormalised PSD using random mixture of B-splines
#' @importFrom Rcpp evalCpp
#' @useDynLib bsplinePsd, .registration = TRUE
#' @keywords internal
qpsd <- function(omega, k, v, w, u, z, recompute, db.list) {
  
  #####
  # Find weights for B-spline mixture using stick-breaking
  #####
  p <- pFromV(v)
  weight <- mixtureWeight(p, w, k)
  
  if (recompute == TRUE) {
    #####
    # Find knots for B-splines using another stick-breaking
    #####
    q <- pFromV(u)  # u here
    newk <- k - 3  # CAUTION: This is the denominator for the second DP.  Here r = 3.
    knot.diffs <- mixtureWeight(q, z, newk)
    knots <- c(0, cumsum(knot.diffs))
    
    # B-spline density matrix
    db.list <- dbspline(omega, knots)  # NOTE: *** the HARD CODED cubic spline ***
  }
  
  #####
  # B-spline mixture
  #####
  
  # Call stored matrix for k mixtures and matrix multiply with weights
  psd <- densityMixture(weight, db.list)
  epsilon <- 1e-20 
  psd <- pmax(psd, epsilon)
  psd <- psd[-c(1, length(psd))]  # COME BACK TO THIS.  Do we want to remove this?
  
  return(list(psd = psd,
              knots = knots,
              db.list = db.list))
  
}

#' Unnormalised log joint prior
#' @keywords internal
lprior <- function(k, v, w, u, z, tau, k.theta, 
                   MG, G0.alpha, G0.beta, 
                   MH, H0.alpha, H0.beta, 
                   tau.alpha, tau.beta) {
  
  logprior <- (MG - 1) * sum(log(1 - v)) +  # log prior for V's - beta(1, MG)
    sum((G0.alpha - 1) * log(w) + (G0.beta - 1) * log(1 - w)) +  # log prior for W's - beta(a, b)
    (MH - 1) * sum(log(1 - u)) +  # log prior for U's - beta(1, MH)
    sum((H0.alpha - 1) * log(z) + (H0.beta - 1) * log(1 - z)) -  # log prior for Z's - beta(a, b)
    k.theta * k ^ 2 -   # log prior for k
    (tau.alpha + 1) * log(tau) - tau.beta / tau  # log prior for tau (Inverse Gamma)
  
  return(logprior)
  
}

#' log Whittle likelihood
#' @keywords internal
llike <- function(omega, FZ, k, v, w, u, z, tau, pdgrm, recompute, db.list) {
  
  # Calculates Whittle or corrected log-likelihood (assume n even) for Gaussian errors
  
  n <- length(FZ)
  m <- n - 2  # *** NOTE: Hard-coded for cubic B-splines ***
  
  # Un-normalised PSD (defined on [0, 1])
  qq.psd <- qpsd(omega, k, v, w, u, z, recompute, db.list)
  q.psd <- qq.psd$psd
  q <- rep(NA, m) 
  q[1] <- q.psd[1]
  q[m] <- q.psd[length(q.psd)]
  q[2 * 1:(m / 2 - 1)] <- q[2 * 1:(m / 2 - 1) + 1] <- q.psd[1:(m / 2 - 1) + 1]
  
  # Normalised PSD (defined on [0, pi])
  f <- tau * q
  
  # Whittle log-likelihood
  llike <- -sum(log(f) + pdgrm[2:(n - 1)] / (f * 2 * pi)) / 2
  
  return(list(llike = llike,
              db.list = qq.psd$db.list))  
  
}

#' Unnormalised log posterior
#' @keywords internal
lpost <- function(omega, FZ, k, v, w, u, z, tau, k.theta, 
                  MG, G0.alpha, G0.beta, 
                  MH, H0.alpha, H0.beta, 
                  tau.alpha, tau.beta,
                  pdgrm, recompute, db.list) {
  
  ll <- llike(omega, FZ, k, v, w, u, z, tau, pdgrm, recompute, db.list)
  
  # Unnormalised log posterior
  lp <- ll$llike + 
    lprior(k, v, w, u, z, tau, k.theta, 
           MG, G0.alpha, G0.beta, 
           MH, H0.alpha, H0.beta, 
           tau.alpha, tau.beta)
  
  return(list(lp = lp,
              db.list = ll$db.list))
  
}
