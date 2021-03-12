######################## Author Information ####################################
# Author: Dhairya Jain, dhairya.jain@iitb.ac.in
# March 08, 2021

######################## File Description ######################################
# Contains the Solutions to the test given at the following link: 
# drive.google.com/file/d/0Bx2D7if2YYptOW1VLXp1bTBMOExZOFhtWWJ3UGhSd0FtUlJj/view
# Exercise 1:
# Uses the functions LogNormalParams and LogNormalPDF
# Exercise 2:
# Uses the functions fit_locdisp_mlfp and CheckConvergence

######################## Function Definitions ##################################
# Answer to Exercise 1.A
LogNormalParams <- function (ex, vx, verbose=TRUE) {
  # Computes the parameters of a log-normal distribution from expectation and
  # variance
  # 
  # Args:
  #   ex: Expectation of random variable
  #   vx: Variance of random variable
  #   verbose: If TRUE, prints the parameters; if not, not. Default is TRUE.
  #
  # Returns:
  #   The parameters mu and sigma.squared of the log normal distribution
  mu            <- log(ex ^ 2 / sqrt(vx + ex ^ 2))
  sigma.squared <- log(1 + (vx / (ex ^ 2)))
  if (verbose)
    cat("Mu = ", mu, "\nSigma Squared = ", sigma.squared, sep="")
  return (list(mu=mu, sigma.squared=sigma.squared))
}

# Used in Exercise 1.B.4
LogNormalPDF <- function (x.line, x.params) {
  # Computes the value of pdf of a log-normal distribution with given parameters 
  # at given points 
  #
  # Args:
  #   x.line: Points at which pdf is calculated
  #   x.params: Parameters mu and sigma.squared of a log-normal distribution
  #
  # Returns:
  #   The values of log-normal pdf, with parameters x.params, at x.line 
  num <- - (log(x.line) - x.params$mu) ^ 2 / (2 * x.params$sigma.squared)
  den <- x.line * sqrt(x.params$sigma.squared * 2 * pi)
  return (exp(num) / den)
}

# Used in Exercise 2.A
CheckConvergence <- function(old_u, new_u, old_s_sq, new_s_sq, threshold){
  # Checks for convergence of output parameters between successful iterations
  # based on a threshold value. The parameters converge iff both relative 
  # euclidean norm of mu and relative frobenius norm of sigma.squared are less 
  # than than the threshold value
  #
  # Args:
  #   old_u: output parameter mu after previous iteration
  #   new_u: output parameter mu after current iteration
  #   old_s_sq: output parameter sigma.squared after previous iteration
  #   new_s_sq: output parameter sigma.squared after previous iteration
  #   threshold: the minimum value of relative norm required for convergence
  #   
  # Returns:
  #   TRUE if converges, FALSE if does not
  
  # Euclidean Norm
  old_u      <- matrix(old_u, nrow=1)
  new_u      <- matrix(new_u, nrow=1)
  norm_old_u <- norm(old_u, type="F")
  norm_e_rel <- norm(new_u - old_u, type="F")
  if (norm_e_rel / norm_old_u > threshold) 
    return (FALSE)
  
  # Frobenius Norm
  norm_old_s <- norm(old_s_sq, type="F")
  norm_f_rel <- norm(new_s_sq - old_s_sq, type="F")
  if (norm_f_rel / norm_old_s > threshold) 
    return (FALSE)
  return (TRUE)
}

# Exercise 2.A
fit_locdisp_mlfp <- function(E, p, dof, threshold, verbose=TRUE){
  # Computes the maximum likelihood estimates of the location and dispersion
  # parameters under the assumption of a multivariate Student t distribution 
  # with 'dof' degrees of freedom. The routine runs iterations and checks for 
  # until convergence between consecutive iterations
  #
  # Args:
  #   E: the observations
  #   p: probability vector expressing the relative weights of the observations
  #   dof: Degree of freedom of Student t distribution
  #   threshold: The convergence threshold of the routine
  #   verbose: If TRUE, prints the parameters and the iterations; if not, not. 
  #            Default is TRUE.
  #   
  # Returns:
  #   The maximum likelihood estimates of the location and dispersion parameters
  #   under the assumption of a multivariate Student t distribution with 'dof' 
  #   degrees of freedom.
  #
  # Matrix meta data 
  t.bar <- nrow(E)  # Number of observations 
  i.bar <- ncol(E)  # Number of variables
  k     <- dof + i.bar  # A constant Used later in calculation
  # Step 0: Initialization of variables
  mu              <- colSums(E * p)
  diff            <- E - matrix(rep(mu, each=t.bar), ncol=i.bar)
  covariance.HFP  <- matrix(0, nrow=i.bar, ncol=i.bar)
  for (t in 1:t.bar)
    covariance.HFP <- covariance.HFP + p[t] * diff[t, ] %*% t(diff[t, ])
    
  if (dof > 2)
    sigma.squared <- covariance.HFP * (dof - 2) / dof
  else 
    sigma.squared <- covariance.HFP
  
  convergence <- FALSE
  iter <- 0
  
  while (!convergence) {
    iter   <- iter + 1
    # Step 1: Update Weights and FP
    sum    <- 0
    w      <- c(rep(0, t.bar))
    q      <- c(rep(0, t.bar))
    sigma.squared.inv = solve(sigma.squared)
    for(s in 1:t.bar){
      w[s] <- k / (dof + t(diff[s, ]) %*% sigma.squared.inv %*% diff[s, ])
      sum  <- sum + p[s] * w[s]
      q[s] <- p[s] * w[s] / sum
    }
    
    # Step 2: Update Output
    old_mu          <- mu
    old_s_sq        <- sigma.squared
    mu              <- colSums(E * q)
    diff            <- E - matrix(rep(mu, each=t.bar), ncol=i.bar)
    sigma.squared   <- matrix(0, nrow=i.bar, ncol=i.bar)
    for (t in 1:t.bar)
      sigma.squared <- sigma.squared + q[t] * diff[t, ] %*% t(diff[t, ])

    # Step 3: Check Convergence
    convergence <- CheckConvergence(old_mu , mu, old_s_sq, sigma.squared, 
                                    threshold)    
  }
  if (verbose)
    cat("\nConverged in ", iter, " Iterations.\nMu = ", mu, 
        "\nSigma Squared = ", sigma.squared, sep="")
  return (list(mu=mu, sigma.squared=sigma.squared))
}
######################## Statements ############################################
#
#
# EXERCISE 1
#
#
# Exercise 1.B
# Part 1: Estimate parameters for E[x] = 3, Var[x] = 5 
x.params <- LogNormalParams(3, 5)

# Part 2: Generate a large sample logn.sample
set.seed(123)
kSampleSize <- 10000
logn.sample <- rlnorm(n=kSampleSize, meanlog=x.params$mu, 
                      sdlog=sqrt(x.params$sigma.squared))
mean(logn.sample)
var(logn.sample)

# Part 3: Plot the points
plot(logn.sample, main="A large sample from LogNormal Distribution", 
     ylab="Values", type='p')
legend(80, 17, legend=c("Point"), col="black", pch=1)

# Part 4: Plot normalized histogram
hist.plot <- hist(logn.sample, breaks=50, prob=TRUE, col="grey", 
                  main="Histogram of Sample vs Exact PDF", 
                  xlab="Sample Value", ylab="Density")
x.line <- seq(min(hist.plot$breaks), max(hist.plot$breaks), length=kSampleSize)
lines(x=x.line, col="blue", lwd=3, y=LogNormalPDF(x.line, x.params))
legend("top", "Exact PDF", pch=NA, lty=1, col="blue", lwd=3)

# Part 5: Plot the empirical CDF and superimpose it on the original CDF
plot(ecdf(logn.sample), main="Empirical CDF vs Exact CDF", xlab="Sample Values", 
     ylab="Probability", col="red")
x.line <- seq(min(logn.sample), max(logn.sample), length=kSampleSize)
lines(x=x.line, 
      y=plnorm(x.line, x.params$mu, sqrt(x.params$sigma.squared)), 
      lty=2, col="blue")
legend(17, 0.9, legend=c("Empirical CDF", "Exact CDF"), lty=1:2, 
       col=c("red", "blue"))

#
#
# EXERCISE 2
#
#
# Exercise 2.B: Test 'fit_locdisp_mlfp' standard bivariate normal distribution
std.normal.bivariate <- cbind(rnorm(1000), rnorm(1000))
locdisp.mlfp.estimates <- fit_locdisp_mlfp(std.normal.bivariate, 
                                           rep(1/1000, 1000), 100, 1e-9)
