
# R package for numerical integration
library(pracma)

# Function to simulate from the beta-Stacy prior law or posterior law given 
# a set of censored observations.
# Inputs: 
## n = number of replicates to generate
## time = vector of observed follow-up times; 
##        if equal to NULL, the function samples from the prior law,
##        otherwise from the posterior law conditional on the observations
## event = vector of event indicators 
##        (event[i]=0 if time[i] is censored, otherwise event[i]=1);
##        if time=NULL, then this input is ignored
## c_bs = R function to compute c(x), the precision function of the 
##        beta-Stacy prior; must be vectorized
## F_bs = R function that computes F(x), the centering cumulative distribution 
##        function of the beta-Stacy prior (assumed continuous); 
##        must be vectorized
## f_bs = R function that computes F(x), the density function of F(x);
##        must be vectorized
## grid = vector of discretization grid of time-points 
##     x[1]=0 < x[2] < x[3] < ... < x[K] 
# Output
## list of two components:
## $grid = final discretization grid (includes both the points in x and
##         the distinct event times in 'time')
## $cumul_prob = matrix of dimension n x length(grid); 
##             the i-th row containes the values of the i-th replicate 
##             S_i of the survival probabilities; s[i,j] = S_i(grid[j])
betastacy_sim <- function(n, c_bs, F_bs, f_bs, grid, time=NULL, event=NULL){
  # Defines the discretization grid
  prior <- is.null(time)
  if(is.null(time)){
    time <- 0
    event <- 0
  } 
  y <- sort(unique(c(grid, time[event==1]))) # Discretization grid
  # Computes the process parameters
  if(!prior){
    c_post <- c_bs_post(y[-1], time, event, c_bs, F_bs, f_bs)
    F_f_d <- F_bs_post_disc(y, time, event, c_bs, F_bs)
    F_f_c <- F_bs_post_cont(y, time, event, c_bs, F_bs, f_bs)
  } else {
    c_post <- c_bs(y[-1])
    F_f_d <- rep(0, times=length(y))
    F_f_c <- F_bs(y)
  }
  # Parameters of the component with fixed jump locations  
  alpha_f <- c_post * diff(F_f_d) * (1-F_f_c[-1])
  beta_f <- c_post * (1-F_f_c[-1]) * (1-F_f_d[-1])
  # Parameters of the discretized component with random jump locations  
  alpha_c <- c_post * diff(F_f_c) 
  beta_c <- c_post * (1-F_f_c[-1]) 
  # Simulates the process
  s <- matrix(0, nrow = n, ncol = length(alpha_f)+1)
  for(i in 1:nrow(s)){
    u_f <- sapply(seq_along(alpha_f), function(i) rbeta(1, alpha_f[i], beta_f[i]))
    u_c <- sapply(seq_along(alpha_c), function(i) rbeta(1, alpha_c[i], beta_c[i]))
    s[i,] <- c(1, cumprod((1-u_f)*(1-u_c)))
  }    
  # Outputs the simulation results
  out <- list(grid = y,
              cumul_probs = 1-s)
  return(out)
}


# Function that computes the continuous component F_c* of the beta-Stacy posterior
# mean function F*= 1- (1-F_c*)(1-F_d*)
# Inputs: 
## x = vector of time-points where to compute F_c*
## time = vector of observed follow-up times
## event = vector of event indicators 
##        (event[i]=0 if time[i] is censored, otherwise event[i]=1)
## c_bs = R function to compute c(x), the precision function of the 
##        beta-Stacy prior; must be vectorized
## F_bs = R function that computes F(x), the centering cumulative distribution 
##        function of the beta-Stacy prior (assumed continuous); 
##        must be vectorized
## f_bs = R function that computes F(x), the density function of F(x);
##        must be vectorized
# Output
## Vector of values of F_c*, computed over the elements of x
F_bs_post_cont <- function(x, time, event, c_bs, F_bs, f_bs){
  # Defines some auxiliary functions
  H <- function(x){ c_bs(x)*f_bs(x)/(c_bs(x)*(1-F_bs(x)) + sum(time >= x)) }
  #L <- function(x){ integrate(Vectorize(H), 0, x, subdivisions = 1000)$value }
  L <- function(x){ ifelse(x>0,pracma::gauss_kronrod(Vectorize(H), 0, x)$value,0) }
  y <- sapply(x, L)
  return( 1-exp(-y) )
}

# Function that computes the discrete component F_d* of the beta-Stacy posterior
# mean function F*= 1- (1-F_c*)(1-F_d*)
# Inputs: 
## x = vector of time-points where to compute F_d*
## time = vector of observed follow-up times
## event = vector of event indicators 
##        (event[i]=0 if time[i] is censored, otherwise event[i]=1)
## c_bs = R function to compute c(x), the precision function of the 
##        beta-Stacy prior; must be vectorized
## F_bs = R function that computes F(x), the centering cumulative distribution 
##        function of the beta-Stacy prior (assumed continuous); 
##        must be vectorized
# Output
## Vector of values of F_d*, computed over the elements of x
F_bs_post_disc <- function(x, time, event, c_bs, F_bs){
  event_times <- sort(unique(time[event == 1]))
  ff <- function(x,time, event, c_bs, F_bs){ sum(time==x & event==1)/(c_bs(x)*(1-F_bs(x)) + sum(time >= x)) }
  Delta_L <- sapply(event_times, ff, time, event, c_bs, F_bs)
  F_disc <-  1-c(1, cumprod(1-Delta_L))
  f <- stepfun(c(0,event_times), c(0,F_disc))
  return( f(x) )
}

# Function that computes the beta-Stacy posterior F*= 1- (1-F_c*)(1-F_d*)
# Inputs: 
## x = vector of time-points where to compute F*
## time = vector of observed follow-up times
## event = vector of event indicators 
##        (event[i]=0 if time[i] is censored, otherwise event[i]=1)
## c_bs = R function to compute c(x), the precision function of the 
##        beta-Stacy prior; must be vectorized
## F_bs = R function that computes F(x), the centering cumulative distribution 
##        function of the beta-Stacy prior (assumed continuous); 
##        must be vectorized
## f_bs = R function that computes F(x), the density function of F(x);
##        must be vectorized
# Output
## Vector of values of F*, computed over the elements of x
F_bs_post <- function(x, time, event, c_bs, F_bs, f_bs){
  Fc <- F_bs_post_cont(x, time, event, c_bs, F_bs, f_bs)
  Fd <- F_bs_post_disc(x, time, event, c_bs, F_bs)
  return( 1-(1-Fc)*(1-Fd) ) 
}

# Function that computes the posterior concentration function c* of the 
# beta-Stcay process
# Inputs: 
## x = vector of time-points where to compute x*
## time = vector of observed follow-up times
## event = vector of event indicators 
##        (event[i]=0 if time[i] is censored, otherwise event[i]=1)
## c_bs = R function to compute c(x), the precision function of the 
##        beta-Stacy prior; must be vectorized
## F_bs = R function that computes F(x), the centering cumulative distribution 
##        function of the beta-Stacy prior (assumed continuous); 
##        must be vectorized
## f_bs = R function that computes F(x), the density function of F(x);
##        must be vectorized
# Output
## Vector of values of x*, computed over the elements of x
c_bs_post <- function(x, time, event, c_bs, F_bs, f_bs){
  num <- c_bs(x)*(1-F_bs(x)) + sapply(x, function(x) sum(time >= x))
  den <- 1-F_bs_post(x, time, event, c_bs, F_bs, f_bs)
  return( num/den )
}

# Function to simulate time-to-event observations from the expected value 
# of the beta-Stacy posterior, i.e. the posterior estimate of the distribution 
# function (F*). It is used to implement STEP 1 of the BSB algorithm.
# The function implements the inverse probability transform integral to 
# generate observations from F*. The uniroot function is used to solve 
# the equation F*(x)=u for x (u is a sample from the uniform distribution).
# Inputs: 
## n = number of replicates to generate
## time = vector of observed follow-up times
## event = vector of event indicators 
##        (event[i]=0 if time[i] is censored, otherwise event[i]=1)
## c_bs = R function to compute c(x), the precision function of the 
##        beta-Stacy prior; must be vectorized
## F_bs = R function that computes F(x), the centering cumulative distribution 
##        function of the beta-Stacy prior (assumed continuous); 
##        must be vectorized
## f_bs = R function that computes F(x), the density function of F(x);
##        must be vectorized
## upper = upper bound for the initial interval where uniroot tries to find
##        x such that F*(x)=u (uniroot is called with the option 
##        extendInt = "upX"). Default value: 50.
## tol = tolerance used by uniroot to test if a solution has been found.
##       Default value: 1e-3
# Output
## A vector with the generated values
samp_bs_fun <- function(n, time, event, c_bs, F_bs, f_bs, 
                        upper=50, tol=1e-3){
  # Samples a from the discrete component F_d (values may be = +Inf)
  event_times <- sort(unique(time[event == 1]))
  probs_cumul <- F_bs_post_disc(event_times, time, event, c_bs, F_bs)
  probs <- diff(c(0,probs_cumul,1))
  X_disc <- sample(c(event_times,Inf), size=n, replace=TRUE, prob = probs)
  # Samples a value from the continuous component F_c using the inverse 
  # probability integral transform algorithm
  X_cont <- vector(mode="numeric", length=n)
  U <- runif(n,0,1)
  I <- order(U, decreasing=TRUE)
  U <- U[I]
  for(i in 1:n){
    X_cont[i] <- uniroot(function(x){ F_bs_post_cont(x, time, event, c_bs, F_bs, f_bs)-U[i] },
                         lower = 0, 
                         upper = ifelse(i==1,upper,X_cont[i-1]), # uses previous estimates to reduce the search interval (gives a slight performance boost)
                         tol=tol,
                         extendInt = "upX")$root
  }
  # Generated sample
  return(pmin(X_disc[I], X_cont))
}

# Function to implement the BSB algorithm
# Inputs: 
## m = number of samples from the beta-Stacy posterior mean F* to use 
##     in the algorithm; higher m provide a more accurate approximation; 
##     default: m=100
## time = vector of observed follow-up times
## event = vector of event indicators 
##        (event[i]=0 if time[i] is censored, otherwise event[i]=1)
## c_bs = R function to compute c(x), the precision function of the 
##        beta-Stacy prior; must be vectorized
## F_bs = R function that computes F(x), the centering cumulative distribution 
##        function of the beta-Stacy prior (assumed continuous); 
##        must be vectorized
## f_bs = R function that computes F(x), the density function of F(x);
##        must be vectorized
## upper, tol =  additional parameters passed to the function samp_bs_fun 
# Output
## bsb_samp = list of two elements: bsb_samp$atoms, a vector 
##            with the m atoms of the BSB distribution function, and 
##            bsb_samp$weights, the vector of the corresponding m 
##            probability weights 
bsb <- function(m, time, event, c_bs, F_bs, f_bs, upper=50, tol=1e-3){
  # Step 1
  bsb_reps <- samp_bs_fun(m, time, event, c_bs, F_bs, f_bs, upper=50, tol=1e-3)
  bsb_reps_unique <- sort(unique(bsb_reps))
  F_m <- ecdf(bsb_reps)
  # Step 2
  cc <- c_bs_post(bsb_reps_unique, time, event, c_bs, F_bs, f_bs)
  alpha <- cc * diff(F_m(c(0, bsb_reps_unique)))
  beta <- cc * (1-F_m(bsb_reps_unique))
  # Step 3
  U <- sapply(seq_along(bsb_reps_unique), function(i) rbeta(1, alpha[i], beta[i]))
  U[length(U)] <- 1
  Z <- U*c(1,cumprod(1-U[-length(U)]))
  # Step 4
  bsb_samp <- list(atoms = bsb_reps_unique, weights = Z)
  return(bsb_samp)
}
