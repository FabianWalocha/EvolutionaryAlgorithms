### Evolutionary Strategy

# Set working directory
setwd("~/Documents/R Projects/EvolutionaryAlgorithms/")
rm(list = ls())

devtools::install_github('IOHprofiler/IOHexperimenter@R-restructure')
library('IOHexperimenter')

mu = 15; #number of parents
rho = 15;
lambda = 100; #size of offspring
sigma <<- 0;

evolutionary_strategy <- function(IOHproblem){
  
  ################################################
  
  # Define functions
  
  target_hit <- function(f, IOHproblem) { # is the target value found, returns bool?
    return(any(f<=IOHproblem$fopt));
  }
  
  evaluate <- function(P, IOHproblem){ # returns the value of the objective function for a population
    assign("evalcount", evalcount+dim(P)[1], envir=.GlobalEnv);
    return(IOHproblem$obj_func(P));  
  }
  
  select <- function(f, mu){ 
    return(order(f, decreasing = TRUE)[seq(mu)]); #return indices of best
  }
  
  recombine <- function(P, f, lambda, rho, discrete=FALSE){ ## add local
    if(discrete){ # each feature is a copy of one of the parents
      if(rho==dim(P)[1]){
        Psel = P;
      } else {
        rho_best = order(f, decreasing=TRUE) # I dont think it should be ordered, but random selection?
        Psel = t(P[rho_best])
      }
      offspring = apply(Psel, 2, function(x) x[sample(1:rho)]);
    } else{ # take average of the parents
      if(rho == dim(P)[1]){
        Psel = P;
        
      } else {
        rho_best = order(f, decreasing=TRUE);
        Psel = P[rho_best];
      }
      offspring = colMeans(Psel); # take the average of all parents for each feature
    }
    return(matrix(rep(offspring,lambda),nrow=lambda)); #does this make sense? returning same off spring 100 times incase of not discrete?
  }
  
  mutate <- function(P, tau, use_corr=FALSE,adapt=FALSE){
    if(length(sigma)>1){ # change boolean to global and local instead?
      rv_global = rnorm(1)
      assign("sigma", sigma * exp(tau[1]*rv_global+tau[2]*rnorm(length(sigma))), envir=.GlobalEnv); # update all sigma's
      P = apply(P,1, function(x) x+sigma*rnorm(length(sigma))) # do mutation with stepsize sigma # can go out if statement?
    } else {
      print(sigma)
      assign("sigma", sigma*exp(tau*rnorm(1))); # update only the global sigma
      print(sigma)
      P = apply(P,1, function(x) x+sigma*rnorm(length(sigma))) # do mutation with stepsize sigma
    }
    P = apply(P, c(1,2), function(x) min(max(P,-5),5)) # bound the possible values between -5 and 5
    return(P)  
  }
  
  ###################################################### 
  
  # Initialization
  
  budget = 10000;
  evalcount <<- 0;
  n = IOHproblem$dimension; # input dimensionality of length n. 
  fopt = -Inf; # optimal value of objective function
  
  xopt = matrix(data=NA, nrow=mu, ncol=1); # find optimal matrix that gives fopt
  P = matrix(runif(mu*n, -5,5),nrow = mu) # Population of size m with dimension (features) n
  Pnew = matrix(data=NA, nrow=mu, ncol=n); # Offspring/children population
  f = evaluate(P,IOHproblem); # fitness values of the population
  
  # sigma is the mutation step size coevolving with x
  sigmaValue <- mean(abs(f-IOHproblem$fopt)/sqrt(n))
  assign("sigma", sigmaValue, envir=.GlobalEnv); #update the sigma globally
  
  # tau is the learning rate
  if(length(sigma)>1){ # if we have individual mutation step sizes 
    tau = c(1/sqrt(2*n), 1/sqrt(2*sqrt(n))); # we have a global and coordinate wise learning rate
  } else { # if we have a global mutation step size
    tau = 1/sqrt(n); 
  }
  
  ########################################################
  
  # Evolution loop
  
  while(evalcount<=(budget-lambda) && !target_hit(f, IOHproblem)){
    Prec = recombine(P, f, lambda, rho);
    Pmut = mutate(Prec, tau);
    f = evaluate(Pmut,IOHproblem); # mu komma lambda selection?
    idxs = select(f,  mu);
    Pnew = Pmut[idxs];
    f = f[idxs];
    if(any(f>fopt)){
      fopt = max(f);
      xopt = P[order(f,decreasing = TRUE)[1]];
    }
  }
}


benchmark_algorithm(user_alg=evolutionary_strategy,
                    instances=seq(5),
                    dimensions=c(2,5,20), 
                    functions=seq(24),
                    repetitions = 5,
                    suite="BBOB",
                    data.dir='./data/es/', 
                    params.track = 'sigma',
                    algorithm.name = paste('ES-(',mu,',',lambda,'),', sep = ""), 
                    algorithm.info = paste('(',mu,',',lambda,') evolutionary strategy',sep = ""))
