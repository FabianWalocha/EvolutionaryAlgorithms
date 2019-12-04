### Evolutionary Strategy

# Set working directory
setwd("~/Documents/R Projects/EvolutionaryAlgorithms/")
rm(list = ls())

devtools::install_github('IOHprofiler/IOHexperimenter@R')
library('IOHexperimenter')

mu = 15; #number of parents
rho = 15;
lambda = 100; #size of offspring
sig <<- 0;
discrete=FALSE;
use_plus = TRUE;
n_sig = FALSE;

if(use_plus){
  alg_name = paste('ES-(',mu,'+',lambda,')',sep="")
}else{
  alg_name = paste('ES-(',mu,'-',lambda,')', sep = "")
}

print(alg_name)

evolutionary_strategy <- function(IOHproblem){
  
  ################################################
  
  # Define functions
  
  target_hit <- function(f, IOHproblem) { # is the target value found, returns bool?
    return(IOHproblem$target_hit());
  }
  
  evaluate <- function(P, IOHproblem){ # returns the value of the objective function for a population
    assign("evalcount", evalcount+dim(P)[1], envir=.GlobalEnv);
    return(IOHproblem$obj_func(P));  
  }
  
  select <- function(f, mu){ 
    return(order(f)[seq(mu)]); #return indices of best
  }
    
  recombine <- function(P, f, lambda, rho, discrete=FALSE){
    if(is.null(dim(P))){
      # If 1 parent vector has dim NULL -> transform into 1xn matrix 
      P = t(P);
    }
    mu = dim(P)[1];
    offsprings = matrix(seq(lambda*dim(P)[2]), nrow=lambda)
    if(rho==mu){
      Psel = P;
      if(discrete){ # each feature is a copy of one of the parents
        for(idx in seq(lambda)){
          offsprings[idx,] = apply(Psel, 2, function(x) x[sample(rho,1)]) #returns recombined offspring
        }
        return(offsprings)
      } else{ # take average of the parents
        offspring = colMeans(Psel); # take the average of all parents for each feature
        return(matrix(rep(offspring,lambda),nrow=lambda)); #does this make sense? returning same off spring 100 times incase of not discrete?
      }
    } else {
      rho_rand = matrix(seq(rho*lambda), nrow= lambda)
      rho_rand = t(apply(rho_rand,1,function(x) sample(mu,rho)))
      for(idx in seq(lambda)){
        Psel = P[rho_rand[idx,],]
        if(discrete){ # each feature is a copy of one of the parents
          offsprings[idx,] = apply(Psel, 2, function(x) x[sample(rho,1)])
        } else{ # take average of the parents
          offsprings[idx,] = colMeans(Psel); # take the average of all parents for each feature
        }
      }
      return(offsprings)
    }
  }
  
  mutate <- function(P, tau){
    if(length(sig)>1){ # change boolean to global and local instead?
      rv_global = rnorm(1)
      assign("sig", sig * exp(tau[1]*rv_global+tau[2]*rnorm(length(sig))), envir=.GlobalEnv); # update all sigma's
      P = apply(P,2, function(x) x+sig*rnorm(length(sig))) # do mutation with stepsize sigma # can go out if statement?
    } else {
      assign("sig", sig*exp(tau*rnorm(1)), envir=.GlobalEnv); # update only the global sigma
      P = apply(P,2, function(x) x+sig*rnorm(length(sig))) # do mutation with stepsize sigma
    }
    if(is.null(dim(P))){
      P = t(P);
    }
    P = apply(P, c(1,2), function(x) min(max(x,-5),5)) # bound the possible values between -5 and 5
    return(P)  
  }
  
  ###################################################### 
  
  # Initialization
  
  evalcount <<- 0;
  n = IOHproblem$dimension; # input dimensionality of length n. 
  budget = 10000*n;
  fopt = Inf; # optimal value of objective function
  
  xopt = matrix(data=NA, nrow=mu, ncol=1); # find optimal matrix that gives fopt
  P = matrix(runif(mu*n, -5,5),nrow = mu) # Population of size m with dimension (features) n
  Pnew = matrix(data=NA, nrow=mu, ncol=n); # Offspring/children population
  f = evaluate(P,IOHproblem); # fitness values of the population
  
  # sigma is the mutation step size coevolving with x
  if(n_sig){
    sigmaValue <- abs(f-IOHproblem$fopt)/sqrt(n);
  }else{
    sigmaValue <- mean(abs(f-IOHproblem$fopt)/sqrt(n))
  }
  assign("sig", sigmaValue, envir=.GlobalEnv); #update the sigma globally
  
  # tau is the learning rate
  if(length(sig)>1){ # if we have individual mutation step sizes 
    tau = c(1/sqrt(2*n), 1/sqrt(2*sqrt(n))); # we have a global and coordinate wise learning rate
  } else { # if we have a global mutation step size
    tau = 1/sqrt(n); 
  }
  
  ########################################################
  
  # Evolution loop
  
  while(evalcount<=(budget-lambda) && !target_hit(f, IOHproblem)){
    Prec = recombine(P, f, lambda, rho, discrete = discrete);
    Pmut = mutate(Prec, tau);
    if(use_plus){
      Pmut = rbind(Pmut, P)
    }
    f = evaluate(Pmut,IOHproblem); 
    idxs = select(f,  mu);
    P = Pmut[idxs,];
    f = f[idxs];
    if(any(f<fopt)){
      fopt = min(f);
      if(is.null(dim(P))){
        P = t(P);
      }
      xopt = P[order(f)[1],];
    }
  }
  print(paste(min(c(fopt,f))-IOHproblem$fopt, evalcount));
}


benchmark_algorithm(user_alg=evolutionary_strategy,
                    instances=seq(5),
                    dimensions=c(2,5,20), 
                    functions=seq(24),
                    repetitions = 5,
                    suite="BBOB",
                    data.dir='./data/es/', 
                    algorithm.name = alg_name,
                    algorithm.info = alg_name)
