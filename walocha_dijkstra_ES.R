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
  
  target_hit <- function(f, IOHproblem) {
    return(any(f<=IOHproblem$fopt));
  }
  evaluate <- function(P, IOHproblem){
    assign("evalcount", evalcount+dim(P)[1], envir=.GlobalEnv);
    return(IOHproblem$obj_func(P));  # return the value of the objective function
  }
  
  select <- function(f, mu){ 
    return(order(f, decreasing = TRUE)[seq(mu)]); #return indices of best
  }
  
  recombine <- function(P, f, lambda, rho, discrete=FALSE){
    if(discrete){
      if(rho==dim(P)[1]){
        Psel = P;
      } else {
        rho_best = order(f, decreasing=TRUE)
        Psel = t(P[rho_best])
      }
      offspring = apply(Psel, 2, function(x) x[sample(1:rho)]);
    } else{
      if(rho == dim(P)[1]){
        Psel = P;
        
      } else {
        rho_best = order(f, decreasing=TRUE);
        Psel = P[rho_best];
      }
      offspring = colMeans(Psel);
    }
    return(matrix(rep(offspring,lambda),nrow=lambda));
  }
  
  mutate <- function(P, tau, use_corr=FALSE,adapt=FALSE){
    if(length(sigma)>1){
      rv_global = rnorm(1)
      assign("sigma", sigma * exp(tau[1]*rv_global+tau[2]*rnorm(length(sigma))), envir=.GlobalEnv);
      P = apply(P,1, function(x) x+sigma*rnorm(length(sigma)))
    } else {
      print(sigma)
      assign("sigma", sigma*exp(tau*rnorm(1)));
      print(sigma)
      P = apply(P,1, function(x) x+sigma*rnorm(length(sigma)))
    }
    P = apply(P, c(1,2), function(x) min(max(P,-5),5))
    return(P)  
  }
  
  ###################################################### Check what IOHproblem$obj_func is
  
  budget = 10000;
  n = IOHproblem$dimension; # input dimensionality of length n. 
  fopt = -Inf;
  
  
  xopt = matrix(data=NA, nrow=mu, ncol=1);
  
  evalcount <<- 0;
  
  # initialization
  P = matrix(runif(mu*n, -5,5),nrow = mu)
  f = evaluate(P,IOHproblem);
  assign("sigma", mean(abs(f-IOHproblem$fopt)/sqrt(n)), envir=.GlobalEnv);
  if(length(sigma)>1){
    tau = c(1/sqrt(2*n), 1/sqrt(2*sqrt(n)));
  } else {
    tau = 1/sqrt(n);
  }
  
  # Evolution loop
  Pnew = matrix(data=NA, nrow=mu, ncol=n); # Offspring/children population
  while(evalcount<=(budget-lambda) && !target_hit(f, IOHproblem)){
    Prec = recombine(P, f, lambda, rho);
    Pmut = mutate(Prec, tau);
    f = evaluate(Pmut,IOHproblem);
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
