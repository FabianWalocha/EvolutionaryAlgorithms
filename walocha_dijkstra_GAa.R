setwd("~/Documents/R Projects/EvolutionaryAlgorithms/")
rm(list = ls())

devtools::install_github('IOHprofiler/IOHexperimenter@R')
library('IOHexperimenter')
library('sigmoid')


mu = 10;
pc = .6;
pm_fac = 0;
pm = sigmoid(pm_fac/100-7);
epsilon = 1.01;



genetic_algorithm <- function(IOHproblem){
  
  target_hit <- function(f, IOHproblem) return(any(f>=IOHproblem$fopt));
  
  encode <- function(P) return(P);  # no decode in this case so just return
  
  evaluate <- function(G, IOHproblem){
    assign("evalcount", evalcount+dim(G)[1], envir=.GlobalEnv);
    return(IOHproblem$obj_func(G));  # return the value of the objective function
  }
  
  select <- function(f, index){ 
    # make fitness positive and transform so we avoid large values
    f = f - min(f)+1.0001; # shift the values to the right to make it positive. we dont want 0 or 1 so plus a constant
    f = log(f); # Downscaling to deal with scores where sum divisor bigger than machine precision
    indices = order(f);
    f = f[indices];
    f = f/sum(f); # normalize to sum to 1
    f = cumsum(f); # to make sure we always pick a parents. the intervals between the cumulative values will still be proportionate to the fitness.
    winner = NA;
    while (is.na(winner)){
      selected_parent = which(f>runif(1))[1] # choose a random parent based on chances of fitness value
      if(!is.na(index)){ # if we have a parent to compare to
        if(index != indices[selected_parent] ){ # check if the selected parent is not the same
          winner = selected_parent
        }
      }else{ # if there is not parent to compare to yet
        winner = selected_parent
      }
    }
    return(indices[winner]);
  }
  
  crossover <- function(p1, p2){
    # we want to randomize the child we 'keep' 
    if(sample(0:1,1)){
      help = p1;
      p1 = p2;
      p2 = help;
    }
    
    n = length(p1)
    crossoverpoint = sample(1:(n-1),1);
    return(c(p1[1:crossoverpoint],p2[(crossoverpoint+1):n])) # swap the tails of the parents at the crossover point
  }
  
  mutate <- function(P,pm){
    return(abs(P-(matrix(data=runif(dim(P)[1]*dim(P)[2]),nrow=dim(P)[1])<pm)))
  }
  
  ###################################################### Check what IOHproblem$obj_func is
  
  budget = 50000;
  n = IOHproblem$dimension; #the binary vector for the genotype has length n. 
  fopt = -Inf;
  
  IOHproblem$set_parameters(pm);
  
  xopt = matrix(data=NA, nrow=mu, ncol=1);
  
  P = matrix(data=NA, nrow=mu, ncol=n); # Genotype. A vector of binary variables
  G = matrix(data=NA, nrow=mu, ncol=n); # Phenotype. A vector of binary variables
  f = rep.int(NA,mu);
  
  evalcount <<- 0;
  
  P = matrix(runif(mu*n),nrow = mu)<0.5
  
  f = evaluate(P,IOHproblem);
  
  # Evolution loop
  Pnew = matrix(data=NA, nrow=mu, ncol=n); # Offspring/children population
  while(evalcount<=(budget-mu) && !target_hit(f, IOHproblem)){
    pm = sigmoid(pm_fac/100-7)
    for(i in 1:mu){
      #select the first parent based on their fitness
      p1 = select(f, NA);
      if(runif(1) < pc) { # with chance of pc, do cross over
        p2 = select(f, p1); # select the other parent to perform cross over with
        Pnew[i,] = matrix(data=crossover(P[p1,],P[p2,]),nrow=1); # update our children with the new crossover result 
      } else {
        Pnew[i,] = P[p1, ]; # if no cross over, the child is the same as the parent
      }
    }
    P = mutate(Pnew,pm); #set our offspring population to be our new parent population
    
    f = evaluate(P,IOHproblem)
    if(any(f>fopt)){
      foptold = fopt
      fopt = max(f);
      xopt = P[order(f,decreasing = TRUE)[1],]
    }
    if(fopt/foptold<epsilon) {
      pm_fac <- min(pm_fac+1,700);
    } else {
      pm_fac <- 0;
    }
  }
  print(c('evaluations:',evalcount,'Optimality:',(IOHproblem$fopt-fopt)))
}

benchmark_algorithm(user_alg=genetic_algorithm,
                    instances=c(1),
                    dimensions=c(100), 
                    functions=seq(23), 
                    data.dir='./data/', 
                    params.track = 'pm',
                    algorithm.name = paste('GA-(',mu,',',mu,'),adaptive,',pc, sep = ""), 
                    algorithm.info = paste('(',mu,',',mu,') genetic algorithm, pm=adaptive, pc=',pc,sep = ""))