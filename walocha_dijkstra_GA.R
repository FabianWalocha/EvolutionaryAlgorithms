# devtools::install_github('IOHprofiler/IOHexperimenter@R')
library('IOHexperimenter')


# function IOH_walocha_dijkstra_GA(IOHproblem) {
# 	walocha_dijkstra_GA(IOHproblem$dimension, IOHproblem$obj_func, 
#                    target = IOHproblem$fopt, budget = budget,
#                    lambda_ = lambda_, mu_ = mu_, 
#                    set_parameters = IOHproblem$set_parameters)
# }
# 
# 
# function walocha_dijkstra_GA(dimension, obj_func, 
# 							 target = NULL, lambda_ = 100, mu_ = 15,  
# 							 budget = NULL, set_parameters = NULL) {
# 
# 
# 	parents = initialize()
# 
# 	while(budget>0) {
# 		if(budget>mu_){
# 			offsprings <- mating_selection(parents)
# 
# 		}
# 		budget <- budget - lambda_
# 		
# 	}
# 
# 
# 
# 	function initialize(){
# 		sapply(1:dimension, function(x) sample(c(0,1),mu_, replace=TRUE));
# 	}
# 
# 	function mutate(mutation_rate){
# 
# 	}
# 
# 
# 	function crossover(crossover_rate){
# 
# 	}
# 
# 	function mating_selection(parents) {
# 		fun_eval = apply(parents,fopt)
# 		return(offsprings[order(fun_eval)]
# 	}
# }

genetic_algorithm <- function(IOHproblem){

  target_hit <- function(f) {
    return(any(f>=IOHproblem$fopt))
  }
  
  select <- function(f, index){
    if(!is.na(index)){
      f = f[-index];
    }
    old_idxs = 1:length(f);
    new_idxs = order(f, decreasing=TRUE);
    old_idxs = old_idxs[new_idxs];
    f = f[new_idxs];
    f = f/sum(f);
    f = cumsum(f);
    winner = which(f>runif(1))[1];
    return(old_idxs[winner]);
  }
  
  crossover <- function(p1, p2){
    
    if(sample(0:1,1)){
      help = p1;
      p1 = p2;
      p2 = help;
    }
    
    n = length(p1)
    intersect = sample(1:(n-1),1);
    return(c(p1[1:intersect],p2[(intersect+1):n]))
  }
  
  mutate <- function(p,pm){
    return(sapply(p,function(x){
      if(runif(1)<pm){
        return(1-x);
      } else{
        return(x);
      }
    })
    )
  }
  
  decode <- function(P){
    return(P);
  }
  
  evaluate <- function(G, IOHproblem) {
    return(IOHproblem$obj_func(G));  
  }

  mu = 10;
  pc = .7;
  pm = .01;
  budget = 50000;
  n = IOHproblem$dimension;
  fopt = -Inf;
  
  xopt = matrix(data=NA, nrow=mu, ncol=1);
  
  P = matrix(data=NA, nrow=mu, ncol=n);
  G = matrix(data=NA, nrow=mu, ncol=n);
  f = rep.int(NA,mu);
  
  evalcount = 0;
  for(i in 1:mu){
    P[i,] = sample(0:1,n,replace=T);
    G[i,] = decode(P[i,]);
    f[i] = evaluate(G[i,], IOHproblem)
    evalcount = evalcount + 1;
  }
  
  # Evolution loop
  
  Pnew = matrix(data=NA, nrow=mu, ncol=n);
  while(evalcount<budget && !target_hit(f)){
    for(i in 1:mu){
      p1 = select(f, NA);
      if(runif(1) < pc) {
        p2 = select(f, p1);
		    Pnew[i,] = matrix(data=crossover(P[p1,],P[p2,]),nrow=1);
      } else {
        Pnew[i,] = P[p1, ];
      }
      Pnew[i,] = matrix(data=mutate(Pnew[i,], pm),nrow=1);
    }
    P = Pnew;
    
    # Decode and evaluate
    for(i in 1:mu){
      G[i,] = decode(P[i,]);
      f[i] = evaluate(G[i,], IOHproblem);
      evalcount = evalcount + 1;
      # Do we have an offspring which has a better fitness than any other before?
      if(any(f>fopt)){
        fopt = max(f);
        xopt = P[order(f,decreasing = TRUE)[1],]
      }
    }
    
  }
  
}

# TODO: what does selection look like? OK
# TODO: How to decode genome? OK
# TODO: How to do evaluation?


