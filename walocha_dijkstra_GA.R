function IOH_walocha_dijkstra_GA(IOHproblem) {
	walocha_dijkstra_GA(IOHproblem$dimension, IOHproblem$obj_func, 
                   target = IOHproblem$fopt, budget = budget,
                   lambda_ = lambda_, mu_ = mu_, 
                   set_parameters = IOHproblem$set_parameters)
}


function walocha_dijkstra_GA(dimension, obj_func, 
							 target = NULL, lambda_ = 100, mu_ = 15,  
							 budget = NULL, set_parameters = NULL) {

	parents = initialize()

	while(budget>0) {
		if(budget>mu_){
			offsprings <- mating_selection(parents)

		}
		budget <- budget - lambda_
		
	}



	function initialize(){
		sapply(1:dimension, function(x) sample(c(0,1),mu_, replace=TRUE));
	}

	function mutate(mutation_rate){

	}


	function crossover(crossover_rate){

	}

	function mating_selection(parents) {
		fun_eval = apply(parents,fopt)
		return(offsprings[order(fun_eval)]
	}
}

function genetic_algorithm(IOHproblem){
  mu = ...;
  pc = ...;
  pm = ...;
  budget = ...;
  n = IOHproblem$number_of_variables();
  
  P = matrix(data=NA, nrow=mu, ncol=n);
  G = matrix(data=NA, nrow=mu, ncol=n);
  f = rep.int(NA,mu);
  
  evalcount = 0;
  for(i in 1:mu){
    P[i,] = sample(0:1,n,replace=T);
    G[i,] = decode(P[i,]);
    f[i] = evaluate(G[i,])
    evalcount = evalcount + 1;
  }
  
  # Evolution loop
  
  Pnew = rep.int(NA, mu);
  while(evalcount<budget && IOHproblem$hit_optimal(problem)){
    for(i in 1:mu){
      p1 = selection(P);
      if(unif(1) < pc) {
        p2 = selection(P);
		    Pnew[i,] = crossover(p1,p2);
      } else {
        Pnew[i,] = copy();
      }
      Pnew[i,] = mutate(Pnew[i,]);
    }
    P = Pnew;
    
    # Decode and evaluate
    for(i in 1:mu){
      G[i,] = decoding(P[i,]);
      f[i] = evaluate(G[i,])
    }
  }
  
}


# TODO: what does selection look like? 
# TODO: How to decode genome?
# TODO: How to do evaluation?


