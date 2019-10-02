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



