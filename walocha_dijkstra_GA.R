devtools::install_github('IOHprofiler/IOHexperimenter@R')
library('IOHexperimenter')

mu = 5;
pc = .6;
pm = .001;

genetic_algorithm <- function(IOHproblem){

  target_hit <- function(f, IOHproblem) return(any(f>=IOHproblem$fopt));
  
  decode <- function(P) return(P); 
  
  evaluate <- function(G, IOHproblem) return(IOHproblem$obj_func(G));  
  
  select <- function(f, index){
    # What if values are negative -> push to positive!
    f = f - min(f)+1.0001;
    # What if values are super big -> log transformationS
    f = log(f);
    f = f/sum(f);
    f = cumsum(f);
    sel_res = which(f>runif(1));
    if(!is.na(index)){
      if(index == sel_res[1]){
        winner = sel_res[2]
      } else {
        winner = sel_res[1];
      }
    } else{
      winner = sel_res[1];
    }
    if(is.na(winner)){
      winner = sample(1:length(f), 1)
    }
    return(winner);
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
  

  ###################################################### Check what IOHproblem$obj_func is
  
  budget = 50000;
  n = IOHproblem$dimension;
  fopt = -Inf;
  
  IOHproblem$set_parameters(pm);
  
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
  while(evalcount<budget && !target_hit(f, IOHproblem)){
    for(i in 1:mu){
      #select the first parent based on their fitness
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

benchmark_algorithm(user_alg=genetic_algorithm,
                    instances=c(1),
                    dimensions=c(100), 
                    functions=c(1), 
                    data.dir='./data/', 
                    params.track = 'pm',
                    algorithm.name = paste('GA-(',mu,',',mu,'),',pm,',',pc, sep = ""), 
                    algorithm.info = paste('(',mu,',',mu,') genetic algorithm, pm=',pm,'pc=',pc,sep = ""))


