/* 
   This code simulates the "supermarket model" a system with N queues 

   The Power of Two Choices on Graphs: the Pair-Approximation is
   Accurate Nicolas Gast Inria ACM MAMA Workshop,

   Copyright : Nicolas Gast (2018) -- nicolas.gast@inria.fr
 */
#include<iostream>
#include<cstdlib>
#include<cmath>

enum {NO_CHOICE=0, LOCAL_CHOICE, RND_CHOICE, SQUARE_CHOICE, ERDOS_RENYI, FIXED_DEGREE, TREE};


class histogram{
  int value_max;
  long long int *my_values;
public:
  histogram() : value_max(5) {
    my_values = (long long int *) malloc(sizeof(long long int)*value_max);
    for(int j=0;j<value_max;j++) my_values[j]=0;
  }
  void add_values(int * myValues, int N){
    for(int i=0;i<N;i++){
      while( value_max <= myValues[i] ) {
	my_values = (long long int *) realloc(my_values,sizeof(long long int)*value_max*2);
	for(int j=value_max;j<2*value_max;j++) my_values[j]=0;
	value_max*=2;
      }
      my_values[myValues[i]]++;
    }
  }
  double average(){
    int t = 0;
    double total = 0;
    for(int j=0;j<value_max;j++) {total += j*my_values[j]; t += my_values[j];}
    return total/t;
  }
  void print(){
    for(int i=0;i<value_max;i++){
      std::cout << i << " " << my_values[i] << "\n";
    }
  }
};

class queue{
  int N;
  int sqrtN;
  const int strategy;
  double lambda; // = probability of an arrival at each time step
  int *stations;
  histogram h;
  int ** neighbors; 
  int * number_neighbors; 
public:
  queue(int N, double rho, int strategy, double k=1): N(N), strategy(strategy) {
    sqrtN = sqrt(N); if (strategy == SQUARE_CHOICE) N = sqrtN*sqrtN;
    stations = (int*) malloc(sizeof(int)*N);
    for (int i=0;i<N;i++) stations[i]=0;
    lambda = rho/(1+rho);
    if (k != 1 && (strategy < ERDOS_RENYI))
      std::cerr << "*** warning: k is set but strategy != RANDOM_GRAPH ***\n";
    else{
      switch(strategy){
      case ERDOS_RENYI: generate_random_graph_erdos_renyi(k/(2*N)); break;
      case FIXED_DEGREE: generate_random_graph_fixed_degree(k); break;
      case TREE: generate_k_tree(k); break;
      }
    }
  }
  void generate_random_graph_erdos_renyi(double p){
    number_neighbors = (int*) malloc(N*sizeof(int));
    neighbors = (int**) malloc(N*sizeof(int*));
    for(int i=0;i<N;i++){
      number_neighbors[i] = 0;
      neighbors[i] = (int*) malloc(N*sizeof(int));
    }
    for(int i=0;i<N;i++)
      for(int j=0;j<N;j++)
	if (rand()/(1.+(double)RAND_MAX) <= p){
	  neighbors[i][number_neighbors[i]] = j;
	  neighbors[j][number_neighbors[j]] = i;
	  number_neighbors[i]++;
	  number_neighbors[j]++;
	}
    for(int i=0;i<N;i++)
      if (number_neighbors[i] == 0){
	number_neighbors[i]=1;
	neighbors[i][0] = i;
      }
  }
  void generate_random_graph_fixed_degree(int k){
    number_neighbors = (int*) malloc(N*sizeof(int));
    neighbors = (int**) malloc(N*sizeof(int*));
    for(int i=0;i<N;i++){
      number_neighbors[i] = 0;
      neighbors[i] = (int*) malloc(k*sizeof(int));
    }
    for(int i=0;i<N;i++){
      while(number_neighbors[i]<k){
	int potential_neighbor = rand()%N;
	if (number_neighbors[potential_neighbor]<k){
	  neighbors[i][number_neighbors[i]] = potential_neighbor;
	  number_neighbors[i]++;
	  neighbors[potential_neighbor][number_neighbors[potential_neighbor]] = i;
	  number_neighbors[potential_neighbor]++;
	}
      }
    }
  }
  void generate_k_tree(int k){
    number_neighbors = (int*) malloc(N*sizeof(int));
    neighbors = (int**) malloc(N*sizeof(int*));
    for(int i=0;i<N;i++){
      number_neighbors[i] = 0;
      neighbors[i] = (int*) malloc(k*sizeof(int));
    }
    for(int i=0;i<N/k;i++)
      for(int j=k*i+1;j<=k*i+k;j++){
	if (j<N) {
	  neighbors[i][number_neighbors[i]] = j;
	  number_neighbors[i]++;
	  neighbors[j][number_neighbors[j]] = i;
	  number_neighbors[j]++;
	}
      }
  }
    
  void iterate(long long int nbSteps){
    for(long long int t=0;t<nbSteps;t++){
      if ( rand()/(1.+(double)RAND_MAX) <= lambda) { // Arrival
	int i = rand()%N, j;
	switch(strategy){
	case NO_CHOICE:     j = i;        break;
	case LOCAL_CHOICE:  j = (i+1)%N;  break;
	case RND_CHOICE:    j = rand()%N; break;
	case SQUARE_CHOICE:
	  if (rand()%2)
	    j = (i+sqrtN)%N;
	  else
	    j = (i/sqrtN)*sqrtN + ((i+1)%sqrtN);
	  break;
	case ERDOS_RENYI: j = neighbors[i][rand()%number_neighbors[i]]; break;
	case FIXED_DEGREE: j = neighbors[i][rand()%number_neighbors[i]]; break;
	case TREE: j = neighbors[i][rand()%number_neighbors[i]]; break;
	default: std::cerr << "strategy unknown\n"; exit(-1);	  
	}
	if(stations[j]<stations[i]) stations[j]++;
	else if (stations[i]<stations[j]) stations[i]++;
	else if (rand()%2) stations[i]++; else stations[j]++;
      }
      else { // Departure
	int i = rand()%N;
	if (stations[i]>0) stations[i]--;
      }
    }
  }
  void iterate_then_hist(long long int nbSteps, int freq){
    double  percent = ((double)nbSteps)/freq/100;
    long long int next_percent = 0;
    for(long long int i=0;i<nbSteps/freq; i++){
      iterate(freq);
      if (i/percent > 50) h.add_values(stations,N);
      if (i>next_percent) {
	next_percent+=percent;
	std::cerr << "\r"<<i/percent<<"\%";
      }
    }
  }
  void print_graph(){
    if (strategy < ERDOS_RENYI) return;
    std::cerr << "on imprime le graphe\n";
    for(int i=0;i<N;i++){
      std::cerr <<"\\draw";
      for(int j=0;j<number_neighbors[i];j++) 
	std::cerr << "("<<i<<") -- ("<< neighbors[i][j] << ") ";
      std::cerr <<";\n";
    }
  }
  void print_hist(){
    h.print();
  }
  double average_queue_length(){
    return h.average();
  }
};


int main(int argc, char ** argv){
  int N=10;
  long long int nbIteration = 10000;
  double rho = 0.5;
  int choice = NO_CHOICE;
  double k = 1;
  for(int i=1;i<argc;i++){
    switch (*(argv[i])){
    case 'N': N = atoi(argv[i]+1); break;
    case 'T': nbIteration = atoll(argv[i]+1); break;
    case 'r': rho = atof(argv[i]+1); break; // note: if rho = -1, we explore values of \rho
    case 'C': choice = atoi(argv[i]+1); break;
    case 'k': k = atof(argv[i]+1); break;
    case 'p': std::cerr<<"on imprime le graphe\n";
      {queue q(N,rho,choice,k); q.print_graph(); exit(-1);} break;
    default: std::cerr << "*** Warning: unknown option " << argv[i] << " ***\n"; break;
    }
  }
  std::cerr << "N= " << N << " nbIteration="<<nbIteration <<" rho="<<rho
	    <<" choice=" << choice << " k=" << k <<"\n";
  queue q(N,rho,choice,k);
  q.iterate_then_hist(nbIteration,N);
  q.print_hist();
  
}
