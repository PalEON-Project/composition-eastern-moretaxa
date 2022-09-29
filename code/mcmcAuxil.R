require(Rcpp)
require(inline)

decide <- function(logMetropolisRatio) {
  if(is.na(logMetropolisRatio)) return(FALSE)
  if(logMetropolisRatio > 0) return(TRUE)
  if(runif(1,0,1) < exp(logMetropolisRatio)) TRUE else FALSE
}

# Ben's adaptive MCMC
adaptJump <- function(n,pjump,pgoal=NULL,max.mult=5,type='simple',i=NULL,K=NULL){
   pjump[pjump==1]=0.999
    if(is.null(pgoal)){
         const=rep(log(.44),length(n))  # optimal acc rate of .44 for d=1
            const[n==2]=log(.35)
            const[n==3]=log(.32)
            const[n>3]=log(.25)} else
                 const=log(pgoal)
    if(type=='simple'){
         if(length(n)==1)
                return(min(max.mult,max(1/max.mult,const/log(pjump)))) else
              return(pmin(max.mult,pmax(1/max.mult,const/log(pjump))))} else
                 if(type=='ben'){
                      c0=10; c1=.8
                         pgoal=exp(const)
                         gamma=c0/((i/K+3)^c1)
                         return(exp(gamma*(pjump-pgoal)))}
 }


##################################################################################
# random sampling of cells to which trees belong, given set of probabilities
##################################################################################

# R version
sampleCells <- function(probs, N, maxCells){
  # given a set of probabilities, such as from township-grid overlap
  # this efficiently does a vectorized multinomial sample
  cumProbs = t(apply(probs, 2, cumsum))
  smp = runif(N)
  id = rep(1, N)
  for(cellIndex in 1:(maxCells-1)) {
    id[smp > cumProbs[ , cellIndex]] <- cellIndex + 1
  }
  # or
  # id = 1 + rowSums(smp > probs)
  id
}

# C++ version
cppFunction('
IntegerVector sampleCells_cpp(NumericMatrix probs) {
  int nTrees = probs.ncol();
  int nCells = probs.nrow();
  IntegerVector id(nTrees);
  double smp;
  int j;
  double accumProb;
  GetRNGstate();

  for(int i = 0; i < nTrees; i++) { 
    smp = ::Rf_runif(0.0,1.0);
    accumProb = probs(0, i);
    j = 0;
    while(smp > accumProb && j < nCells - 1) {
     j++;
     accumProb += probs(j, i);
    }
    id(i) = j+1;  // convert back to 1-based indexing
  }
  PutRNGstate();
  
  return id;
}
')


# drop-in replacement for R's table() for speed
cppFunction('
IntegerVector table_cpp(IntegerVector cell, int I) {
  IntegerVector n(I);
  for(int j = 0; j < cell.size(); j++) {
    n(cell[j] - 1)++;
  }
  return n;
}
')

##################################################################################
# calculate posterior probabilities of cell membership by tree
##################################################################################

cppFunction('
NumericMatrix calcProbs(NumericMatrix W, NumericMatrix alpha_next, NumericMatrix prior, IntegerMatrix possIds, IntegerVector maxCellsByTown, IntegerVector town) {
  int nTrees = W.ncol();
  int nTaxa = W.nrow();
  int i, j, cell;
  double tmp, likSum;
  NumericVector alphaVec;
  int maxCells = possIds.nrow();
  NumericMatrix probs(maxCells, nTrees);

  for(i = 0; i < nTrees; i++) {
    tmp = 0.0;
    for(j = 0; j < maxCellsByTown[town[i]-1]; j++) {
      cell = possIds(j, i);
      alphaVec = alpha_next( _ , cell-1);
      //likSum = sum(dnorm(W( _ , i) - alphaVec, true));
      likSum = sum(pow(W( _ , i) - alphaVec, 2));
      // why was I doing this before?:
      //lik(j, i) = exp(log(lik(j, i)) / nTaxa);
      probs(j, i) = prior(j, i) * exp(-0.5*likSum);
      tmp += probs(j, i);
    }
    for(j = 0; j < maxCellsByTown[town[i] - 1]; j++) {
      probs(j, i) /= tmp;
    }
  }
  return probs;
}
')


##################################################################################
# sampling from a truncated random normal
##################################################################################

# fastest R version
rtruncnorm <- function(n, mean, lower=-Inf, upper=Inf) {
  plower <- (pnorm(lower, mean, 1))
  #pupper <- (pnorm(upper, mean, sd))
#  return(qnorm(plower+(pupper-plower)*runif(n))+mean)
  return(qnorm(plower+(pnorm(upper, mean, 1) - plower)*runif(n))+mean)
}

# single core C++ version
cppFunction('
NumericVector rtruncnorm_cpp(int n, NumericVector mean, NumericVector lower, NumericVector upper) {
  NumericVector plower(n), tmp(n);
  GetRNGstate();
  // handles either lower or upper are scalars but not both
  if(lower.size() == 1) {
    double value = lower[0];
    plower = pnorm(value - mean, true, false); 
    tmp = mean + qnorm(plower + (pnorm(upper - mean, true, false) - plower)*runif(n), true, false);
  }
  if(upper.size() == 1) {
    double value = upper[0];
    plower = pnorm(lower - mean, true, false); 
    tmp = mean + qnorm(plower + (pnorm(value - mean, true, false) - plower)*runif(n), true, false);
  }
  if(upper.size() == 1 && lower.size() == 1) {
    std::cerr << "Error: rtruncnorm_cpp not set up to handle both lower and upper truncation points as scalars. Proceeding anyway." << std::endl;
  }
  PutRNGstate();
  return tmp;
}
')

# multi-core C++ version
# note that if use this, RNG is different than all of the above approaches, which give exactly the same results
cppFunction('
NumericVector rtruncnorm_cpp_mp(int n, NumericVector mean, NumericVector lower, NumericVector upper) {
  NumericVector tmp(n);
  GetRNGstate();
  double plower;
  if(upper.size() == 1) {
    double value = upper[0];
    #pragma omp parallel for private(plower)
    for(int i = 0; i < n; ++i) {
      plower = ::Rf_pnorm5(lower[i]-mean[i], 0.0, 1.0, 1, 0);
      tmp[i] = mean[i] + ::Rf_qnorm5(plower +
      (::Rf_pnorm5(value-mean[i], 0.0, 1.0, 1, 0) - plower)*
      ::Rf_runif(0.0,1.0),  0.0, 1.0, 1, 0);
    }
    PutRNGstate();
    return tmp;
  }
  if(lower.size() == 1) {
    double value = lower[0];
    #pragma omp parallel for  private(plower)
    for(int i = 0; i < n; ++i) { 
      plower = ::Rf_pnorm5(value-mean[i], 0.0, 1.0, 1, 0);
      tmp[i] = mean[i] + ::Rf_qnorm5(plower +
      (::Rf_pnorm5(upper[i]-mean[i], 0.0, 1.0, 1, 0) - plower)*
      ::Rf_runif(0.0,1.0),  0.0, 1.0, 1, 0);
    }
    PutRNGstate();
    return tmp;
  }
  #pragma omp parallel for  private(plower)
  for(int i = 0; i < n; ++i)  {
    plower = ::Rf_pnorm5(lower[i]-mean[i], 0.0, 1.0, 1, 0);
    tmp[i] = mean[i] + ::Rf_qnorm5(plower +
    (::Rf_pnorm5(upper[i]-mean[i], 0.0, 1.0, 1, 0) - plower)*
    ::Rf_runif(0.0,1.0),  0.0, 1.0, 1, 0);
  }
  PutRNGstate();
  return tmp;
}
', plugins = c("openmp"))


##################################################################################
# functions to find row maxes
##################################################################################

# formerly rowmax2
# simple C++ implementation 
cppFunction('
NumericVector rowmax_cpp(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector maxvals(nrow);
  double max=0.0;
  int i=0;
  for(int j = 0; j < nrow; j++){
    max = x(j,ncol-1);
    for(i = 0; i < ncol; ++i){
      if(x(j,i)>max){
        max = x(j,i);
      } 
    }
    maxvals[j] = max;
  }  
  return maxvals;
}
')

# formerly rowmax2_new
# allows for excluding rows and excluding a single column here in C++ so don't need to do in R
cppFunction('
 NumericVector rowmax2_cpp(NumericMatrix x, IntegerVector inds, int excVar){
    int nrow = x.nrow(), ncol = x.ncol();
    int n = inds.size();
    int j, startCol;
    excVar -= 1;
    NumericVector maxvals(n);
    startCol = 0;
    double tmp;
    if(excVar == 0) {
      startCol = 1;
    }
    for(j = 0; j < n; j++){
     maxvals[j] = x(inds[j]-1, startCol);
    }
    startCol += 1;
    for(int i = startCol; i < ncol; i++) {
       if(i != excVar) {
         for(j = 0; j < n; j++) {
         tmp = x(inds[j]-1, i);
         if(tmp > maxvals[j]) {
           maxvals[j] = tmp;
         }
       }
     }
   }
   return maxvals;
}')

# formerly rowmax2_new_mp
cppFunction('
 NumericVector rowmax2_cpp_mp(NumericMatrix x, IntegerVector inds, int excVar){
    int nrow = x.nrow(), ncol = x.ncol();
    int n = inds.size();
    int j, startCol;
    excVar -= 1;
    NumericVector maxvals(n);
    startCol = 0;
    double tmp;
    if(excVar == 0) {
      startCol = 1;
    }
    for(j = 0; j < n; j++){
     maxvals[j] = x(inds[j]-1, startCol);
    }
    startCol += 1;
    for(int i = startCol; i < ncol; i++) {
       if(i != excVar) {
        #pragma omp parallel for private(tmp) shared(i, maxvals)
        for(j = 0; j < n; j++) {
         tmp = x(inds[j]-1, i);
         if(tmp > maxvals[j]) {
           maxvals[j] = tmp;
         }
       }
     }
   }
   return maxvals;
}', plugins = c("openmp"))





##################################################################################
# compute sums over W
##################################################################################

cppFunction('
  NumericMatrix compute_cell_sums_cpp(NumericMatrix W, NumericVector cell, int I, int P){
    int N = W.nrow();    
    int ind = 0;
    NumericMatrix returnmat(I,P);
    for(int i = 0; i < N; ++i){
    	ind = cell[i];
    	for(int j = 0; j < P; ++j){
    	  returnmat(ind-1,j) = returnmat(ind-1,j) + W(i,j);     	
    	}
    }	
    return returnmat;
  }
')

##################################################################################
# estimate cell probabilities based on latent alpha values
##################################################################################

cppFunction('
  NumericMatrix compute_cell_probabilities_cpp(NumericMatrix alpha, int numrv, int I, int P){

    GetRNGstate();
    NumericVector rvVals(P);
    NumericMatrix probs(I, P);

    double max;    
    int i, k, l;
    int maxind;


    for(i = 0; i < I; ++i){
      for(l = 0; l < P; ++l){
         probs(i,l) = 0.0;
      }
      for(k = 0; k < numrv; ++k){
        for(l = 0; l < P; ++l){
          rvVals(l) = Rf_rnorm(alpha(i,l),1);
        }
        maxind = P-1;
        max = rvVals(P-1);
        for(l = 0; l < (P-1); ++l){
          if(rvVals(l) > max){
            maxind = l;
            max = rvVals(l);
          } 
        }
        probs(i,maxind) += 1.0;
      }
      for(l = 0; l < P; ++l) {
        probs(i,l) /= numrv;
      }
    }

    PutRNGstate();
    return probs;
  }
')


cppFunction('
  NumericMatrix compute_cell_probabilities_cpp_mp(NumericMatrix alpha, int numrv, int I, int P){

    GetRNGstate();
    NumericMatrix probs(I, P);
    int i;

    #pragma omp parallel for  
    for(i = 0; i < I; ++i){
      NumericVector rvVals(P);
      double max;    
      int k, l;
      int maxind;

      for(l = 0; l < P; ++l){
         probs(i,l) = 0.0;
      }
      for(k = 0; k < numrv; ++k){
        for(l = 0; l < P; ++l){
          rvVals(l) = Rf_rnorm(alpha(i,l),1);
        }
        maxind = P-1;
        max = rvVals(P-1);
        for(l = 0; l < (P-1); ++l){
          if(rvVals(l) > max){
            maxind = l;
            max = rvVals(l);
          } 
        }
        probs(i,maxind) += 1.0;
      }
      for(l = 0; l < P; ++l) {
        probs(i,l) /= numrv;
      }
    }

    PutRNGstate();
    return probs;
  }
', plugins = c("openmp"))

           
