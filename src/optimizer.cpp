#include "optimizer.hpp"

OptimProblem::OptimProblem() {
    primalbraidapp  = NULL;
    adjointbraidapp = NULL;
    objective = 0.0;
    fidelity = 0.0;
    regul = 0.0;
    mpirank_braid = 0;
    mpisize_braid = 0;
    mpirank_space = 0;
    mpisize_space = 0;
    mpirank_optim = 0;
    mpisize_optim = 0;
    mpirank_world = 0;
    mpisize_world = 0;
    printlevel = 0;
    diag_only = false;
}

OptimProblem::OptimProblem(myBraidApp* primalbraidapp_, myAdjointBraidApp* adjointbraidapp_, Gate* targate_, MPI_Comm comm_hiop_, MPI_Comm comm_init_, const std::vector<double> optim_bounds_, double optim_regul_, std::string optiminit_, std::string datadir_, int optim_printlevel_, int ilower_, int iupper_, std::string initial_cond_type_){
    primalbraidapp  = primalbraidapp_;
    adjointbraidapp = adjointbraidapp_;
    targetgate = targate_;
    comm_hiop = comm_hiop_;
    comm_init = comm_init_;
    regul = optim_regul_;
    optiminit_type = optiminit_;
    bounds = optim_bounds_;
    datadir = datadir_;
    printlevel = optim_printlevel_;
    ilower = ilower_;
    iupper = iupper_;
    initcond_type = initial_cond_type_;

    MPI_Comm_rank(primalbraidapp->comm_braid, &mpirank_braid);
    MPI_Comm_size(primalbraidapp->comm_braid, &mpisize_braid);
    MPI_Comm_rank(PETSC_COMM_WORLD, &mpirank_space);
    MPI_Comm_size(PETSC_COMM_WORLD, &mpisize_space);
    MPI_Comm_rank(comm_hiop, &mpirank_optim);
    MPI_Comm_size(comm_hiop, &mpisize_optim);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank_world);
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize_world);
    MPI_Comm_rank(comm_init, &mpirank_init);
    MPI_Comm_size(comm_init, &mpisize_init);

    if      (initcond_type.compare("all") == 0 )      diag_only = false;
    else if (initcond_type.compare("diagonal") == 0 ) diag_only = true;
    else {
      printf("Wrong initial condition type: %s \n", initcond_type.c_str());
      exit(1);
    }

    /* Open optim file */
    if (mpirank_world == 0 && printlevel > 0) {
      char filename[255];
      sprintf(filename, "%s/optim.dat", datadir.c_str());
      optimfile = fopen(filename, "w");
      fprintf(optimfile, "#iter    obj_value           fidelity              ||grad||              inf_du               ls trials \n");
    }
}

OptimProblem::~OptimProblem() {
  /* Close optim file */
  if (mpirank_world == 0 && printlevel > 0) fclose(optimfile);
}



void OptimProblem::setDesign(int n, const double* x) {

  MasterEq* mastereq = primalbraidapp->mastereq;

  /* Pass design vector x to oscillator */
  int nparam;
  double *paramRe, *paramIm;
  int j = 0;
  /* Iterate over oscillators */
  for (int ioscil = 0; ioscil < mastereq->getNOscillators(); ioscil++) {
      /* Get number of parameters of oscillator i */
      nparam = mastereq->getOscillator(ioscil)->getNParam();
      /* Get pointers to parameters of oscillator i */
      paramRe = mastereq->getOscillator(ioscil)->getParamsRe();
      paramIm = mastereq->getOscillator(ioscil)->getParamsIm();
      /* Design storage: x = (ReParams, ImParams)_iOscil */
      /* Set Re parameters */
      for (int iparam=0; iparam<nparam; iparam++) {
          paramRe[iparam] = x[j]; j++;
      }
      /* Set Im parameters */
      for (int iparam=0; iparam<nparam; iparam++) {
          paramIm[iparam] = x[j]; j++;
      }
  }
}


void OptimProblem::getDesign(int n, double* x){

  double *paramRe, *paramIm;
  int nparam;
  int j = 0;
  /* Iterate over oscillators */
  MasterEq* mastereq = primalbraidapp->mastereq;
  for (int ioscil = 0; ioscil < mastereq->getNOscillators(); ioscil++) {
      /* Get number of parameters of oscillator i */
      nparam = mastereq->getOscillator(ioscil)->getNParam();
      /* Get pointers to parameters of oscillator i */
      paramRe = mastereq->getOscillator(ioscil)->getParamsRe();
      paramIm = mastereq->getOscillator(ioscil)->getParamsIm();
      /* Design storage: x = (ReParams, ImParams)_iOscil */
      /* Set Re params */
      for (int iparam=0; iparam<nparam; iparam++) {
          x[j] = paramRe[iparam]; j++;
      }
      /* Set Im params */
      for (int iparam=0; iparam<nparam; iparam++) {
          x[j] = paramIm[iparam]; j++;
      }
  }
}

bool OptimProblem::get_prob_sizes(long long& n, long long& m) {

  // n - number of design variables 
  n = 0;
  MasterEq* mastereq = primalbraidapp->mastereq;
  for (int ioscil = 0; ioscil < mastereq->getNOscillators(); ioscil++) {
      n += 2 * mastereq->getOscillator(ioscil)->getNParam(); // Re and Im params for the i-th oscillator
  }
  
  // m - number of constraints 
  m = 0;          

  return true;
}


bool OptimProblem::get_vars_info(const long long& n, double *xlow, double* xupp, NonlinearityType* type) {

  /* Iterate over oscillators */
  int j = 0;
  MasterEq* mastereq = primalbraidapp->mastereq;
  for (int ioscil = 0; ioscil < mastereq->getNOscillators(); ioscil++) {
      /* Get number of parameters of oscillator i */
      int nparam = mastereq->getOscillator(ioscil)->getNParam();
      /* Iterate over real and imaginary part */
      for (int i = 0; i < 2 * nparam; i++) {
          xlow[j] = - bounds[ioscil];
          xupp[j] =   bounds[ioscil]; 
          j++;
      }
  }

  for (int i=0; i<n; i++) {
    type[i] =  hiopNonlinear;
  }

  return true;
}

bool OptimProblem::get_cons_info(const long long& m, double* clow, double* cupp, NonlinearityType* type){
  assert(m==0);
  return true;
}


bool OptimProblem::eval_f(const long long& n, const double* x_in, bool new_x, double& obj_value){
// bool OptimProblem::eval_f(Index n, const Number* x, bool new_x, Number& obj_value){

  if (mpirank_world == 0) printf(" EVAL F... \n");
  MasterEq* mastereq = primalbraidapp->mastereq;
  int dim = mastereq->getDim();
  double Re_local = 0.0;
  double Im_local = 0.0;
  double obj_local = 0.0;
  Vec finalstate = NULL;
  Vec initstate = NULL;
  int ninit;
  if (diag_only) ninit = sqrt(dim);
  else ninit = dim;
  /* Run simulation, only if x_in is new. Otherwise, f(x_in) has been computed already and stored in fidelity. */
  // this is fishy. check if fidelity is computed correctly in grad_f
  // if (new_x) { 

    /* Pass design vector x to oscillator */
    setDesign(n, x_in);

    /*  Iterate over initial condition */
    objective = 0.0;
    for (int iinit = ilower; iinit <= iupper; iinit++) {
      
      /* Only diagonal elements, if requested */
      if ( diag_only && iinit % (ninit+1) != 0 ) continue; 

      if (mpirank_braid == 0) printf("%d: %d FWD. ", mpirank_init, iinit);
      /* Run forward with initial condition iinit */
      initstate = primalbraidapp->PreProcess(iinit);
      if (initstate != NULL) initialCondition(iinit, initstate);
      primalbraidapp->Drive();
      finalstate = primalbraidapp->PostProcess(); // this return NULL for all but the last time processor

      /* Add to objective function */
      if (finalstate != NULL) {
        targetgate->compare(iinit, finalstate, obj_local);
        objective += obj_local;
      }
      if (mpirank_braid == 0) printf("%d: local objective: %1.14e\n", mpirank_init, obj_local);
    }
  // }

  /* Broadcast objective from last to all time processors */
  MPI_Bcast(&objective, 1, MPI_DOUBLE, mpisize_braid-1, primalbraidapp->comm_braid);

  /* Sum up objective from all initial conditions */
  double myobj = objective;
  MPI_Allreduce(&myobj, &objective, 1, MPI_DOUBLE, MPI_SUM, comm_init);

  // if (mpirank_init == 0) printf("%d: global sum objective: %1.14e\n\n", mpirank_init, objective);

  /* Compute objective 1/(2*N^2) ||W-G||_F^2 */
  objective = 1./(2. * ninit) * objective;

  /* Compute fidelity 1. - objective */
  fidelity = 1. - objective; 

  /* Add regularization objective += gamma/(2n) * ||x||^2*/
  for (int i=0; i<n; i++) {
    objective += regul / (2.0*n) * pow(x_in[i], 2.0);
  }

  /* Return objective value */
  obj_value = objective;

  return true;
}


bool OptimProblem::eval_grad_f(const long long& n, const double* x_in, bool new_x, double* gradf){
  if (mpirank_world == 0) printf(" EVAL GRAD F...\n");

  MasterEq* mastereq = primalbraidapp->mastereq;
  double obj_Re_local, obj_Im_local;
  int dim = mastereq->getDim();
  double Re_local = 0.0;
  double Im_local = 0.0;
  double obj_local = 0.0;
  Vec initstate = NULL;
  Vec finalstate = NULL;
  Vec initadjoint = NULL;
  int ninit;
  if (diag_only) ninit = sqrt(dim);
  else ninit = dim;

  /* Pass x to Oscillator */
  setDesign(n, x_in);

  /* Derivative of regularization gamma * ||x||^2 (ADD ON ONE PROC ONLY!) */
  for (int i=0; i<n; i++) {
    if (mpirank_init == 0 && mpirank_braid == 0) gradf[i] = regul / n * x_in[i];
    else gradf[i] = 0.0;
  }

  /* Derivative objective 1/(2N^2) J */
  double obj_bar = 1./(2.*ninit);

  /* Iterate over initial conditions */
  objective = 0.0;
  for (int iinit = ilower; iinit <= iupper; iinit++) {

    /* Only diagonal elements, if requested */
    if ( diag_only && iinit % (ninit+1) != 0 ) continue;

    /* --- Solve primal --- */
    if (mpirank_braid == 0) printf("%d: %d FWD. ", mpirank_init, iinit);
    initstate = primalbraidapp->PreProcess(iinit); // returns NULL if not stored on this proc
    if (initstate != NULL) initialCondition(iinit, initstate);
    primalbraidapp->Drive();
    finalstate = primalbraidapp->PostProcess(); // returns NULL if not stored on this proc

    /* Add to objective function */
    if (finalstate != NULL) {
      targetgate->compare(iinit, finalstate, obj_local);
      objective += obj_local;
    }
    if (mpirank_braid == 0) printf("%d: local objective: %1.14e\n", mpirank_init, obj_local);

    /* --- Solve adjoint --- */
    if (mpirank_braid == 0) printf("%d: %d BWD.", mpirank_init, iinit);
    initadjoint = adjointbraidapp->PreProcess(iinit); // return NULL if not stored on this proc
    if (initadjoint != NULL) 
       targetgate->compare_diff(iinit, finalstate, initadjoint, obj_bar);
    adjointbraidapp->Drive();
    adjointbraidapp->PostProcess();

    /* Add to Ipopt's gradient */
    const double* grad_ptr = adjointbraidapp->getReducedGradientPtr();
    for (int i=0; i<n; i++) {
        gradf[i] += grad_ptr[i]; 
    }
  }
  
  /* Broadcast objective from last to all processors */
  MPI_Bcast(&objective, 1, MPI_DOUBLE, mpisize_braid-1, primalbraidapp->comm_braid);

  /* Sum up objective from all initial conditions */
  double myobj = objective;
  MPI_Allreduce(&myobj, &objective, 1, MPI_DOUBLE, MPI_SUM, comm_init);
  // if (mpirank_init == 0) printf("%d: global sum objective: %1.14e\n\n", mpirank_init, objective);

  /* Compute objective 1/(2*N^2) ||W-G||_F^2 */
  objective = 1./(2.*ninit) * objective;

  /* Compute fidelity 1/N^2 |trace|^2 */
  fidelity = 1. - objective;

  /* Add regularization: objective += gamma/(2n)*||x||^2 */
  for (int i=0; i<n; i++) {
    objective += regul / (2.0*n) * pow(x_in[i], 2.0);
  }

  /* Sum up the gradient from all braid processors */
  double* mygrad = new double[n];
  for (int i=0; i<n; i++) {
    mygrad[i] = gradf[i];
  }
  MPI_Allreduce(mygrad, gradf, n, MPI_DOUBLE, MPI_SUM, primalbraidapp->comm_braid);

  /* Sum up the gradient from all initial condition processors */
  for (int i=0; i<n; i++) {
    mygrad[i] = gradf[i];
  }
  MPI_Allreduce(mygrad, gradf, n, MPI_DOUBLE, MPI_SUM, comm_init);

  /* Compute gradient norm */
  double gradnorm = 0.0;
  for (int i=0; i<n; i++) {
    gradnorm += pow(gradf[i], 2.0);
  }
  // if (mpirank_world == 0) printf("%d: ||grad|| = %1.14e\n", mpirank_init, gradnorm);

  delete [] mygrad;
    
  return true;
}


bool OptimProblem::eval_cons(const long long& n, const long long& m, const long long& num_cons, const long long* idx_cons, const double* x_in, bool new_x, double* cons) {
    assert(m==0);
    /* No constraints. Nothing to be done. */
    return true;
}


bool OptimProblem::eval_Jac_cons(const long long& n, const long long& m, const long long& num_cons, const long long* idx_cons, const double* x_in, bool new_x, double** Jac){
    assert(m==0);
    /* No constraints. Nothing to be done. */
    return true;
}

bool OptimProblem::get_starting_point(const long long &global_n, double* x0) {

  /* Set initial parameters. */
  // Do this on one processor only, then broadcast, to make sure that every processor starts with the same initial guess. 
  if (mpirank_world == 0) {
    if (optiminit_type.compare("zero") == 0)  { // init with zero
      for (int i=0; i<global_n; i++) {
        x0[i] = 0.0;
      }
    } else if ( optiminit_type.compare("random") == 0 || optiminit_type.compare("random_seed") == 0)  { // init random

      /* Set the random seed */
      if ( optiminit_type.compare("random") == 0) srand(1);  // fixed seed
      else srand(time(0)); // random seed

      /* Set to random initial guess. between [-1:1] */
      for (int i=0; i<global_n; i++) {
        x0[i] = (double) rand() / ((double)RAND_MAX);
        x0[i] = 2.*x0[i] - 1.;
      }
      /* Trimm back to the box constraints */
      int j = 0;
      MasterEq* mastereq = primalbraidapp->mastereq;
      for (int ioscil = 0; ioscil < mastereq->getNOscillators(); ioscil++) {
          int nparam = mastereq->getOscillator(ioscil)->getNParam();
          for (int i = 0; i < 2 * nparam; i++) {
              x0[j] = x0[j] * bounds[ioscil];
              j++;
          }
      }
    } else {
      /* read from file */
      read_vector(optiminit_type.c_str(), x0, global_n); 
    }
  }

  /* Broadcast the initial guess */
  MPI_Bcast(x0, global_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  /* Pass to oscillator */
  setDesign(global_n, x0);
  
  /* Flush initial control functions */
  if (mpirank_world == 0 && printlevel > 0) {
    int ntime = primalbraidapp->ntime;
    double dt = primalbraidapp->total_time / ntime;
    char filename[255];
    MasterEq* mastereq = primalbraidapp->mastereq;
    for (int ioscil = 0; ioscil < mastereq->getNOscillators(); ioscil++) {
        sprintf(filename, "%s/control_init_%02d.dat", datadir.c_str(), ioscil+1);
        mastereq->getOscillator(ioscil)->flushControl(ntime, dt, filename);
    }
  }

 
  return true;
}

/* This is called after HiOp finishes. x is LOCAL to each processor ! */
void OptimProblem::solution_callback(hiop::hiopSolveStatus status, int n, const double* x, const double* z_L, const double* z_U, int m, const double* g, const double* lambda, double obj_value) {
  
  if (mpirank_world == 0 && printlevel > 0) {
    char filename[255];
    FILE *paramfile;

    /* Print optimized parameters */
    sprintf(filename, "%s/param_optimized.dat", datadir.c_str());
    paramfile = fopen(filename, "w");
    for (int i=0; i<n; i++){
      fprintf(paramfile, "%1.14e\n", x[i]);
    }
    fclose(paramfile);

    /* Print out control functions */
    setDesign(n, x);
    int ntime = primalbraidapp->ntime;
    double dt = primalbraidapp->total_time / ntime;
    MasterEq* mastereq = primalbraidapp->mastereq;
    for (int ioscil = 0; ioscil < mastereq->getNOscillators(); ioscil++) {
        sprintf(filename, "%s/control_optimized_%02d.dat", datadir.c_str(), ioscil+1);
        mastereq->getOscillator(ioscil)->flushControl(ntime, dt, filename);
    }
  }
}


/* This is called after each iteration. x is LOCAL to each processor ! */
bool OptimProblem::iterate_callback(int iter, double obj_value, int n, const double* x, const double* z_L, const double* z_U, int m, const double* g, const double* lambda, double inf_pr, double inf_du, double mu, double alpha_du, double alpha_pr, int ls_trials) {

  /* Output */
  if (mpirank_world == 0 && printlevel > 0) {

    /* Compute current gradient norm. */
    const double* grad_ptr = adjointbraidapp->getReducedGradientPtr();
    double gnorm = 0.0;
    for (int i=0; i<n; i++) {
      gnorm += pow(grad_ptr[i], 2.0);
    }

    /* Print to optimization file */
    fprintf(optimfile, "%05d  %1.14e  %1.14e  %1.14e  %1.14e  %02d\n", iter, obj_value, fidelity, gnorm, inf_du, ls_trials);
    fflush(optimfile);

    /* Print parameters and controls to file */
    if (printlevel > 1) {
      char filename[255];

      /* Print optimized parameters */
      FILE *paramfile;
      sprintf(filename, "%s/param_iter%04d.dat", datadir.c_str(), iter);
      paramfile = fopen(filename, "w");
      for (int i=0; i<n; i++){
        fprintf(paramfile, "%1.14e\n", x[i]);
      }
      fclose(paramfile);

      /* Print control functions */
      setDesign(n, x);
      int ntime = primalbraidapp->ntime;
      double dt = primalbraidapp->total_time / ntime;
      MasterEq* mastereq = primalbraidapp->mastereq;
      for (int ioscil = 0; ioscil < mastereq->getNOscillators(); ioscil++) {
          sprintf(filename, "%s/control_iter%04d_%02d.dat", datadir.c_str(), iter, ioscil+1);
          mastereq->getOscillator(ioscil)->flushControl(ntime, dt, filename);
      }
    }
  }

  return true;
}


bool OptimProblem::get_MPI_comm(MPI_Comm& comm_out){
  comm_out = comm_hiop;
  return true;
}


int OptimProblem::initialCondition(int iinit, Vec x){

  /* Set x to i-th unit vector */
  VecZeroEntries(x); 
  VecSetValue(x, iinit, 1.0, INSERT_VALUES);
  // VecView(x, PETSC_VIEWER_STDOUT_WORLD);
  
  return 0;
}

