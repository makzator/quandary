#include <ctf.hpp>
#include <float.h>
#include <math.h>
using namespace CTF;

/* Function to cast a tensor into a matrix. Here we 
   assume the ordering is such that the input tensor is symmetric in 
   blocks of 2 (i.e. ordered by qubit) */
Matrix<double> TensorToMat(Tensor<double> &input, World &dw){
  int order = input.order;
  int sz = 1;
  int64_t lens[2];
  Tensor<double> inter;

  for(int i =0;i<order;i = i+2){
    sz *= input.lens[i];
  }

  inter = input;
  lens[0] = sz;
  lens[1] = sz;
  inter.reshape(2,lens);
  Matrix<double> output = Matrix<double>(inter);
  return output;
}

/* If we assume the Hamiltonian is of the form H(t) = -i(a+a)' then
   we have an explicit form for the solution which we use to compute the
   error via timestepping. */
double compute_error_real(Tensor<double> &u, Tensor<double> &v, const double t, World &dw){
  Tensor<double> re, im;
  double err = 0.0;
  int64_t * inds;
  double * vals;
  char indname[3] = "ij";

  inds = new int64_t[2];
  vals = new double[2];

  re  = Tensor<double> (u);
  im  = Tensor<double> (v);
  re = 0.0;
  im = 0.0;

  inds[0] = 0;
  inds[1] = 3;
  vals[0] = pow(cos(t),2);
  vals[1] = 1-vals[0];
  re.write(2,inds,vals);
  inds[0] = 1;
  inds[1] = 2;
  vals[0] = -0.5*sin(2*t);
  vals[1] = 0.5*sin(2*t);;
  im.write(2,inds,vals);
  re[indname] = re[indname] - u[indname];
  im[indname] = im[indname] - v[indname];
  err += re.norm2();
  err += im.norm2();

  delete[] inds;
  delete[] vals;
  return err;
}

/* If we assume the Hamiltonian is of the form H(t) = a-a' then
   we have an explicit form for the solution which we use to compute the
   error via timestepping. */
double compute_error_im(Tensor<double> &u, Tensor<double> &v, const double t, World &dw){
  Tensor<double> re, im;
  double err = 0.0;
  int64_t * inds;
  double * vals;
  char indname[3] = "ij";

  inds = new int64_t[4];
  vals = new double[4];

  re  = Tensor<double> (u);
  im  = Tensor<double> (v);
  re = 0.0;
  im = 0.0;

  inds[0] = 0;
  inds[1] = 1;
  inds[2] = 2;
  inds[3] = 3;
  vals[0] = pow(cos(t),2);
  vals[1] = -0.5*sin(2*t);
  vals[2] = -0.5*sin(2*t);
  vals[3] = pow(sin(t),2);

  re.write(4,inds,vals);
  im = 0.0;
  // inds[0] = 1;
  // inds[1] = 2;
  // vals[0] = -0.5*sin(2*t);
  // vals[1] = 0.5*sin(2*t);;
  // im.write(2,inds,vals);
  re[indname] = re[indname] - u[indname];
  im[indname] = im[indname] - v[indname];
  err += re.norm2();
  err += im.norm2();

  delete[] inds;
  delete[] vals;
  return err;
}

/* This routine initializes tensors for the state variables and 
   the Hamiltonian matrices.*/
void Init_Tensors(int nOsc,int* nLevels,Matrix<double> * &Hsym,Matrix<double> * &Hanti,
                  Tensor<double> &u,Tensor<double> &v,
                  Tensor<double> &up,Tensor<double> &vp,
                  Tensor<double> &tmp_u,Tensor<double> &tmp_v,
                  Tensor<double> &Au, Tensor<double> &Av, 
                  World & dw){
  char ind_names[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  char psi_ind[100];
  int64_t * inds;
  int * sym_arr;
  int * asym_arr;
  int * outsym_arr;
  int maxLevel = 1;
  int n;
  int np = dw.np;
  int offset;
  double s;
  double * vals;
  int* rho_Levels;

  // Determine the largest number of energy levels for a given oscillator in the system.
  for(int i=0;i<nOsc;i++){
    maxLevel = (nLevels[i] >maxLevel) ? nLevels[i] : maxLevel;
  }

  vals = new double[maxLevel-1];
  inds = new int64_t[(maxLevel-1)];
  rho_Levels = new int[2*nOsc];

  /* Allocate and initialize Hamiltonian for each oscillator*/
  Hsym  = new Matrix<>[nOsc];
  Hanti = new Matrix<>[nOsc];
  for(int i=0;i<nOsc;i++){
    n = nLevels[i];
    Hsym[i] = Matrix<>(n,n,SY,dw);
    Hanti[i] = Matrix<>(n,n,AS,dw);
    for(int j=0;j<n-1;j++){
      vals[j] = sqrt((double)(j+1))/np;
      inds[j] = (n+1)*j + n;
    }
    Hsym[i].write(n-1,inds,vals);    // a+a'
    Hsym[i].sparsify();
    Hanti[i].write(n-1,inds,vals);   // a-a'
    Hanti[i].sparsify();
    rho_Levels[2*i] = n;
    rho_Levels[2*i+1] = n;
  }
  sym_arr = new int[2*nOsc];
  asym_arr = new int[2*nOsc];
  outsym_arr = new int[2*nOsc];
  for(int i=0;i<nOsc;i++){ 
    sym_arr[2*i] = SY;
    sym_arr[2*i+1] = NS;
    asym_arr[2*i] = AS;
    asym_arr[2*i+1] = NS;
    outsym_arr[2*i] = NS;
    outsym_arr[2*i+1] = NS;

    psi_ind[2*i] = ind_names[2*i];
    psi_ind[2*i+1] = ind_names[2*i+1];
  }

  /* Initialize density tensors */
  u      = Tensor<double> (2*nOsc, rho_Levels, outsym_arr, dw);
  v      = Tensor<double> (2*nOsc, rho_Levels, outsym_arr, dw);
  tmp_u  = Tensor<double> (2*nOsc, rho_Levels, outsym_arr, dw);
  tmp_v  = Tensor<double> (2*nOsc, rho_Levels, outsym_arr, dw);
  up     = Tensor<double> (2*nOsc, rho_Levels, outsym_arr, dw);
  vp     = Tensor<double> (2*nOsc, rho_Levels, outsym_arr, dw);
  Au     = Tensor<double> (2*nOsc, rho_Levels, outsym_arr, dw);
  Av     = Tensor<double> (2*nOsc, rho_Levels, outsym_arr, dw);


  for(int i=0;i<2*nOsc;i++){ 
    sym_arr[i] = NS;
  }
  // Tensor<double> tmp  = Tensor<double> (2*nOsc, rho_Levels, sym_arr, dw);
  // /* Initial Condition, normalized to have probability 1. */
  // u.fill_random(0.0,1.0);
  // s = u.norm2();
  // s = 1.0/(s*np);
  // u[psi_ind] = u[psi_ind]*s;
  // v = 0.0;
  // v.fill_random(0.0,1.0);
  // s = v.norm2();
  // s = 0.5/(s*np);
  // v[psi_ind] = v[psi_ind]*s;

  /* Initial condition: a canonical unit vector */
  inds[0] = 0;
  inds[1] = 2;
  // inds[2] = 2;
  vals[0] = 1.0/np;
  // vals[1] = sqrt(0.5)/np;
  // vals[2] = sqrt(0.5)/np;
  u.write(1,inds,vals);
  // v.write(1,inds,vals);

  // // Bell state
  // inds[0] = 0;
  // inds[1] = 3;
  // vals[0] = 0.5/np;
  // vals[1] = 0.5/np;
  // u.write(2,inds,vals);

  // inds[0] = 1;
  // inds[1] = 2;
  // vals[0] = 0.5/np;
  // vals[1] = -0.5/np;
  // v.write(2,inds,vals);

  // tmp = u;
  // tmp.prnt();

  // Deallocate working arrays
  delete[] rho_Levels;
  delete[] vals;
  delete[] sym_arr;
  delete[] inds;
}


/* Perform a left multiplaction by the usual control Hamiltonian for the 
 real-valued formulation. Here we try to built in functions for contractions. */
void LeftMult(int nOsc, int* nLevels, Matrix<double>* &Hsym,
             Matrix<double>* &Hanti, Tensor<double> &u,
             Tensor<double> &v,Tensor<double> &Au,
             Tensor<double> &Av, World & dw){
  char ind_names[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  char psi_ind[100];
  char Apsi_ind[100];
  char Ham_ind[2];


  for(int i=0;i<2*nOsc;i++){ 
    psi_ind[i] = ind_names[i];
    Apsi_ind[i] = ind_names[i];
  }
  // For left multiplication always sum along the second dimension
  Ham_ind[1] = ind_names[2*nOsc];

  /* Compute left multiplication by [K -S;S K] */
  for(int i=0;i<nOsc;i++){
    
    // Manipulate strings to contact along proper dimensions for each operator
    Ham_ind[0] = Apsi_ind[2*i];
    psi_ind[2*i] = ind_names[2*nOsc];
    
    Au.contract(1.0, Hanti[i],Ham_ind,u,psi_ind,1.0,Apsi_ind);
    Au.contract(1.0, Hsym[i],Ham_ind,v,psi_ind,1.0,Apsi_ind);
    Av.contract(-1.0, Hsym[i],Ham_ind,u,psi_ind,1.0,Apsi_ind);
    Av.contract(1.0, Hanti[i],Ham_ind,v,psi_ind,1.0,Apsi_ind);
    // Au[Apsi_ind] = Au[Apsi_ind] + Hanti[i][Ham_ind]*u[psi_ind] + Hsym[i][Ham_ind]*v[psi_ind];
    // Av[Apsi_ind] = Av[Apsi_ind] + Hanti[i][Ham_ind]*v[psi_ind] - Hsym[i][Ham_ind]*u[psi_ind];

    // Reset string
    psi_ind[2*i] = ind_names[2*i];
  }

  // printf("Now we cast to a matrix.\n");
  // Matrix<double> mat = TensorToMat(Au,dw);
  // mat.prnt();

}



/* Perform a right multiplication by the usual control Hamiltonian for the 
 real-valued formulation.*/
void RightMul(int nOsc, int* nLevels, Matrix<double>* &Hsym,
             Matrix<double>* &Hanti, Tensor<double> &u,
             Tensor<double> &v,Tensor<double> &Au,
             Tensor<double> &Av, World & dw){
  char ind_names[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  char psi_ind[100];
  char Apsi_ind[100];
  char Ham_ind[2];

  for(int i=0;i<2*nOsc;i++){ 
    psi_ind[i] = ind_names[i];
    Apsi_ind[i] = ind_names[i];
  }
  // For right multiplication always sum along the first dimension
  Ham_ind[0] = ind_names[2*nOsc+1];

  /* Compute left multiplication by [K -S;S K] */
  for(int i=0;i<nOsc;i++){
    // Manipulate strings to contact along proper dimensions for each operator
    Ham_ind[1] = Apsi_ind[2*i+1];
    psi_ind[2*i+1] = ind_names[2*nOsc+1];

    // Order for Lindblad
    Au.contract(-1.0, Hanti[i],Ham_ind,u,psi_ind,1.0,Apsi_ind);
    Au.contract(-1.0, Hsym[i],Ham_ind,v,psi_ind,1.0,Apsi_ind);
    Av.contract(1.0, Hsym[i],Ham_ind,u,psi_ind,1.0,Apsi_ind);
    Av.contract(-1.0, Hanti[i],Ham_ind,v,psi_ind,1.0,Apsi_ind);
    // Au[Apsi_ind] = Au[Apsi_ind] - Hanti[i][Ham_ind]*u[psi_ind] - Hsym[i][Ham_ind]*v[psi_ind];
    // Av[Apsi_ind] = Av[Apsi_ind] - Hanti[i][Ham_ind]*v[psi_ind] + Hsym[i][Ham_ind]*u[psi_ind];

    // Reset string
    psi_ind[2*i+1] = Apsi_ind[2*i+1];
  }

}



/* Attempt to do a forward solve of the Libdblad equation */
double forward_euler(double dt, double t, int nOsc,int* nLevels,Matrix<double> * &Hsym,Matrix<double> * &Hanti,
                  Tensor<double> &u,Tensor<double> &v,
                  Tensor<double> &up,Tensor<double> &vp,
                  Tensor<double> &Au,Tensor<double> &Av, 
                  World & dw){
  char ind_names[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  char psi_ind[2];
  char Apsi_ind[2];

  up = u;
  vp = v;

  // First we only wish to add in the left multiplication (i.e. ignore second index)
  for(int i=0;i<nOsc;i++){
    psi_ind[2*i] = ind_names[2*i];
    psi_ind[2*i+1] = ind_names[2*nOsc+1];
    Apsi_ind[2*i] = ind_names[2*i];
    Apsi_ind[2*i+1] = ind_names[2*i+1];
  }

  // Perform left multiplication
  LeftMult(nOsc,nLevels,Hsym,Hanti,u,v,Au,Av,dw);

  for(int i=0;i<2*nOsc;i++){
    psi_ind[i] = ind_names[i];
    Apsi_ind[i] = ind_names[i];
  }
  // Now only wish to add in the right multiplication (i.e. ignore first index)
  RightMul(nOsc,nLevels,Hsym,Hanti,u,v,Au,Av,dw);

  up.sum(dt,Au,psi_ind,1.0,Apsi_ind);
  vp.sum(dt,Av,psi_ind,1.0,Apsi_ind);

  t+=dt;

  // Swap time levels
  u = up;
  v = vp;

  return t;
}

/* For implicit time stepping, use a Neumann series to 
  compute the action of the inverse */
void Neumann(double dt, double tol, int nOsc,int* nLevels, 
             Matrix<double> * &Hsym, Matrix<double> * &Hanti,
             Tensor<double> &u,Tensor<double> &v,
             Tensor<double> &tmp_u,Tensor<double> &tmp_v,
             Tensor<double> &Au,Tensor<double> &Av, 
             World & dw){

  char ind_names[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  char psi_ind[100];
  char Apsi_ind[100];
  double coeff = 0.5*dt;
  double res_r = 1.0;
  double res_i = 1.0;
  double res   = 1.0;
  int iter = 0;

  for(int i=0;i<2*nOsc;i++){
    psi_ind[i] = ind_names[i];
    Apsi_ind[i] = ind_names[i];
  }


  u = Au;
  v = Av;
  tmp_u = 0.0;
  tmp_v = 0.0;
  // while(abs(res) > tol){
  while(iter < 2){
    LeftMult(nOsc,nLevels,Hsym,Hanti,Au,Av,tmp_u,tmp_v,dw);
    RightMul(nOsc,nLevels,Hsym,Hanti,Au,Av,tmp_u,tmp_v,dw);
    u.sum(coeff,tmp_u,psi_ind,1.0,Apsi_ind);
    v.sum(coeff,tmp_v,psi_ind,1.0,Apsi_ind);
    res_r = coeff*tmp_u.norm2();
    res_i = coeff*tmp_v.norm2();
    res = std::max(res_r,res_i);
    coeff *= 0.5*dt;
    Au = tmp_u;
    Av = tmp_v;
    tmp_u = 0.0;
    tmp_v = 0.0;
    iter++;
  }


}
/* Implicit midpoint rule time stepper */
double implicit_midpoint(double dt, double t, int nOsc,int* nLevels,
                  Matrix<double> * &Hsym,Matrix<double> * &Hanti,
                  Tensor<double> &u,Tensor<double> &v,
                  Tensor<double> &tmp_u,Tensor<double> &tmp_v,
                  Tensor<double> &up,Tensor<double> &vp,
                  Tensor<double> &Au,Tensor<double> &Av, 
                  World & dw){
  char ind_names[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  char psi_ind[100];
  char Apsi_ind[100];
  double tol = min(pow(dt,2),1e-3);

  for(int i=0;i<2*nOsc;i++){
    psi_ind[i] = ind_names[i];
    Apsi_ind[i] = ind_names[i];
  }

  // Build rhs 
  LeftMult(nOsc,nLevels,Hsym,Hanti,u,v,Au,Av,dw);
  RightMul(nOsc,nLevels,Hsym,Hanti,u,v,Au,Av,dw);
  tmp_u = Au;
  tmp_v = Av;
  up = u;
  vp = v;

  // Use a Neumann series to invert system
  Neumann(dt,tol,nOsc,nLevels,Hsym,Hanti,u,v,tmp_u,tmp_v,Au,Av,dw);

  up.sum(dt,u,psi_ind,1.0,Apsi_ind);
  vp.sum(dt,v,psi_ind,1.0,Apsi_ind);

  t+=dt;

  // Swap time levels
  u = up;
  v = vp;

  return t;
}


int main(int argc, char ** argv){
  int rank, np, pass;
  int nOsc = 3;
  int maxOsc = 3;
  int * nLevels;
  int nsteps;
  double st_time, end_time;
  double timings[4];
  double dt = 1e-4;
  double t = 0.0;
  double finalTime = 1.0;
  double err;

  FILE *extFile;
  char fName[100];

  Matrix<double> * Hsym;
  Matrix<double> * Hanti;
  Tensor<double> u;
  Tensor<double> v;
  Tensor<double> tmp_u;
  Tensor<double> tmp_v;
  Tensor<double> up;
  Tensor<double> vp;
  Tensor<double> Au;
  Tensor<double> Av;

  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  // Define ctf world object
  World dw(MPI_COMM_WORLD,argc, argv);

  if(rank == 0){
    sprintf(fName, "Time_%2.2i.txt", np);
    extFile = fopen(fName, "w");
  }

  // Loop over increasing number of "oscillators"
  for(int iter=0;iter<maxOsc;iter++){
    // nOsc = 2*(iter+1);
    // nOsc = iter+1;
    // nOsc = 3;
    nLevels = new int[nOsc];
    // for(int i=0;i<nOsc;i++) nLevels[i] = 30;
    nLevels[0] = 3;
    nLevels[1] = 3;
    nLevels[2] = 30;

    Init_Tensors(nOsc,nLevels,Hsym,Hanti,u,v,up,vp,tmp_u,tmp_v,Au,Av,dw);
    nsteps = (int) ceil(finalTime/dt);
    dt     = finalTime/((double) nsteps);
    printf("The time step is %3.2e.\n",dt);

    // nsteps = 1;
    st_time = MPI_Wtime();
    for(int k = 0;k<nsteps;k++){
      // t = forward_euler(dt, t, nOsc, nLevels, Hsym, Hanti, u, v, up, vp, Au, Av, dw);
      t = implicit_midpoint(dt,t,nOsc,nLevels,Hsym,Hanti,u,v,tmp_u,tmp_v,up,vp,Au,Av,dw);
      Au = 0.0;
      Av = 0.0;
    }
    end_time = MPI_Wtime();
    printf("Total runtime : %6.2e.\n",end_time);

    // printf("We stepped forward in time.\n");
    // u.prnt();
    // v.prnt();
    // printf("Final time : %6.2e.\n",t);
    // err = compute_error_im(u,v, t, dw);
    // // err = compute_error_real(u,v, t, dw);
    // printf("The error is %16.15e.\n",err);
    // exit(1);

    // Left multiplication (fused)
    st_time = MPI_Wtime();
    // MatMul1(nOsc,nLevels, Hsym,Hanti, u,v,Au,Av,dw);
    end_time = MPI_Wtime();
    timings[0] = end_time-st_time;

    // Left multiplication (contraction)
    st_time = MPI_Wtime();
    LeftMult(nOsc,nLevels, Hsym,Hanti, u,v,Au,Av,dw);
    end_time = MPI_Wtime();
    timings[1] = end_time-st_time;

    // Right multiplication (fused)
    st_time = MPI_Wtime();
    // RMv1(nOsc,nLevels, Hsym,Hanti, u,v,Au,Av,dw);
    end_time = MPI_Wtime();
    timings[2] = end_time-st_time;

    // Right multiplication (contraction)
    st_time = MPI_Wtime();
    RightMul(nOsc,nLevels, Hsym,Hanti, u,v,Au,Av,dw);
    end_time = MPI_Wtime();
    timings[3] = end_time-st_time;

    if(rank == 0){
      fprintf(extFile, "%4.4i %10.9e %10.9e %10.9e %10.9e\n", nOsc, timings[0], timings[1], timings[2], timings[3]);
    }
    u.~Tensor();
    v.~Tensor();
    Au.~Tensor();
    Av.~Tensor();
    tmp_u.~Tensor();
    tmp_v.~Tensor();
    delete[] nLevels;
    delete[] Hanti;
    delete[] Hsym;
  }

  if(rank==0) fclose(extFile);

  MPI_Finalize();
  return 0;
}
