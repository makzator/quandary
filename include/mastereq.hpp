#include "defs.hpp"
#include "oscillator.hpp"
#include "util.hpp"
#include <petscts.h>
#include <vector>
#include <assert.h>
#include <iostream> 
#include "gate.hpp"
#pragma once


/* Define a matshell context containing pointers to data needed for applying the RHS matrix to a vector */
typedef struct {
  std::vector<int> nlevels;
  IS *isu, *isv;
  Oscillator*** oscil_vec;
  std::vector<double> xi;
  std::vector<double> detuning_freq;
  std::vector<double> collapse_time;
  std::vector<double> control_Re, control_Im;
  Mat** Ac_vec;
  Mat** Bc_vec;
  Mat *Ad, *Bd;
  Vec *Acu, *Acv, *Bcu, *Bcv;
  double time;
} MatShellCtx;


/* Define the Matrix-Vector products for the RHS MatShell */
int myMatMult_matfree_2osc(Mat RHS, Vec x, Vec y);              // Matrix free solver for 2 oscillators
int myMatMultTranspose_matfree_2Osc(Mat RHS, Vec x, Vec y);
int myMatMult_matfree_4osc(Mat RHS, Vec x, Vec y);              // Matrix free solver for 4 oscillators
int myMatMultTranspose_matfree_4Osc(Mat RHS, Vec x, Vec y);
int myMatMult_sparsemat(Mat RHS, Vec x, Vec y);                 // Sparse matrix solver
int myMatMultTranspose_sparsemat(Mat RHS, Vec x, Vec y);


/*
 * Implements the Lindblad master equation
 */
class MasterEq{

  protected:
    int dim;                 // Dimension of vectorized system = N^2
    int noscillators;        // Number of oscillators
    std::vector<int> nlevels; 
    Oscillator** oscil_vec;  // Vector storing pointers to the oscillators

    Mat RHS;                // Realvalued, vectorized systemmatrix (2N^2 x 2N^2)
    MatShellCtx RHSctx;     // MatShell context that contains data needed to apply the RHS

    Mat* Ac_vec;  // Vector of constant mats for time-varying Hamiltonian (real) 
    Mat* Bc_vec;  // Vector of constant mats for time-varying Hamiltonian (imag) 
    Mat  Ad, Bd;  // Real and imaginary part of constant system matrix

    std::vector<double> xi;             // Constants for frequencies of drift Hamiltonian
    std::vector<double> detuning_freq;  // detuning frequencies of drift Hamiltonian
    std::vector<double> collapse_time;  // Time-constants for decay and dephase operators

    /* Auxiliary stuff */
    int mpirank_petsc;   // Rank of Petsc's communicator
    int nparams_max;     // Maximum number of design parameters per oscilator 
    IS isu, isv;         // Vector strides for accessing u=Re(x), v=Im(x) 

    double *dRedp;
    double *dImdp;
    Vec Acu, Acv, Bcu, Bcv;
    int* cols;           // holding columns when evaluating dRHSdp
    PetscScalar* vals;   // holding values when evaluating dRHSdp
 
  public:
    bool usematfree;  // Flag for using matrix free solver

  public:
    MasterEq();
    MasterEq(std::vector<int> nlevels, Oscillator** oscil_vec_, const std::vector<double> xi_, std::vector<double> detuning_freq_,
             LindbladType lindbladtype_, const std::vector<double> collapse_time_, bool usematfree_);
    ~MasterEq();

    /* initialize matrices needed for applying sparse-mat solver */
    void initSparseMatSolver(LindbladType lindbladtype);

    /* Return the i-th oscillator */
    Oscillator* getOscillator(const int i);

    /* Return number of oscillators */
    int getNOscillators();

    /* Return dimension of vectorized system */
    int getDim();

    /* 
     * Uses Re and Im to build the vectorized Hamiltonian operator M = vec(-i(Hq-qH)+Lindblad). 
     * This should always be called before applying the RHS matrix.
     */
    int assemble_RHS(const double t);

    /* Access the right-hand-side matrix */
    Mat getRHS();

    /* 
     * Compute gradient of RHS wrt control parameters:
     * grad += alpha * RHS(x)^T * x_bar  
     */
    void computedRHSdp(const double t,const Vec x,const Vec x_bar, const double alpha, Vec grad);

    /* Compute reduced density operator for a sub-system defined by IDs in the oscilIDs vector */
    void createReducedDensity(const Vec rho, Vec *reduced, const std::vector<int>& oscilIDs);
    /* Derivative of reduced density computation */
    void createReducedDensity_diff(Vec rhobar, const Vec reducedbar, const std::vector<int>& oscilIDs);

    /* Set the oscillators control function parameters from global design vector x */
    void setControlAmplitudes(const Vec x);

    /* Set initial conditions 
     * In:   iinit -- index in processors range [rank * ninit_local .. (rank+1) * ninit_local - 1]
     *       ninit -- number of initial conditions 
     *       initcond_type -- type of initial condition (pure, fromfile, diagona, basis)
     *       oscilIDs -- ID of oscillators defining the subsystem for the initial conditions  
     * Out: initID -- Idenifyier for this initial condition: Element number in matrix vectorization. 
     *       rho0 -- Vector for setting initial condition 
     */
    int getRhoT0(const int iinit, const int ninit, const InitialConditionType initcond_type, const std::vector<int>& oscilIDs, Vec rho0);

};


/* Tensor contraction inlines for 2 oscillators */
inline double Hd(const double xi0, const double xi01, const double xi1, const double detuning0, const double detuning1, const int a, const int b) {
  return - xi0*M_PI * a * (a-1) - xi01*M_PI*2 * a * b - xi1*M_PI * b * (b-1) + detuning0*2*M_PI*a + detuning1*2*M_PI*b; 
};

inline double L2(double dephase0, double dephase1, const int i0, const int i1, const int i0p, const int i1p){
  return dephase0 * ( i0*i0p - 1./2. * (i0*i0 + i0p*i0p) ) + dephase1 * ( i1*i1p - 1./2. * (i1*i1 + i1p*i1p) );
};



inline double L1diag(double decay0, double decay1, const int i0, const int i1, const int i0p, const int i1p){
  return - decay0 / 2.0 * ( i0 + i0p ) - decay1 / 2.0 * ( i1 + i1p );
};


inline int TensorGetIndex(const int nlevels0, const int nlevels1,const  int i0, const int i1, int i0p, const int i1p){
  return i0*nlevels1 + i1 + (nlevels0 * nlevels1) * ( i0p * nlevels1 + i1p);
};



/* Tensor contraction inlines for 4 oscillators */
inline double Hd(const double xi0, const double xi01, const double xi02,  const double xi03, const double xi1, const double xi12, const double xi13, const double xi2, const double xi23, const double xi3, const double detuning0, const double detuning1, const double detuning2, const double detuning3, const int i0, const int i1, const int i2, const int i3) {
  return - xi0*M_PI * i0 * (i0-1)
         - xi1*M_PI * i1 * (i1-1)
         - xi2*M_PI * i2 * (i2-1)
         - xi3*M_PI * i3 * (i3-1)
         - xi01*M_PI*2 * i0 * i1
         - xi02*M_PI*2 * i0 * i2
         - xi03*M_PI*2 * i0 * i3
         - xi12*M_PI*2 * i1 * i2
         - xi13*M_PI*2 * i1 * i3
         - xi23*M_PI*2 * i2 * i3
         + detuning0*2*M_PI*i0
         + detuning1*2*M_PI*i1
         + detuning2*2*M_PI*i2
         + detuning3*2*M_PI*i3;
};

inline double L1diag(const double decay0, const double decay1, const double decay2, const double decay3,
                     const int i0,  const int i1,  const int i2,  const int i3,
                     const int i0p, const int i1p, const int i2p, const int i3p){
  return - decay0 / 2.0 * ( i0 + i0p )
         - decay1 / 2.0 * ( i1 + i1p )
         - decay2 / 2.0 * ( i2 + i2p )
         - decay3 / 2.0 * ( i3 + i3p );
};


inline double L2(const double dephase0, const double dephase1, const double dephase2, const double dephase3,
                 const int i0,  const int i1,  const int i2,  const int i3,
                 const int i0p, const int i1p, const int i2p, const int i3p){
  return   dephase0 * ( i0*i0p - 1./2. * (i0*i0 + i0p*i0p) )
         + dephase1 * ( i1*i1p - 1./2. * (i1*i1 + i1p*i1p) )
         + dephase2 * ( i2*i2p - 1./2. * (i2*i2 + i2p*i2p) )
         + dephase3 * ( i3*i3p - 1./2. * (i3*i3 + i3p*i3p) );
};

inline int TensorGetIndex(const int nlevels0, const int nlevels1, const int nlevels2, const int nlevels3,
                          const int i0,  const int i1,  const int i2,  const int i3,
                          const int i0p, const int i1p, const int i2p, const int i3p){
  // return i0*nlevels1 + i1 + (nlevels0 * nlevels1) * ( i0p * nlevels1 + i1p);
  return i0*nlevels1*nlevels2*nlevels3 + i1*nlevels2*nlevels3 + i2*nlevels3 + i3 + (nlevels0*nlevels1*nlevels2*nlevels3) * ( i0p * nlevels1*nlevels2*nlevels3 + i1p*nlevels2*nlevels3 + i2p*nlevels3 + i3p);
};
