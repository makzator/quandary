#include "oscillator.hpp"
#include "util.hpp"
#include <petscts.h>
#pragma once

/* 
 * Abstract base class for Hamiltonian systems 
 */
class Hamiltonian{
  protected:
    int dim;                 // Dimension of vectorized system 
    int nlevels;             // Number of levels
    int noscillators;        // Number of oscillators 
    Oscillator** oscil_vec;  // Vector storing pointers to the oscillators

    Mat Re, Im;             // Real and imaginary part of Hamiltonian operator
    Mat M;                  // Realvalued, vectorized Hamiltonian operator vec(-i(Hq-qH))

  public:
    /* Default constructor sets zero */
    Hamiltonian();
    /* This constructor sets the variables and allocates Re, Im and M */
    Hamiltonian(int nlevels_, int noscillators_, Oscillator** oscil_vec_);
    ~Hamiltonian();

    /* Return dimension of vectorized system */
    int getDim();

    /* 
     * Builds up the Hamiltonian operator vectorized M = vec(-i(Hq-qH)), from Re, Im
     */
    virtual int apply(double t);

    /* Access the Hamiltonian */
    Mat getM();
};

/*
 * Hamiltonian with two oscillators 
 */
class TwoOscilHam : public Hamiltonian {

  Mat A1, A2;  // Building blocks for real part of Hamiltonian
  Mat B1, B2;  // Building blocks for imaginary part of Hamiltonian
  Mat Hd;      // Constant part of Hamiltonian matrix ("drift Hamiltonian")

  double* xi;   // xi = xi1, xi2, xi12

  public:
    TwoOscilHam(int nlevels_, double* xi, Oscillator** oscil_vec_); 
    ~TwoOscilHam(); 

    /* Helper function for constructing building blocks */
    int BuildingBlock(Mat C, int sign, int k, int m);

    /* Apply the Hamiltonian operator */
    virtual int apply(double t);
 
};


/* 
 * Ander's testcase with analytic solution 
 */
class AnalyticHam : public Hamiltonian {


};

/* 
 * Compute the analytic solution for the 2-oscillator, 2-levels test case.
 */
PetscErrorCode ExactSolution(PetscReal t,Vec s, PetscReal freq);


/*
 *  Set the initial condition at time t_0 to the analytic solution 
 *  of the 2-level, 2-oscillator case.
 */
PetscErrorCode InitialConditions(Vec x,PetscReal freq);



/* Real part for Oscillator 1 of analytic solution */
PetscScalar F1(PetscReal t,PetscReal freq);


/* Imaginary part for Oscillator 2 of analytic solution */
PetscScalar G2(PetscReal t,PetscReal freq);

