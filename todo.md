[x] Choose a time-horizon [0,T] that involve about 5-10 periods (e.g. 10us)
[ ] Find an appropriate time-step size to resolve that interval. 
     [ ] Create a reference solution using a very very small \Delta t, store the solution
     [ ] For a sequence of \Delta t from 'big' to 'small', compute the (relative) errors of the solution at final time T 
     [ ] Plot the errors over \Delta t, should be decreasing with second order 
     [ ] Choose error value you like, and fix the time-step size that corresponds to this.
[ ] Perform a systematic parameter study for XBraid
    [ ]  Choose a range of coarsening factors, testing in particular the difference between even and odd ones, since we assume there might be a significant different due to the complex-only eigenvalues. Choose a range of maximum number of braid levels. Use F-relaxation, later also test FCF-relaxation to see if it improves things
    [ ] Choose a stopping criterion for XBraid so that the absolute error is of the same order as the error of the time-stepper as above.
    [ ] Run XBraid with all possible combination, comparing the number of iterations XBraid takes. 
    [ ] Visualize those in a table, or another format of your choice, and document
    [ ] Plot XBraids convergence (residual over number of iterations) for the best one
    [ ] Run scaling tests for the best one, i.e. run this case on 2^x cores, for x=1,2,\dots, while measuring the runtime. Make sure to allocate entire nodes when submitting those jobs, so that results are not bias by other computations on that same node. 
[ ] Create plots according to Ben's 2-level convergence theory (see email Friday 7/10, 3:05pm).
    [ ] Choose maxlevels=2, choose a coarsening factor c from above, or just test e.g. c=2, c=5
    [ ] Evaluate the eigenvalues of Phi and Psi, for decreasing \Delta t (all being small enough to get down to the chosen time-stepping error though), and constant controls p,q. Here, Psi is the same as Phi but only uses a time-step size c*dt
    [ ] Evaluate Ben's formula for computing the convergence factor of two-level error propagation using the above eigenvalues, plot that over dt*(eigenvalue), for all eigenvalues


We will do the above steps for all of the three cases:
* Alice-only with na=2 levels
* Alice-Cavity with 3x20 levels, constant controls
* Alice-Cavity with 3x20 levels, using time-varying controls given by the optimized parameters provided by Stefanie's optimization
