
using Optim
	N = 4
	#N = 2
	
	Nguard = 3

	Ntot = N + Nguard
	
	Ident = Matrix{Float64}(I, Ntot, Ntot)   
	utarget = Ident[1:Ntot,1:N]

	utarget[:,3] = Ident[:,4]
	utarget[:,4] = Ident[:,3]

	cfl = 0.05

	T = 150.0

	testadjoint = 0
	maxpar =0.09
	
	params = objfunc.parameters(N,Nguard,T,testadjoint,maxpar,cfl, utarget)
	

	pcof0  = zeros(351) 
	#pcof0 = (rand(250) .- 0.5).*maxpar*0.1
	order = 2
	weight = 2
	penaltyweight = 2

    function f(pcof)
    #@show(pcof)
     f =objfunc.traceobjgrad(pcof,params,order,false,false,weight,penaltyweight )
    # @show(f)
     return f[1]
     end

    function g!(G,pcof,params,order)
    	objf, Gtemp = objfunc.traceobjgrad(pcof,params,order,false, true, weight, penaltyweight )
    	
    	Gtemp = vcat(Gtemp...) 
    	for i in 1:length(Gtemp)
    	  G[i] = Gtemp[i]
    	end
    end

   gopt!(G,pcof) = g!(G,pcof,params,order)

   res = optimize(f, gopt!, pcof0, LBFGS(),Optim.Options(show_trace =true))

   @time pcof = Optim.minimizer(res)
   display(res)

   pl1, pl2, objv, grad = objfunc.traceobjgrad(pcof,params,order, true, true, weight)

   println("Objfunc = ", objv)
   pl1


	
  	 
	
