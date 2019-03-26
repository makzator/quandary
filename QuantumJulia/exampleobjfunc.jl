function exampleobjfunc()
	N = 4
	#3
	Nguard = 3
	Ntot = N + Nguard
	
	Ident = Matrix{Float64}(I, Ntot, Ntot)   
	utarget = Ident[1:Ntot,1:N]
	utarget[:,3] = Ident[:,4]
	utarget[:,4] = Ident[:,3]
	
	#utarget[:,1] = Ident[:,2]
	#utarget[:,2] = Ident[:,1]
	
	cfl = 0.05
	T = 0.2
	testadjoint = 0
	maxpar =0.09
	
	params = objfunc.parameters(N,Nguard,T,testadjoint,maxpar,cfl, utarget)
	#pcof = rand(4)
	pcof = [0.001, 0.0, 0.0]
	order = 2
	
	#pl1, pl2, objv, grad = objfunc.traceobjgrad(pcof,params,order, true)

    objv, grad = objfunc.traceobjgrad(pcof, params, order, false, true)
	
	println("objv: ", objv)
	println("objgrad: ", grad)
	
#	pl1
end