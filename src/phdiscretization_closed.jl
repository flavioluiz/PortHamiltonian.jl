#remarks: the J matrix is skewsymmetric in two cases:
#        1) both the energy and co-energy variables are discretized
#           at Gauss-Lagrange collocation points
#        2) at least the energy variables are discretized at
#          Gauss-Lagrange... why?!
function discrete_phs_closed(N,a,b)
	Nx = N   	#    energy variables: x
	Ne = Nx   # co-energy variables: e
	
	ei,we = lgwt(Ne,a,b)    # discretization of co-energy variables
	ei1 = [a; ei; b]
	#ei1 = collect(linspace(a,b,Ne,))
	#ei = ei1[2:end-1]
	if Nx > 1
		xi,wx = lgwt(Nx,a,b)
		#xi = collect(linspace(a,b,Nx))
	else
		xi,wx = (a+b)/2, (b-a);
	end
	#M = massmatrix(ei,xi,we);
	xquad, wquad = lgwt(Ne+2,a,b)
	M2 = massmatrix(ei1,xi, xquad, wquad);
	M2 = M2[:,2:end-1]
	M1 = massmatrix(ei,xi, xquad, wquad);
	D2 = dermatrix(ei,xi,1);
	D1 = dermatrix(ei1,xi,1);
	D1 = D1[:,2:end-1]
	
    phJ = -[zeros(N,N) D1 / (M2);
		          D2 / (M1) zeros(N,N)];
    phB = zeros(2*N,1);
    phD = [];
		
	if N > 1
		Q = massmatrix(xi,xi,xquad,wquad)
		#Q = sparse(diagm(wx[:],0));
		Q = blkdiag(Q,Q);
		else
		#Q = diagm([wx; wx][:]);
	end
	
	pflow = Collocation_points(xi, wx)
	peffort = Collocation_points(ei, we)
	
	p = Phs(phJ, phB, phD, Q)
	#p.Q = Q
	p.disc_data = Discrete_domain(pflow, peffort)
	#func(x :: Array) = (0.5*x'*p.Q*x)[1]
	#p.Hamiltonian = func
	#p.GradHam = ForwardDiff.gradient(p.Hamiltonian)
	p
end