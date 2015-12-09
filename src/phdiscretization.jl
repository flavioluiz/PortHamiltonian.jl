function discrete_phs(N,a,b)
	Nx = N   	#    energy variables: x
	Ne = Nx+1   # co-energy variables: e
	ei,we = lgwt(Ne,a,b)    # discretization of co-energy variables
	if Nx > 1
		xi,wx = lgwt(Nx,a,b)
	else
		xi,wx = (a+b)/2, (b-a);
	end
	M = massmatrix(ei,xi,we);
	D = dermatrix(ei,xi,1);
	p0 = map(i->leg_pol(a,ei,i), 1:length(ei))
	pL = map(i->leg_pol(b,ei,i), 1:length(ei))
	
	F = [M 0*M;
	     pL' pL'*0;
		 0*M M;
		 0*p0' p0'];
	E = -[D*0 D;
	      pL'*0 -pL';
		  D D*0;
		  p0' p0'*0];
		  J = E/F
    phJ = -[zeros(N,N) J[1:N,(N+2):end-1];
		          J[N+2:end-1,1:N] zeros(N,N)];
    phB = -[J[[1:N;N+2:2*N+1],N+1] J[[1:N;(N+2):(2*N+1)],end]];
    phD = [J[N+1,N+1] J[N+1,end];
        J[end,N+1] J[end,end]];
		
	if N > 1
		Q = sparse(diagm(wx[:],0));
		Q = blkdiag(Q,Q);
		else
		Q = diagm([wx; wx][:]);
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


function discrete_phs2(N,a,b)
	Nx = N
	Ne = Nx+2
	ei,wi = lgwt(Ne,a,b)
	if Nx > 2
		xi,wx = lgwt(Nx,a,b)
	else
		xi,wx = (a+b)/2, (b-a);
	end
	M = massmatrix(ei,xi,wi);
	D = dermatrix(ei,xi,2);
	leg_pol_der = i -> (algo_diff(x->leg_pol(x,ei,i), 1))
	p0 = map(i->leg_pol(a,ei,i), 1:length(ei))
	p0d = map(i->leg_pol_der(i)(a), 1:length(ei))
	pL = map(i->leg_pol(b,ei,i), 1:length(ei))
	pLd = map(i->leg_pol_der(i)(b), 1:length(ei))
	
	F = [M 0*M;
	     pL' pL'*0;
	     pLd' pL'*0
		 0*M M;
	     0*pLd' p0d'
		 0*p0' p0'];
	E = [D*0 D;
	      pL'*0 -pLd';
	      pL'*0  pL';
	       -D D*0;
	       p0' p0'*0;
	      -p0d' p0'*0];
	J = E/F
	neworder = [1:N; (1:N)+N+2; N+1:N+2; 2*N+3:2*N+4];
    Jnew = J[neworder,neworder]
    phJ = -Jnew[1:2*N,1:2*N];
    phB = - Jnew[1:2*N, 2*N+(1:4)];
    phD = Jnew[2*N+(1:4), 2*N+(1:4)];
    
	if N > 1
		Q = sparse(diagm(wx[:],0));
		Q = blkdiag(Q,Q);
		else
		Q = diagm([wx; wx][:]);
	end
	
	pflow = Collocation_points(xi, wx)
	peffort = Collocation_points(ei, wi)
	
	p = Phs(phJ, phB, phD, Q)
	#p.Q = Q
	p.disc_data = Discrete_domain(pflow, peffort)
	#func(x :: Array) = (0.5*x'*p.Q*x)[1]
	#p.Hamiltonian = func
	#p.GradHam = ForwardDiff.gradient(p.Hamiltonian)
	p
end
