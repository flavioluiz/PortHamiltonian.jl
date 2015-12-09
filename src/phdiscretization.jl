function discrete_phs(N,a,b)
	Nz = N
	Nx = Nz+1
	xi,wi = lgwt(Nx,a,b)
	if Nz > 1
		zi,wiz = lgwt(Nz,a,b)
	else
		zi,wiz = (a+b)/2, (b-a);
	end
	M = massmatrix(xi,zi,wi);
	D = dermatrix(xi,zi,1);
	p0 = map(i->leg_pol(a,xi,i), 1:length(xi))
	pL = map(i->leg_pol(b,xi,i), 1:length(xi))
	
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
	phw = wiz;
	phzi = zi;
	if N > 1
		Q = sparse(diagm(wiz[:],0));
		Q = blkdiag(Q,Q);
		else
		Q = diagm([wiz; wiz][:]);
	end
	
	pflow = Collocation_points(zi, wiz)
	peffort = Collocation_points(xi, wi)
	
	p = Phs(phJ, phB, phD, Q)
	#p.Q = Q
	p.disc_data = Discrete_domain(pflow, peffort)
	#func(x :: Array) = (0.5*x'*p.Q*x)[1]
	#p.Hamiltonian = func
	#p.GradHam = ForwardDiff.gradient(p.Hamiltonian)
	p
end


function discrete_phs2(N,a,b)
	Nz = N
	Nx = Nz+2
	xi,wi = lgwt(Nx,a,b)
	if Nz > 2
		zi,wiz = lgwt(Nz,a,b)
	else
		zi,wiz = (a+b)/2, (b-a);
	end
	M = massmatrix(xi,zi,wi);
	D = dermatrix(xi,zi,2);
	leg_pol_der = i -> (algo_diff(x->leg_pol(x,xi,i), 1))
	p0 = map(i->leg_pol(a,xi,i), 1:length(xi))
	p0d = map(i->leg_pol_der(i)(a), 1:length(xi))
	pL = map(i->leg_pol(b,xi,i), 1:length(xi))
	pLd = map(i->leg_pol_der(i)(b), 1:length(xi))
	
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
		  print (norm(J+J'))
	neworder = [1:N, (1:N)+N+2, N+1:N+2, 2*N+3:2*N+4];
    Jnew = J[neworder,neworder]
    phJ = -Jnew[1:2*N,1:2*N];
    phB = - Jnew[1:2*N, 2*N+(1:4)];
    phD = Jnew[2*N+(1:4), 2*N+(1:4)];
    
	phw = wiz;
	phzi = zi;
	if N > 1
		Q = sparse(diagm(wiz[:],0));
		Q = blkdiag(Q,Q);
		else
		Q = diagm([wiz; wiz][:]);
	end
	
	pflow = Collocation_points(zi, wiz)
	peffort = Collocation_points(xi, wi)
	
	p = Phs(phJ, phB, phD, Q)
	#p.Q = Q
	p.disc_data = Discrete_domain(pflow, peffort)
	#func(x :: Array) = (0.5*x'*p.Q*x)[1]
	#p.Hamiltonian = func
	#p.GradHam = ForwardDiff.gradient(p.Hamiltonian)
	p
end
