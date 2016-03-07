function discrete_phs(N,a,b)
	Nx = N   	#    energy variables: x
	Ne = Nx+1   # co-energy variables: e
	
	ei,we = lgwt(Ne,a,b)    # discretization of co-energy variables
	#ei = collect(linspace(a,b,Ne,))
	if Nx > 1
		xi,wx = lgwt(Nx,a,b)
		#xi = collect(linspace(a,b,Nx))
	else
		xi,wx = (a+b)/2, (b-a);
	end
	#M = massmatrix(ei,xi,we);
	xquad, wquad = lgwt(Ne,a,b)
	M = massmatrix(ei,xi, xquad, wquad);
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


function discrete_phs2_distports(N,K,a,b,ad,bd)
	Nx = N
	Ne = Nx+2
	zei,we = lgwt(Ne,a,b)
	zvi,wv = lgwt(K,ad,bd)
	if Nx > 2
		zxi,wx = lgwt(Nx,a,b)
	else
		zxi,wx = (a+b)/2, (b-a);
	end
	M = massmatrix(zei,zxi,we);
	
	leg_pol_der2 = (z,zi,i) -> algo_diff(x->leg_pol(x,zi,i),2)(z);
	# M_theta_phizz can be exactly obtained by quadrature if K >= N-3:
	if K < N-3
		warn("K < N-3: The distributed ports matrix is not exactly computed by quadrature")
	end
	M_theta_phizz = massmatrix(zvi,zxi,wv, leg_pol_der2);
	M_phi = diagm(wx[:]);
	B_bar = M_phi \ M_theta_phizz;
	D = dermatrix(zei,zxi,2);
	leg_pol_der = i -> (algo_diff(x->leg_pol(x,zei,i), 1))
	p0 = map(i->leg_pol(a,zei,i), 1:length(zei))
	p0d = map(i->leg_pol_der(i)(a), 1:length(zei))
	pL = map(i->leg_pol(b,zei,i), 1:length(zei))
	pLd = map(i->leg_pol_der(i)(b), 1:length(zei))
	
	F = [M 0*M;
	     p0' pL'*0;
	     p0d' pL'*0
		 0*M M;
	     0*pLd' pLd'
		 0*p0' pL'];
	F = blkdiag(F,eye(K));
	E = [D*0 D;
	      pL'*0 p0d';
	      pL'*0  -p0';
	       -D D*0;
	       -pL' p0'*0;
	      pLd' p0'*0];
	E = blkdiag(E, zeros(K,K));
	E[1:size(D,1),end-K+1:end] = - B_bar;
    E[end-K+1:end,1:size(D,2)] =  B_bar' * M;
	J = E/F
	neworder = [1:N; (1:N)+N+2; N+1:N+2; 2*N+3:2*N+4; 2*N+4+(1:K)];
    Jnew = J[neworder,neworder]
    phJ = -Jnew[1:2*N,1:2*N];
    phB = - Jnew[1:2*N, 2*N+(1:4)];
	phBdistributed = - Jnew[1:2*N, 2*N+4+(1:K)];
    phD = Jnew[2*N+(1:4), 2*N+(1:4)];
    
	if N > 1
		Q = sparse(diagm(wx[:],0));
		Q = blkdiag(Q,Q);
		else
		Q = diagm([wx; wx][:]);
	end
	
	pflow = Collocation_points(zxi, wx)
	peffort = Collocation_points(zei, we)
	pdistributed = Collocation_points(zvi, wv)
	
	p = Phs(phJ, phB, phD, Q)
	p.Bd = phBdistributed
	#p.Q = Q
	p.disc_data = Discrete_domain(pflow, peffort)
	#func(x :: Array) = (0.5*x'*p.Q*x)[1]
	#p.Hamiltonian = func
	#p.GradHam = ForwardDiff.gradient(p.Hamiltonian)
	p
end
