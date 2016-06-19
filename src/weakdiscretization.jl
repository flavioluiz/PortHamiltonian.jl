function weak_phs1(N,a,b)
	N1 = N+1   # variables x1, e1 and v1
	N2 = N   	#    variables x2, e2 and v2
	
	x1,w1,P = lglnodes(N1-1,a,b)    # discretization of x1 variables
	if N2 > 1
		x2,w2 = lgwt(N2,a,b)
	else
		x2,w2 = (a+b)/2, (b-a);
	end
	#M = massmatrix(ei,xi,we);
	#xquad, wquad = lgwt(Ne,a,b)\
	#M = massmatrix(ei,xi, xquad, wquad);
	M1 = massmatrix(x1,x1, w1);
	M2 = massmatrix([x2],[x2], [w2]);
	
	D = dermatrix(x1,x2,1)'*M2 ;
	p0 = map(i->leg_pol(a,x1,i), 1:length(x1))
	pL = map(i->leg_pol(b,x1,i), 1:length(x1))
	
	
    phJ = [zeros(N1,N1) D;
	   -D' zeros(N2, N2)];
    phB = [p0 -pL; zeros(N2,2)];
    phD = zeros(2,2);
		
	Q = blkdiag(inv(M1),inv(M2));
	#pflow = Collocation_points(x1, w1)
	#peffort = Collocation_points(x2, w2)
	
	p = Phs(phJ, phB, phD, Q)
	#p.Q = Q
	#p.disc_data = Discrete_domain(pflow, peffort)
	#func(x :: Array) = (0.5*x'*p.Q*x)[1]
	#p.Hamiltonian = func
	#p.GradHam = ForwardDiff.gradient(p.Hamiltonian)
	p
end

function weak_phs2(N,a,b)
	N1 = N   # variables x1, e1 and v1
	N2 = N+1   	#    variables x2, e2 and v2
	
	x2,w2,P = lglnodes(N2-1,a,b)    # discretization of x1 variables
	if N1 > 1
		x1,w1 = lgwt(N1,a,b)
	else
		x1,w1 = (a+b)/2, (b-a);
	end
	#M = massmatrix(ei,xi,we);
	#xquad, wquad = lgwt(Ne,a,b)\
	#M = massmatrix(ei,xi, xquad, wquad);
	M1 = massmatrix([x1],[x1], [w1]);
	M2 = massmatrix(x2,x2, w2);
	
	D = dermatrix(x2,x1,1)'*M1 ;
	p0 = map(i->leg_pol(a,x2,i), 1:length(x2))
	pL = map(i->leg_pol(b,x2,i), 1:length(x2))
	
	
    phJ = [zeros(N1,N1) -D';
	   D zeros(N2, N2)];
    phB = [zeros(N1,2); p0 -pL];
    phD = zeros(2,2);
		
		Q = blkdiag(inv(M1),inv(M2));
	
	#pflow = Collocation_points(x1, w1)
	#peffort = Collocation_points(x2, w2)
	
	p = Phs(phJ, phB, phD, Q)
	#p.Q = Q
	#p.disc_data = Discrete_domain(pflow, peffort)
	#func(x :: Array) = (0.5*x'*p.Q*x)[1]
	#p.Hamiltonian = func
	#p.GradHam = ForwardDiff.gradient(p.Hamiltonian)
	p
end

#coupled_gyrator(ph1 :: Phs, ports1, ph2 :: Phs, ports2, couple_matrix)
function finelem(Nelem,Nint,a,b)
	dx = (b-a)/Nelem;
	
	p = weak_phs1(Nint,0,dx)
	for i = 2:Nelem
		if (i%2 == 1)
			pn = weak_phs1(Nint,0,dx)
		else
			pn = weak_phs2(Nint,0,dx)
		end
		p = PortHamiltonian.coupled_gyrator(p,[2],pn,[1],[-1])
	end
	p
end