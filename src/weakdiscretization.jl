function weak_phs_FEM(Ne,N1, N2,a,b, weightQ1, weightQ2)
	# this function uses the partioned FEM method, finding a single PHS model
	# for the whole (1D) domain.
	# TO DO: - correct the input matrix, and allow for higher-order polynomials
	#	N1  (old = Norder+1)  variables x1, e1 and v1
	#	N2  (old = Norder)  #    variables x2, e2 and v2
		
		x1,w1,P = lglnodes(N1-1,(b-a)/Ne,0)    # discretization of x1 variables
		x1i,w1i,Pi = lglnodes(N1,0,(b-a)/Ne)  # these points are used to integrate the mass matrix using quadrature
		if N2 > 1
			x2,w2 = lglnodes(N2-1,(b-a)/Ne,0)
		else
			x2,w2 = [(a+b)/2], [(b-a)/Ne];
		end
		#M = massmatrix(ei,xi,we);
		#xquad, wquad = lgwt(Ne,a,b)\
		#M = massmatrix(ei,xi, xquad, wquad);
		M1 = massmatrix(x1,x1, x1i,w1i);
		x1i2,w1i2,Pi2 = lglnodes(N1*4,0,(b-a)/Ne)  # these points are used to integrate the mass matrix using quadrature
		L1 = ELid -> massmatrix(x1,x1, x1i2,w1i2, (x, zi, i) -> PortHamiltonian.leg_pol(x,zi,i)*weightQ1(x+(ELid-1)*(b-a)/Ne));
		
#		if N2 == 1
#			M2 = massmatrix(x2,x2, w2);
#		else
		x2i,w2i,P2i = lglnodes(N2,0,(b-a)/Ne)  # these points are used to integrate the mass matrix using quadrature
		M2 = massmatrix(x2,x2, x2i,w2i);
		x2i2,w2i2,P2i2 = lglnodes(N2*4,0,(b-a)/Ne)
		L2 = ELid -> massmatrix(x2,x2, x2i2,w2i2, (x, zi, i) -> PortHamiltonian.leg_pol(x,zi,i)*weightQ2(x+(ELid-1)*(b-a)/Ne));
#		end
		

		D = dermatrix(x1,x2,1)'*M2 ;		
		
		p0 = map(i->leg_pol(0,x1,i), 1:length(x1))
		pL = map(i->leg_pol((b-a)/Ne,x1,i), 1:length(x1))
		B = [p0 pL]
	
		M1full = zeros(Ne*(N1-1)+1,Ne*(N1-1)+1);
		L1full = zeros(Ne*(N1-1)+1,Ne*(N1-1)+1);
		if N2 == 1
			M2full = zeros(Ne,Ne);
			L2full = zeros(Ne,Ne);
		else
			M2full = zeros((Ne)*(N2-1)+1,(Ne)*(N2-1)+1);
			L2full = zeros((Ne)*(N2-1)+1,(Ne)*(N2-1)+1);
		end
		if N2 == 1
			Dfull = zeros(Ne*(N1-1)+1,Ne);
		else
			Dfull = zeros(Ne*(N1-1)+1,(Ne)*(N2-1)+1);
		end
		Bfull = zeros(Ne*(N1-1)+1,2);
	
		for i = 1:Ne
			M1full[((i-1)*(N1-1)+1):((i)*(N1-1)+1),((i-1)*(N1-1)+1):((i)*(N1-1)+1)] = M1full[((i-1)*(N1-1)+1):((i)*(N1-1)+1),((i-1)*(N1-1)+1):((i)*(N1-1)+1)] + M1;
			L1full[((i-1)*(N1-1)+1):((i)*(N1-1)+1),((i-1)*(N1-1)+1):((i)*(N1-1)+1)] = L1full[((i-1)*(N1-1)+1):((i)*(N1-1)+1),((i-1)*(N1-1)+1):((i)*(N1-1)+1)] + L1(i);
			if N2 == 1
				M2full[i:i,i] = M2;
				L2full[i:i,i] = L2(i);
			else
				M2full[((i-1)*(N2-1)+1):((i)*(N2-1)+1),((i-1)*(N2-1)+1):((i)*(N2-1)+1)] = M2full[((i-1)*(N2-1)+1):((i)*(N2-1)+1),((i-1)*(N2-1)+1):((i)*(N2-1)+1)] + M2;
				L2full[((i-1)*(N2-1)+1):((i)*(N2-1)+1),((i-1)*(N2-1)+1):((i)*(N2-1)+1)] = L2full[((i-1)*(N2-1)+1):((i)*(N2-1)+1),((i-1)*(N2-1)+1):((i)*(N2-1)+1)] + L2(i);
			end
			if N2 == 1
				Dfull[((i-1)*(N1-1)+1):((i)*(N1-1)+1),i] = D;
			else
				Dfull[((i-1)*(N1-1)+1):((i)*(N1-1)+1),((i-1)*(N2-1)+1):((i)*(N2-1)+1)] = Dfull[((i-1)*(N1-1)+1):((i)*(N1-1)+1),((i-1)*(N2-1)+1):((i)*(N2-1)+1)]+D;
			end
			if i == 1
				Bfull[((i-1)*(N1-1)+1):((i)*(N1-1)+1),1] =Bfull[((i-1)*(N1-1)+1):((i)*(N1-1)+1),1]+B[:,1];
			end
			if i == Ne
				Bfull[((i-1)*(N1-1)+1):((i)*(N1-1)+1),2] =Bfull[((i-1)*(N1-1)+1):((i)*(N1-1)+1),2]+B[:,2];
			end
		end
	
		# interconnection matrix J:
		if N2 == 1
		phJ = [zeros(Ne*(N1-1)+1,Ne*(N1-1)+1) Dfull;
			-Dfull' zeros(Ne,Ne)];
		else
			phJ = [zeros(Ne*(N1-1)+1,Ne*(N1-1)+1) Dfull;
			-Dfull' zeros((Ne)*(N2-1)+1,(Ne)*(N2-1)+1)];
		end
		phD = zeros(2,2)
		# Hamiltonian matrix Q:
		Q = blkdiag((M1full\L1full)/(M1full), (M2full\L2full)/(M2full));		
	
		# input matrix
		if N2 == 1
			phB = [Bfull; zeros(Ne,2)];
		else
			phB = [Bfull; zeros((Ne)*(N2-1)+1,2)];
		end
	
	
		p = Phs(phJ, phB, phD, Q)
	
		p
end


function weak_phs_FEM(Ne,N1, N2,a,b)
	# this function uses the partioned FEM method, finding a single PHS model
	# for the whole (1D) domain.
	# TO DO: - correct the input matrix, and allow for higher-order polynomials
	#	N1  (old = Norder+1)  variables x1, e1 and v1
	#	N2  (old = Norder)  #    variables x2, e2 and v2
		
		x1,w1,P = lglnodes(N1-1,(b-a)/Ne,0)    # discretization of x1 variables
		x1i,w1i,Pi = lglnodes(N1,0,(b-a)/Ne)  # these points are used to integrate the mass matrix using quadrature
		if N2 > 1
			x2,w2 = lglnodes(N2-1,(b-a)/Ne,0)
		else
			x2,w2 = [(a+b)/2], [(b-a)/Ne];
		end
		#M = massmatrix(ei,xi,we);
		#xquad, wquad = lgwt(Ne,a,b)\
		#M = massmatrix(ei,xi, xquad, wquad);
		M1 = massmatrix(x1,x1, x1i,w1i);
		if N2 == 1
			M2 = massmatrix(x2,x2, w2);
		else
			x2i,w2i,P2i = lglnodes(N2,0,(b-a)/Ne)  # these points are used to integrate the mass matrix using quadrature
			M2 = massmatrix(x2,x2, x2i,w2i);
		end
		
		D = transpose(dermatrix(x1,x2,1))*M2 ;		
		
		p0 = map(i->leg_pol(0,x1,i), 1:length(x1))
		pL = map(i->leg_pol((b-a)/Ne,x1,i), 1:length(x1))
		B = [p0 pL]
	
		M1full = zeros(Ne*(N1-1)+1,Ne*(N1-1)+1);
		if N2 == 1
			M2full = zeros(Ne,Ne);
		else
			M2full = zeros((Ne)*(N2-1)+1,(Ne)*(N2-1)+1);
		end
		if N2 == 1
			Dfull = zeros(Ne*(N1-1)+1,Ne);
		else
			Dfull = zeros(Ne*(N1-1)+1,(Ne)*(N2-1)+1);
		end
		Bfull = zeros(Ne*(N1-1)+1,2);
	
		for i = 1:Ne
			M1full[((i-1)*(N1-1)+1):((i)*(N1-1)+1),((i-1)*(N1-1)+1):((i)*(N1-1)+1)] = M1full[((i-1)*(N1-1)+1):((i)*(N1-1)+1),((i-1)*(N1-1)+1):((i)*(N1-1)+1)] + M1;
			if N2 == 1
				M2full[i:i,i] = M2;
			else
				M2full[((i-1)*(N2-1)+1):((i)*(N2-1)+1),((i-1)*(N2-1)+1):((i)*(N2-1)+1)] = M2full[((i-1)*(N2-1)+1):((i)*(N2-1)+1),((i-1)*(N2-1)+1):((i)*(N2-1)+1)] + M2;
			end
			if N2 == 1
				Dfull[((i-1)*(N1-1)+1):((i)*(N1-1)+1),i] = D;
			else
				Dfull[((i-1)*(N1-1)+1):((i)*(N1-1)+1),((i-1)*(N2-1)+1):((i)*(N2-1)+1)] = Dfull[((i-1)*(N1-1)+1):((i)*(N1-1)+1),((i-1)*(N2-1)+1):((i)*(N2-1)+1)]+D;
			end
			if i == 1
				Bfull[((i-1)*(N1-1)+1):((i)*(N1-1)+1),1] =Bfull[((i-1)*(N1-1)+1):((i)*(N1-1)+1),1]+B[:,1];
			end
			if i == Ne
				Bfull[((i-1)*(N1-1)+1):((i)*(N1-1)+1),2] =Bfull[((i-1)*(N1-1)+1):((i)*(N1-1)+1),2]+B[:,2];
			end
		end
	
		# interconnection matrix J:
		if N2 == 1
		phJ = [zeros(Ne*(N1-1)+1,Ne*(N1-1)+1) Dfull;
			-transpose(Dfull) zeros(Ne,Ne)];
		else
			phJ = [zeros(Ne*(N1-1)+1,Ne*(N1-1)+1) Dfull;
			-transpose(Dfull) zeros((Ne)*(N2-1)+1,(Ne)*(N2-1)+1)];
		end
		phD = zeros(2,2)
		# Hamiltonian matrix Q:
		Q = blkdiag(inv(M1full), inv(M2full));
	
	
		# input matrix
		if N2 == 1
			phB = [Bfull; zeros(Ne,2)];
		else
			phB = [Bfull; zeros((Ne)*(N2-1)+1,2)];
		end
	
	
		p = Phs(phJ, phB, phD, Q)
	
		p
	end
	
function weak_phs_FEM(Ne,Norder,a,b)
	weak_phs_FEM(Ne,Norder+1,Norder+1,a,b)
end

#= 
function weak_phs_FEM(Ne,Norder,a,b)
# this function uses the partioned FEM method, finding a single PHS model
# for the whole (1D) domain.
# TO DO: - correct the input matrix, and allow for higher-order polynomials
	N1 = Norder+1   # variables x1, e1 and v1
	N2 = Norder   	#    variables x2, e2 and v2
	
	x1,w1,P = lglnodes(N1-1,(b-a)/Ne,0)    # discretization of x1 variables
	x1i,w1i,Pi = lglnodes(N1,0,(b-a)/Ne)  # these points are used to integrate the mass matrix using quadrature
	if N2 > 1
		x2,w2 = lglnodes(N2-1,(b-a)/Ne,0)
	else
		x2,w2 = [(a+b)/2], [(b-a)/Ne];
	end
	#M = massmatrix(ei,xi,we);
	#xquad, wquad = lgwt(Ne,a,b)\
	#M = massmatrix(ei,xi, xquad, wquad);
	M1 = massmatrix(x1,x1, x1i,w1i);
	if Norder == 1
		M2 = massmatrix(x2,x2, w2);
	else
		x2i,w2i,P2i = lglnodes(N2,0,(b-a)/Ne)  # these points are used to integrate the mass matrix using quadrature
		M2 = massmatrix(x2,x2, x2i,w2i);
	end
	
	D = dermatrix(x1,x2,1)'*M2 ;
	
	p0 = map(i->leg_pol(0,x1,i), 1:length(x1))
	pL = map(i->leg_pol((b-a)/Ne,x1,i), 1:length(x1))
	B = [p0 pL]
	display(B)

	M1full = zeros(Ne*Norder+1,Ne*Norder+1);
	if Norder == 1
		M2full = zeros(Ne,Ne);
	else
		M2full = zeros((Ne)*(Norder-1)+1,(Ne)*(Norder-1)+1);
	end
	if Norder == 1
		Dfull = zeros(Ne*Norder+1,Ne);
	else
		Dfull = zeros(Ne*Norder+1,(Ne)*(Norder-1)+1);
	end
	Bfull = zeros(Ne*Norder+1,2);

	for i = 1:Ne
		M1full[((i-1)*Norder+1):((i)*Norder+1),((i-1)*Norder+1):((i)*Norder+1)] = M1full[((i-1)*Norder+1):((i)*Norder+1),((i-1)*Norder+1):((i)*Norder+1)] + M1;
		if Norder == 1
			M2full[i:i,i] = M2;
		else
			M2full[((i-1)*(Norder-1)+1):((i)*(Norder-1)+1),((i-1)*(Norder-1)+1):((i)*(Norder-1)+1)] = M2full[((i-1)*(Norder-1)+1):((i)*(Norder-1)+1),((i-1)*(Norder-1)+1):((i)*(Norder-1)+1)] + M2;
		end
		if Norder == 1
			Dfull[((i-1)*Norder+1):((i)*Norder+1),i] = D;
		else
			Dfull[((i-1)*Norder+1):((i)*Norder+1),((i-1)*(Norder-1)+1):((i)*(Norder-1)+1)] = Dfull[((i-1)*Norder+1):((i)*Norder+1),((i-1)*(Norder-1)+1):((i)*(Norder-1)+1)]+D;
		end
		if i == 1
			Bfull[((i-1)*Norder+1):((i)*Norder+1),1] =Bfull[((i-1)*Norder+1):((i)*Norder+1),1]+B[:,1];
		end
		if i == Ne
			Bfull[((i-1)*Norder+1):((i)*Norder+1),2] =Bfull[((i-1)*Norder+1):((i)*Norder+1),2]+B[:,2];
		end
	end

	# interconnection matrix J:
	if Norder == 1
	phJ = [zeros(Ne+1,Ne+1) Dfull;
		-Dfull' zeros(Ne,Ne)];
	else
		phJ = [zeros(Ne*Norder+1,Ne*Norder+1) Dfull;
		-Dfull' zeros((Ne)*(Norder-1)+1,(Ne)*(Norder-1)+1)];
	end
	phD = zeros(2,2)
	# Hamiltonian matrix Q:
	Q = blkdiag(inv(M1full), inv(M2full));


	# input matrix
	if Norder == 1
		phB = [Bfull; zeros(Ne,2)];
	else
		phB = [Bfull; zeros((Ne)*(Norder-1)+1,2)];
	end


	p = Phs(phJ, phB, phD, Q)

	p
end =#

function weak_phs1(N,a,b)
# this function uses the partioned FEM method, finding a PHS model for a
# part of the domain. This could be used with higher order polynomials
# (using a single element), or coupling several elements (using lower-order polynomials)
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
	M2 = massmatrix(x2,x2, w2);
	
	D = transpose(dermatrix(x1,x2,1))*M2 ;
	p0 = map(i->leg_pol(a,x1,i), 1:length(x1))
	pL = map(i->leg_pol(b,x1,i), 1:length(x1))
	
	
    phJ = [zeros(N1,N1) D;
	   -transpose(D) zeros(N2, N2)];
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
	M1 = massmatrix(x1,x1, w1);
	M2 = massmatrix(x2,x2, w2);
	
	D = transpose(dermatrix(x2,x1,1))*M1 ;
	p0 = map(i->leg_pol(a,x2,i), 1:length(x2))
	pL = map(i->leg_pol(b,x2,i), 1:length(x2))
	
	
    phJ = [zeros(N1,N1) -transpose(D);
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