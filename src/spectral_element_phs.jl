function new_spectral_element_phs(Nel,Npol,a,b)
	dx = (b-a)/Nel
	ph = discrete_phs(Npol,0.,dx);
	zi = ph.disc_data.flow.xi;
	wi = ph.disc_data.flow.wi;
	zvec = zeros(length(zi)*Nel)
	wvec = zeros(length(wi)*Nel)
	
	phfull = ph
	neworder = (Nel-1)*Npol*2 + collect(1:Npol)
	for e = 2:Nel
		phfull = PortHamiltonian.coupled_gyrator(phfull,[2],ph,[1],1)
		for ii = (Nel-e)*Npol*2 + collect(1:Npol)
			push!(neworder,ii)
		end
	end
	neworder = [neworder;neworder+Npol]
	phfull.J = phfull.J[neworder,neworder]
	phfull.B = phfull.B[neworder,:]
	return phfull

end

function new_spectral_element_phs_golo(Nel,alfa,a,b)
	dx = (b-a)/Nel
	ph = Phs([0 -1/alfa;1/alfa 0], [1/alfa 0;0. -1/alfa], [0 -(1-alfa)/alfa;-(alfa-1)/alfa 0], eye(2)/dx)
	Npol = 1
	phfull = ph
	neworder = (Nel-1)*Npol*2 + collect(1:Npol)
	for e = 2:Nel
		phfull = PortHamiltonian.coupled_gyrator(phfull,[2],ph,[1],1)
		for ii = (Nel-e)*Npol*2 + collect(1:Npol)
			push!(neworder,ii)
		end
	end
	neworder = [neworder;neworder+Npol]
	phfull.J = phfull.J[neworder,neworder]
	phfull.B = phfull.B[neworder,:]
	return phfull

end


function new_block_spectral_element_phs(Nel,Npol,a,b)
	dx = (b-a)/Nel
	ph = discrete_phs(Npol,0.,dx);
	zi = ph.disc_data.flow.xi;
	wi = ph.disc_data.flow.wi;
	zvec = zeros(length(zi)*Nel)
	wvec = zeros(length(wi)*Nel)
	blockph = Phs([0 -1;1 0], [1 0;0. -1], zeros(2,2), eye(2)*1e5)
	
	phfull = ph
	neworder = (Nel-1)*Npol*2 + collect(1:Npol)
	for e = 2:Nel
		phfull = PortHamiltonian.coupled_gyrator(phfull,[2],blockph,[1],1)
		phfull = PortHamiltonian.coupled_gyrator(phfull,[2],ph,[1],1)
		for ii = (Nel-e)*Npol*2 + collect(1:Npol)
			push!(neworder,ii)
		end
	end
	#neworder = [neworder;neworder+Npol]
	#phfull.J = phfull.J[neworder,neworder]
	#phfull.B = phfull.B[neworder,:]
	return phfull

end


type spectral_element_phs
	phs::Phs
	
	
	function spectral_element_phs(Nel,Npol,a,b)
	# creates a "spectral element" phs model for
	# first order derivative PHS
	# Npol: is the polynomial approximation order
	# Nel: is the number of elements
	# (a,b): is the domain interval
	
	# to do: 
	# * matrices B and D
	# * create a Hamiltonian/grad method or something like that
		this = new();
		dx = (b-a)/Nel;
		xpoints = linspace(a,b,Nel+1);
		ph = discrete_phs(Npol,0.,dx);
		zi = ph.disc_data.flow.xi;
		wi = ph.disc_data.flow.wi;
		zvec = zeros(length(zi)*Nel)
		wvec = zeros(length(wi)*Nel)
		J = zeros(2*length(wi)*Nel, 2*length(wi)*Nel)
		alpha1_index = Array(Int, length(wi)*Nel,1)
		alpha2_index = Array(Int, length(wi)*Nel,1)
		
		for e = 1:Nel
			ze = zi + xpoints[e]
			zvec[1+length(zi)*(e-1):length(zi)*e] = ze;
			println(wi, 1+length(zi)*(e-1):length(zi)*e)
			wvec[1+length(zi)*(e-1):length(zi)*e] = wi;
			J[1+2*length(zi)*(e-1):2*length(zi)*e,1+2*length(zi)*(e-1):2*length(zi)*e] = ph.J;			
			for ei = (e+1):Nel
				J[1+2*length(zi)*(e-1):2*length(zi)*e,1+2*length(zi)*(ei-1):2*length(zi)*(ei)] = -(-ph.D[2,1])^(ei-e+1)*ph.B[:,1]*ph.B[:,2]';
			end
			for ei = 1:(e-1)
				J[1+2*length(zi)*(e-1):2*length(zi)*e,1+2*length(zi)*(ei-1):2*length(zi)*(ei)] = (ph.D[1,2])^(ei-e+1)*ph.B[:,2]*ph.B[:,1]';
			end
			alpha1_index[1+length(zi)*(e-1):length(zi)*e] = [(1:length(zi))+length(zi)*2*(e-1)];
			alpha2_index[1+length(zi)*(e-1):length(zi)*e] = (1:length(zi))+length(zi)*(2*e-1);
		end
		
		a = alpha1_index;
		b = alpha2_index;
		Jn = J[[a[:];b[:]],[a[:];b[:]]]
		Qn = sparse(diagm(wvec[:],0));
		Qn = blkdiag(Qn,Qn);
		return Jn,Qn,a,b,ph.D
	end
end