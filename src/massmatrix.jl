
function massij(xi,zi,i,j,wi, int_polynomial)
	# to compute the mass matrix by quadrature,
	# we only need to compute the value of the
	# polynomials at point xi[j]!
	
	# xi has enough points for exactly interpolating
	# the polynomial given by the product of phi_i(x)*psi_j(x)
	
	# but phi_i(x) = 1 at xi[i] and 0 otherwise.
	# for this reason:
	return int_polynomial(xi[j],zi,i)*wi[j]
end

function massij(xi,zi,i,j,xquad, wquad, int_polynomial)
	# case where the collocation points are not Gauss-Legendre
	# an auxiliary set of quadrature points need to be evaluated
	# to compute the integral!
	mij = 0.;
	for k = 1:length(xquad)
		mij += int_polynomial(xquad[k],xi,j)*int_polynomial(xquad[k],zi,i)*wquad[k]
	end
	return mij
end

function massmatrix(N :: Integer ,a :: Number,b :: Number, order = 1)
	xi,wi = lgwt(N+order,a,b)
	if N > 1
		zi,wiz = lgwt(N,a,b)
	else
		zi,wiz = (a+b)/2, (b-a);
	end
	Mass = massmatrix(xi,zi,wi)
	Mass, wiz
end

function massmatrix(xi :: Array, zi :: Number, wi :: Array)
	massmatrix(xi, [zi], wi)
end

function massmatrix(xi :: Array, zi :: Array, wi :: Array, int_polynomial = leg_pol)
	Mass = zeros(length(zi),length(xi));
	for i = 1:length(zi)
		for j = 1:length(xi)
			Mass[i:i,j] = massij(xi,zi,i,j,wi, int_polynomial)
		end
	end
	Mass
end

function massmatrix(xi :: Array, zi :: Array, xquad :: Array, wquad :: Array, int_polynomial = leg_pol)
	Mass = zeros(length(zi),length(xi));
	for i = 1:length(zi)
		for j = 1:length(xi)
			Mass[i:i,j] = massij(xi,zi,i,j,xquad, wquad, int_polynomial)
		end
	end
	Mass
end

function Q_partial(p :: Phs, a, b)
	xi = p.disc_data.flow.xi;
	N = length(xi);
	zi,wzi = lgwt(N*2,a,b);
	Q_parts = zeros(N,N);
	for i = 1:N
		for j = 1:N
			for k = 1:length(zi)
				Q_parts[i,j] += wzi[k]*leg_pol(zi[k],xi,i)*leg_pol(zi[k],xi,j)
			end
		end
	end
	return Q_parts
end