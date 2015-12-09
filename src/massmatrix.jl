
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