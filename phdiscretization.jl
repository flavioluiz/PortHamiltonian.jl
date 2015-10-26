include("lgwt.jl")
include("../ForwardDiff.jl/src/ForwardDiff.jl")
function leg_pol(x,xi,j)
	# evaluate the "j" Lagrange polynomial at point "x"
	P = 1.;
        for k = 1:length(xi)
	        if (k != j)
		        P = P*(x[1]-xi[k])/(xi[j]-xi[k])
            end
        end
    return P
end

function pol_vec(x,xi)
	# evaluate each Lagrange polynomial at point 'x'
	P = zeros(length(xi),1);
	for i = 1:length(xi)
		P[i] = leg_pol(x,xi,i)
	end
	return P
end

function massij(xi,zi,i,j,wi)
	# to compute the mass matrix by quadrature,
	# we only need to compute the value of the
	# polynomials at point xi[j]!
	
	# xi has enough points for exactly interpolating
	# the polynomial given by the product of phi_i(x)*psi_j(x)
	
	# but phi_i(x) = 1 at xi[i] and 0 otherwise.
	# for this reason:
	return leg_pol(xi[j],zi,i)*wi[j]
end

function massmatrix(N,a,b)
	xi,wi = lgwt(N+1,a,b)
	zi,wiz = lgwt(N,a,b)
	Mass = zeros(length(zi),length(xi));
	for i = 1:length(zi)
		for j = 1:length(xi)
			Mass[i:i,j] = massij(xi,zi,i,j,wi)
		end
	end
	return Mass
end

type Phs
	J;
	B;
	D;
end

function moulla(N,a,b)
	M = massmatrix(N,a,b);
	D,xi,zi = dermatrix(N,a,b);
	p0 = pol_vec(a,xi);
	pL = pol_vec(b,xi);
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

	return Phs(phJ, phB, phD)
end

function dermatrix(N,a,b)
	#xi = cos(pi*(1:1:(N+1))/(N+2))
	#zi = cos(pi*(1:1:(N))/(N+1))
	xi,wi = lgwt(N+1,a,b)
	zi,wiz = lgwt(N,a,b)
	D = zeros(length(zi),length(xi))
	for j = 1:length(xi)
		g = ForwardDiff.derivative(x-> leg_pol(x,xi,j))
		for i = 1:length(zi)		
			D[i:i,j] = g(zi[i])
		end
	end
	return D,xi,zi
end

function dermatrix2(N)
	xi = cos(pi*(1:1:(N+2))/(N+3))
	zi = cos(pi*(1:1:(N))/(N+1))
	D = zeros(length(zi),length(xi))
	for j = 1:length(xi)
		g = ForwardDiff.hessian(x-> leg_pol(x,xi,j))
		for i = 1:length(zi)		
			D[i:i,j] = g([zi[i]])
		end
	end
	return D,xi,zi
end
