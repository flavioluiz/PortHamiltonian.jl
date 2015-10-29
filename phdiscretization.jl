# * higher order? rdiff? ctaylor? audi?
# * moulla2
# * test cases: fluid, beam
# * Hamiltonian function -> use forwarddiff for gradient and Hessian?
# * automatic coupling phs?
# * better way to simulate?
# * automatic subdivide in elements?


include("lgwt.jl")
include("derivative.jl")
include("massmatrix.jl")
include("../ForwardDiff.jl/src/ForwardDiff.jl")
function leg_pol(x :: Array,xi :: Array,j :: Integer)
	leg_pol(x[1],xi,j)
end
function leg_pol(x :: Number,xi :: Array,j :: Integer)
	# evaluate the "j" Lagrange polynomial defined using
    # the collocation points "xi" at point "x"
	P = 1.;
        for k = 1:length(xi)
	        if (k != j)
		        P = P*(x[1]-xi[k])/(xi[j]-xi[k])
            end
        end
    return P
end

type Phs
	J;
	B;
	D;
	wi;
    zi;
	Q;
end

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
		Q = sparse(diagm(wi[:],0));
		Q = blkdiag(Q,Q);
		else
		Q = diagm([wi; wi][:]);
	end
	
	return Phs(phJ, phB, phD, phw,phzi,Q)
end

