# * higher order? rdiff? ctaylor? audi?
# * moulla2
# * test cases: fluid, beam
# * Hamiltonian function -> use forwarddiff for gradient and Hessian?
# * automatic coupling phs?
# * better way to simulate?
# * automatic subdivide in elements?

module PortHamiltonian

using ForwardDiff
import Base.show
import Base.eig

include("types.jl")
include("lgwt.jl")
include("derivative.jl")
include("massmatrix.jl")
include("phdiscretization.jl")
include("eigenphs.jl")

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
				#println((x[1]-xi[k])/(xi[j]-xi[k]))
            end
        end
    return P
end

end
