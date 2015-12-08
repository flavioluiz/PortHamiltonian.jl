# * higher order? rdiff? ctaylor? audi?
# * moulla2
# * test cases: fluid, beam
# * Hamiltonian function -> use forwarddiff for gradient and Hessian?
# * automatic coupling phs?
# * better way to simulate?
# * automatic subdivide in elements?
__precompile__()

module PortHamiltonian

using ForwardDiff


###########
# imports # 
########### 
# additional functionality is added to the following Base functions:
import Base: show, eig, blkdiag

###########
# exports #
###########
export set_constraint!,
	   frequencies,
	   constraint_elimination,
	   coupled_gyrator,
	   discrete_phs

function blkdiag(x :: Array, y :: Array)
	res = zeros(size(x,1)+size(y,1),size(x,2)+size(y,2));
	res[1:size(x,1), 1:size(x,2)] = x;
	res[1+size(x,1):end, 1+size(x,2):end] = y;
	return res;
end

include("types.jl")
include("lgwt.jl")
include("derivative.jl")
include("massmatrix.jl")
include("phdiscretization.jl")
include("eigenphs.jl")
include("constraint_elimination.jl")
include("spectral_element_phs.jl")
include("coupled_phs.jl")

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
