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
	   removeport!,
	   frequencies,
	   constraint_elimination,
	   coupled_gyrator,
	   discrete_phs,
	   discrete_phs_closed,
	   discrete_phs2,
	   discrete_phs2_distports,
	   finelem

function blkdiag(x :: Array, y :: Array)
	res = zeros(size(x,1)+size(y,1),size(x,2)+size(y,2));
	res[1:size(x,1), 1:size(x,2)] = x;
	res[1+size(x,1):end, 1+size(x,2):end] = y;
	return res;
end

include("types.jl")
function blkdiag(x :: Phs, y :: Phs)
	warn("only linear phs with constraint can be concatenated yet")
	phnew = Phs(blkdiag(x.J,y.J),blkdiag(x.B,y.B),blkdiag(x.D,y.D), blkdiag(x.Q,y.Q))
	phnew.TransfMatrix = blkdiag(x.TransfMatrix,y.TransfMatrix)
	phnew.R = blkdiag(x.R,y.R)	
	Nx = size(x.TransfMatrix,1)
	for key in keys(y.StatesNames)
		y.StatesNames[key] += Nx
	end
	phnew.StatesNames = merge(x.StatesNames, y.StatesNames)
	phnew.InputsNames = [x.InputsNames, y.InputsNames]
	return phnew
end

function removeport!(p :: Phs, portnumber :: Integer)
	keep = 1:size(p.B,2) .!= portnumber
	p.B = p.B[:,keep]
	p.D = p.D[keep,keep]
	return p
end
include("lgwt.jl")
include("lglnodes.jl")
include("derivative.jl")
include("massmatrix.jl")
include("phdiscretization.jl")
include("phdiscretization_closed.jl")
include("eigenphs.jl")
include("constraint_elimination.jl")
include("spectral_element_phs.jl")
include("coupled_phs.jl")
include("weakdiscretization.jl")

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
