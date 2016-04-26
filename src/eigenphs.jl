
#function blkd(A :: Array,B :: Array)
#	# block diagonal matrix of A and B
#	al = size(A,1)
#	ac = size(A,2)
#	bl = size(B,1)
#	bc = size(B,2)
#	C = zeros(al+bl, ac + bc)
#	C[1:al, 1:ac] = A
#	C[al+1:end, ac+1:end] = B
#	C
#end
function eig(p :: Phs)
	# computes the eigenvalues and eigenvectors of a linear PHS
	#  (which haves a Q matrix!)
	#
	#	constrained (G != 0) cases can be addressed using
	#  generalized eigenvalues
	#
	#
	#  TODO: eigenvalues around a linearized condition:
	#       just need to compute the Hamiltonian Hessian at a given point!
	#       For e.g.: hess = ForwardDiff.hessian(p.Hamiltonian)
	#                 Q = hess(x)
	#               --> this can be useful for validating our nonlinear system!
	# 

	if isdefined(p, :Q)
		if isdefined(p,:G)
			nconst = size(p.G,2)
			I = eye(size(p.J,1))
			z = zeros(nconst,nconst)
			E = blkdiag(I,z)
			A = [p.J*p.Q p.G;
			     p.G'*p.Q p.G_D]
			e = eigfact(A,E)
			ind = sortperm(imag(e.values))
			return e.values[ind], e.vectors[ind,ind]
		else
			a,v = eig(p.J*p.Q)
			ind = sortperm(imag(a))
			return a[ind], v[:,ind]
		end
	else
		fprintln("Undefined Q matrix (nonlinear system?)!")
		return
	end
end

function eigdamp(p :: Phs)
	# computes the frequency and damping of PHS with damping
	#  (which haves a Q matrix!)
	#
	#	constrained (G != 0) cases can be addressed using
	#  generalized eigenvalues
	#

	if isdefined(p, :Q)
		if isdefined(p,:G)
			nconst = size(p.G,2)
			I = eye(size(p.J,1))
			z = zeros(nconst,nconst)
			E = blkdiag(I,z)
			A = [(p.J-p.R)*p.Q p.G;
			     p.G'*p.Q p.G_D]
			e = eigfact(A,E)
			ind = sortperm(imag(e.values))
			return e.values[ind], e.vectors[ind,ind]
		else
			a,v = eig((p.J-p.R)*p.Q)
			ind = sortperm(imag(a))
			return a[ind], v[:,ind]
		end
	else
		fprintln("Undefined Q matrix (nonlinear system?)!")
		return
	end
end

function damp(p :: Phs)
	a,v = eigdamp(p)
	return [abs(a) -cos(angle(a))]
end
function frequencies(a :: Array)
		imag(a[imag(a).>=0])
end

function frequencies(ph :: Phs)
	a,v = eig(ph)
	return frequencies(a)
end