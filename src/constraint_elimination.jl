function constraint_elimination(p :: Phs)
# Given a constrained Phs with N state variables
#     and Nconst constraints lambda:
#                          Xdot = J GradH + B u + G lambda
#                          y = B' GradH + D u
#                          0 = G' GradH
# this function finds an equivalent Phs with reduced order (N-Nconst)
# without the constraints lambda. In the linear case, it is equivalent
# to:
#                  z1dot = Jnew Qnew z1 + Bnew u
#                  y = Bnew' Qnew z1 + D u
# 
	if isdefined(p,:Q)
		nstates = size(p.Q,1)
		nconst = size(p.G,2)
		Gnull = nullspace(p.G')'
		M = [Gnull; inv((p.G'*p.G)[1])*p.G']
		Jtilde = M*p.J*M'
		Minv = inv(M)
		Qtilde = Minv'*p.Q*Minv
		J11 = Jtilde[1:nstates-nconst,1:nstates-nconst]
		Q11 = Qtilde[1:nstates-nconst,1:nstates-nconst]
		Q12 = Qtilde[1:nstates-nconst,nstates+1:end]
		Q21 = Qtilde[nstates+1:end,1:nstates-nconst]
		Q22 = Qtilde[nstates+1:end,nstates+1:end]
		Qn = Q11-Q12*inv(Q22)*Q21
		Phs(J11, Gnull * p.B, p.D, Qn)
	else
	
		println("Nonlinear system. Not working yet.")
	end
end