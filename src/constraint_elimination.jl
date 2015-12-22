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
		if (p.G_D .< eps()) == trues(size(p.G_D))
			nstates = size(p.Q,1)
			nconst = size(p.G,2)
			Gnull = nullspace(p.G')'
			M = [Gnull; inv((p.G'*p.G)[1])*p.G']
			Jtilde = M*p.J*M'
			Minv = inv(M)
			Qtilde = Minv'*p.Q*Minv
			J11 = Jtilde[1:nstates-nconst,1:nstates-nconst]
			Q11 = Qtilde[1:nstates-nconst,1:nstates-nconst]
			Q12 = Qtilde[1:nstates-nconst,(end-nconst+1):end]
			Q21 = Qtilde[(end-nconst+1):end,1:nstates-nconst]
			Q22 = Qtilde[(end-nconst+1):end,(end-nconst+1):end]
			Qn = Q11-Q12*inv(Q22)*Q21
			return Phs(J11, Gnull * p.B, p.D, Qn)
		elseif (det(p.G_D) > 1e8*eps(norm(p.J)))  # G_D is invertible
		## Constrained systems with non-zero, invertible, direct link
			# Given a constrained Phs with N state variables
			#     and Nconst constraints lambda:
			#                          Xdot = J GradH + B u + G lambda
			#                          y = B' GradH + D u
			#                          0 = G' GradH  + G_D lambda
			# this function finds an equivalent Phs without the constraints lambda. 
			# In the linear case, it is equivalent
			# to:
			#                  zdot = Jnew Qnew z + Bnew u
			#                  y = Bnew' Qnew z + Dnew u
			nconst = size(p.G,2)
			println("non-zero, invertible G_D:");
			Gnull = nullspace(p.G')';
			M = [Gnull; inv((p.G'*p.G)[1])*p.G'];
			Jtilde = M*p.J*M';
			Minv = inv(M);
			Qtilde = Minv'*p.Q*Minv;
			JD = zeros(size(p.J));
			JD[(end-nconst+1):end,(end-nconst+1):end] = inv(p.G_D);
			return Phs(Jtilde-JD, M*p.B, p.D, Qtilde)
		else
			println("non-zero, non-invertible G_D! don't know how to deal with it yet.");
		end
	else
	
		println("Nonlinear system. Not working yet.")
	end
end

## TO DO:
## more general case of direct links and constraints
	# Given a constrained Phs with N state variables
	#     and Nconst constraints lambda:
	#                          Xdot = J GradH + B u + G lambda
	#                          y = B' GradH + D u + D_lambda lambda
	#                          0 = G' GradH  + G_D lambda + G_u u
	# this function finds an equivalent Phs without the constraints lambda. 
	# In the linear case, it is equivalent
	# to:
	#                  zdot = Jnew Qnew z + Bnew u
	#                  y = Bnew' Qnew z + Dnew u
	#
	
## TO DO: system inversion:
## invert input/ouput of systems with non-zero, invertible D matrix:
