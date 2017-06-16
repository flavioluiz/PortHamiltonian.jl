#
# this file shows an example of pseudo-spectral discretization for PHS systems for wave equation
#  three cases are presented. It should be used for verifying the discretization, eigenvalues,
#  and constraint elimination methods. Three cases are tested:
#         1) PHS with u = 0 (Neumann-Dirichlet B.C)
#         2) PHS with y = 0 (Neumann-Dirichlet B.C, but with a constraint)
#         3) PHS with y = 0 after removing the constraints.
#

# auxiliary function: used for printing results

function print_table(table)
	num_y = size(table,1)
	num_x = size(table,2)
	for i = 1:num_y
		for j = 1:num_x
			if j < 3
				@printf("%.2f\t", table[i,j])
			else
				@printf("%.2e", table[i,j])
			end
		end
		@printf("\n")
	end
end

# modules:
using PortHamiltonian
#using PyPlot


N = 10;
a = 0.;
b = 1.;

# load the phs discretization structure 
## the Hamiltonian is considered linear with Q = I
ph = discrete_phs(N, a, b); 

eigval, eigvec = eig(ph);
num_freq = frequencies(eigval)

exact_freq = pi/2*collect(1:2:2*N);
err = exact_freq - num_freq;

# comparison between numerical and exact results:
comp_table = [num_freq exact_freq err./exact_freq];
comp_table = comp_table[1:10,:]
@printf("\nNatural frequencies of the Dirichlet-Neumann B.C. wave equation (with u = 0)\n")
@printf("exact\tnumeric\terror \n")
print_table(comp_table)
@printf("\t From a total %d frequencies, %d are represented with error less than 0.01\n", length(err), sum(abs(err) .< 0.01)) 
@printf("\n\t %d are smaller than %.1e\n",sum(abs(err) .< eps()*1e6),eps()*1e6) 
#plot the first five modes
#plot(ph.disc_data.flow.xi, real(eigvec[1:N,N+2]))

# constrained version
ph_constrained = deepcopy(ph);
set_constraint!(ph_constrained, ph_constrained.B[:,1], ph_constrained.D[1:1,1:1]);

## eigenvalues of the constrained
eigval, eigvec = eig(ph_constrained)
num_freq = frequencies(eigval)
exact_freq = pi*collect(-0:1:N-1);
err = exact_freq - num_freq;

# comparison between numerical and exact results:
comp_table = [num_freq exact_freq err./exact_freq]
comp_table = comp_table[1:10,:]
@printf("\nNatural frequencies of the Neumann-Dirichlet B.C. wave equation (with y = 0)\n")
@printf("exact\tnumeric\terror \n")
print_table(comp_table)
@printf("\t\t From a total %d frequencies, %d are represented with error less than 0.01\n", length(err), sum(abs(err) .< 0.01)) 
@printf("\n\t %d are smaller than %.1e\n",sum(abs(err) .< eps()*1e6),eps()*1e6) 
## find an explicit equivalent system
ph_c_elim = constraint_elimination(ph_constrained)
eigval, eigvec = eig(ph_c_elim);
num_freq = frequencies(eigval)

exact_freq = pi*collect(0:1:N-1);
err = exact_freq - num_freq;
# comparison between numerical and exact results:
comp_table = [num_freq exact_freq err./exact_freq]
comp_table = comp_table[1:10,:]
@printf("\nNatural frequencies of the Neumann-Dirichlet B.C. wave equation (with y=0 after constraint elimination)\n")
@printf("exact\tnumeric\terror \n")
print_table(comp_table)
@printf("\t\t From a total %d frequencies, %d are represented with error less than 0.01", length(err), sum(abs(err) .< 0.01)) 
@printf("\n\t %d are smaller than %.1e\n",sum(abs(err) .< eps()*1e6),eps()*1e6) 
