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
using PyPlot


N = 10;
a = 0.;
b = 1.;

# load the phs discretization structure 
## the Hamiltonian is considered linear with Q = I
ph = discrete_phs(N, a, b); 

eigval, eigvec = eig(ph);
num_freq = frequencies(eigval)

exact_freq = pi/2*[1:2:2*N];
err = exact_freq - num_freq;

# comparison between numerical and exact results:
comp_table = [num_freq exact_freq err./exact_freq];
comp_table = comp_table[1:10,:]
@printf("\nNatural frequencies of the Dirichlet-Neumann B.C. wave equation (with u = 0)\n\n")
@printf("exact\tnumeric\terror \n")
print_table(comp_table)

#plot the first five modes
plot(ph.disc_data.flow.xi, real(eigvec[1:N,N+2]))

# constrained version
ph_constrained = deepcopy(ph);
set_constraint!(ph_constrained, ph_constrained.B, ph_constrained.D);

## eigenvalues of the constrained
eigval, eigvec = eig(ph_constrained)
num_freq = frequencies(eigval)
exact_freq = pi/2*[1:2:2*N];
err = exact_freq - num_freq;

# comparison between numerical and exact results:
comp_table = [num_freq exact_freq err./exact_freq]
@printf("\nNatural frequencies of the Neumann-Dirichlet B.C. wave equation (with y = 0)\n\n")
@printf("exact\tnumeric\terror \n")
print_table(comp_table)



## find an explicit equivalent system
ph_c_elim = constraint_elimination(ph_constrained)
eigval, eigvec = eig(ph_c_elim);
num_freq = frequencies(eigval)

exact_freq = pi/2*[1:2:2*N];
err = exact_freq - num_freq;
# comparison between numerical and exact results:
comp_table = [num_freq exact_freq err./exact_freq]
comp_table = comp_table[1:10,:]
@printf("\nNatural frequencies of the Neumann-Dirichlet B.C. wave equation (with y=0 after constraint elimination)\n\n")
@printf("exact\tnumeric\terror \n")
print_table(comp_table)
