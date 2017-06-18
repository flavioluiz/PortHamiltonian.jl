N = 30
ph_constrained = discrete_phs(N, 0, 1); 
set_constraint!(ph_constrained, ph_constrained.B[:,1], ph_constrained.D[1:1,1:1]);
eigval, eigvec = eig(ph_constrained)
num_freq = frequencies(eigval)
exact_freq = pi*collect(-0:1:N-1);
err = exact_freq - num_freq;

@test abs(err[2]) < 1e-13

ph_c_elim = constraint_elimination(ph_constrained)
eigval, eigvec = eig(ph_c_elim);
num_freq = frequencies(eigval)
exact_freq = pi*collect(0:1:N-1);
err = exact_freq - num_freq;

@test abs(err[2]) < 1e-13