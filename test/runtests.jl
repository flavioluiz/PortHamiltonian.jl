using PortHamiltonian
using Base.Test

# write your own tests here
N = 30
ph = discrete_phs(N, 0, 1); 
eigval, eigvec = eig(ph);
num_freq = frequencies(eigval)
exact_freq = pi/2*collect(1:2:2*N);
err = exact_freq - num_freq;

@test abs(err[2]) < 1e-13

