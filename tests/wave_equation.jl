using PortHamiltonian
using PyPlot

N = 20;
a = 0.;
b = 1.;

# load the phs discretization structure 
## the Hamiltonian is considered linear with Q = I
ph = PortHamiltonian.discrete_phs(N, a, b); 
show(ph)

eigval, eigvec = eig(ph);
num_freq = PortHamiltonian.frequencies(eigval)

exact_freq = pi/2*[1:2:2*N];
err = exact_freq - num_freq;

# comparison between numerical and exact results:
[num_freq exact_freq err./exact_freq]

#plot the first five modes
plot(ph.disc_data.flow.xi, imag(eigvec[1:N,N+1:N+5]))

# constrained version

## eigenvalues of the constrained

## find an explicit equivalent system

## eigenvalues of the explicit system