using PortHamiltonian

eigval, eigvec = eig(ph);
num_freq = frequencies(eigval)
exact_freq = pi/2*collect(1:2:2*N);
err = exact_freq - num_freq;

@test abs(err[2]) < 1e-13