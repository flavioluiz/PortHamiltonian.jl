# write your own tests here
N = 10
ph = PortHamiltonian.weak_phs1(N, 0, 1); 
eigval, eigvec = eig(ph);
num_freq = frequencies(eigval)

@test abs(num_freq[2]-pi) < 1e-13

N = 10
ph = PortHamiltonian.weak_phs2(N, 0, 1); 
eigval, eigvec = eig(ph);
num_freq = frequencies(eigval)

@test abs(num_freq[2]-pi) < 1e-13


ph = PortHamiltonian.finelem(N+1,5, 0, 1)
eigval, eigvec = eig(ph);
num_freq = frequencies(eigval)
@test abs(num_freq[2]-pi) < 1e-13

ph = PortHamiltonian.finelem(N,5, 0, 1)
eigval, eigvec = eig(ph);
num_freq = frequencies(eigval)
@test abs(num_freq[2]-pi*3/2) < 1e-13