# write your own tests here
N = 10
ph = discrete_phs2(N, 0, 1); 
num_freq = frequencies(ph)

@test abs(num_freq[1]- 3.51602) < 1e-3

N = 10
K = 7
a = 0
b = 1
ad = 0
bd = 0.1
phport = discrete_phs2_distports(N,K,a,b,ad,bd)
Bd = phport.Bd*ones(K,1)
@test abs(Bd[2] - 110.799) < 1e-2