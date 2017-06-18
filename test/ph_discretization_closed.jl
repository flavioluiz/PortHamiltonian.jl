
N = 10
ph = discrete_phs_closed(N, 0, 1); 
num_freq = frequencies(ph)

@test abs(num_freq[2]-pi) < 1e-13