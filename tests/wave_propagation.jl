using PortHamiltonian

# mixed finite elements (Golo)
#Nel = 10; alfa = 0.5; a = 0; b = 1;
#p = PortHamiltonian.new_spectral_element_phs_golo(Nel,alfa,a,b)

# pseudo-spectral discretization (Moulla)
p = discrete_phs(10,0,1)

# pseudo-spectral "elements"
#p = PortHamiltonian.new_spectral_element_phs(1,40,0,1)

# pseudo-spectral elements mixed with golo alpha=1 elements:
#p = PortHamiltonian.new_block_spectral_element_phs(Nel,Npol,a,b)

Ns = size(p.J,1)
freq = 1.*pi/2

function dynamics(t,x, xd)
	xd[:] = p.J * p.GradHam(x) + p.B[:,2] * sin(freq*t) - p.B[:,1] * output(t,x)
end
function output(t,x)
	return (p.B.' * p.GradHam(x) + p.D[:,2] * sin(freq*t))[1]
end

using Sundials
to = collect(linspace(0,10/freq,1000))
x0 = zeros(Ns);
res = Sundials.cvode(dynamics, x0, to)

y = zeros(length(to))
for i = 1:length(y)
	y[i] = output(to[i], res[i,:][:])
end

using PyPlot

figure()
plot(to+1, sin(freq*to),to, y)
figure()
plot(to, y-(sin(freq*(to-1)).*(to.>=1)  ))