using PortHamiltonian
using PyPlot
using PyCall
@pyimport matplotlib.patches as patch

cfig = figure()
ax = cfig[:add_subplot](1,1,1)

ax[:set_aspect]("equal")
Nint = 5
# weak form
Nelem = 20
p = finelem(Nelem,Nint,0,1)

for i = 1:Nelem
	if (i%2 == 1)
	c = patch.Rectangle(
		(-0.5+(i-1)*(Nint*2+1), -0.5+(i-1)*(Nint*2+1)), Nint*2+1, Nint*2+1,
		alpha=0.3,fc="blue")
	else
	c = patch.Rectangle(
		(-0.5+(i-1)*(Nint*2+1), -0.5+(i-1)*(Nint*2+1)), Nint*2+1, Nint*2+1,
		alpha=0.3,fc="red")
	end
	ax[:add_artist](c)
end
spy(p.J)

pM = PortHamiltonian.new_spectral_element_phs(Nelem,Nint+1,0,1)

figure();
spy(pM.J)
freqp = frequencies(p)
freqpM = frequencies(pM)
[freqp[1:10] freqp[1:10]-[1:2:20]*pi/2 freqpM[1:10] freqpM[1:10]-[1:2:20]*pi/2]