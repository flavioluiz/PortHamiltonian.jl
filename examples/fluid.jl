using PortHamiltonian
using PyPlot

include("dataexperiment.jl")

# load structural parameters
pslosh, pbeam, prigid = dataexperiment(0.4);

L = pslosh["a"];
N = 20;
disc = discrete_phs_closed(N, 0, L)

Br = [eye(N); zeros(N,N);zeros(1,N)];
B = blkdiag(zeros(2*N,1),[1]);
D = blkdiag(disc.D,[0]);
Jrb = [0];
invM = disc.Q[N+1:end,N+1:end]
R = 000000.00*Br*norm(invM)*Br';
J = blkdiag(disc.J,Jrb)-R;

# defining Q for saint-venant equation (no rotation)    
rho = pslosh["rho"];
g = pslosh["g"];
b = pslosh["b"];
hbar = pslosh["h"];
xi,w = PortHamiltonian.lgwt(N,0,L);
Q = blkdiag(disc.Q,[0])
pfluid = PortHamiltonian.Phs(J, B, D, Q);
#PortHamiltonian.set_constraint!(pfluid, B[:,1],[0])

function Hamiltonian(X)
	alpha1 = X[1:N]
	alpha2 = X[N+1:2*N]
	p = X[end]
	M = w' * (alpha1 .* alpha2)
	H = (w'*((alpha1 .* (alpha2.^2))/rho + rho*g*((alpha1).^2)/b)/2 + ((p-M).^2)/(2*prigid["mass"]))[1]
	return H
end

grad = ForwardDiff.gradient(x->Hamiltonian(x))
hess = ForwardDiff.hessian(x->Hamiltonian(x))
xeq = [ones(N,1)*b*hbar;zeros(N,1);0.001*prigid["mass"]*0][:]
Qlin = hess(xeq[:])
pfluid.Q = Qlin
pfluid.Hamiltonian = x-> Hamiltonian(x)
pfluid.GradHam = x-> grad(x)
pfluid.hessian = x-> hess(x)

#pfluid.GradHam = x -> Q*x
function dynamics(t, X, Xd)
	Xd[1:end]=  pfluid.J * pfluid.GradHam(X[1:end]) + 0.2*pfluid.B[:,2] * sin(2*pi*0.5*t) 
	#r[end:end] = pfluid.G' * pfluid.hessian(X[1:end-1]) * (pfluid.J * pfluid.GradHam(X[1:end-1]) + pfluid.B[:,3] * sin(2*pi*t) *0.5+ pfluid.G * X[end])
end

#
# using Sundials
#
# t = (linspace(0,8,1000))
# yout = Sundials.cvode(dynamics,xeq, [t])
# #yout, ypout = Sundials.idasol(dynamics2, [xeq;l0;0], [xeq;0;0]*0, [t])
# pfluid.GradHam = x-> Qlin*x
# youtlin = Sundials.cvode(dynamics,xeq, [t])
# #youtlin, ypoutlin = Sundials.idasol(dynamics2, [xeq;l0; 0], [xeq;0; 0]*0, [t])
# #pn = PortHamiltonian.constraint_elimination(pfluid)
# using PyPlot
#
# figure(1);
# surf(xi,t,yout[:,N+1:end-1])
# #surf(xi,t,yout[:,1:N])
# #include("anima.jl")