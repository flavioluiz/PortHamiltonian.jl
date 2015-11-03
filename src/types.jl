

type Collocation_points
	xi         #  collocation points
	wi         #  quadrature integration weights
end
function show(io ::IO, object :: Collocation_points)
	println(io, "Collocation points and quadrature weights")
	println(io, "...", length(object.xi), " collocation points")
	println(io, ".xi"); println(io, ".wi")
end

type Discrete_domain
	flow     :: Collocation_points # zi
	effort   :: Collocation_points # xi
end
function show(io ::IO, object :: Discrete_domain)
	println(io, "Discretization domain points and quadrature weights")
	println(io, "...", length(object.flow.xi), " collocation points for the flow variables")
	println(io, "...", length(object.effort.xi), " collocation points for the effort variables")
	println(io, ".flow"); println(io, ".effort")
end

type Phs
	J :: Array; # interconnection matrix
	B :: Array; # control/output matrix
	D :: Array; # direct matrix
	Q  # Q matrix, optional
	G  :: Array # constraint matrix (optional)
	Hamiltonian :: Function; # Hamiltonian
	GradHam :: Function; # Hamiltonian gradient
	hessian :: Function # Hamiltonian hessian
	disc_data ::  Discrete_domain # discretization data
	
	function Phs(J :: Array, B :: Array, D :: Array, Ham :: Function)
		this = new()
		this.J = J
		this.B = B
		this.D = D
		this.Hamiltonian = Ham
		this
	end
	function Phs(J :: Array, B :: Array, D :: Array, Q)
		this = new()
		this.J = J;	this.B = B;	this.D = D;	this.Q = Q;
		Ham(x :: Array) = (0.5* x'*Q*x)[1]
		this.Hamiltonian = Ham
		GHam(x :: Array) = Q*x
		this.GradHam = GHam
		this.hessian = x->Q
		this
	end
end
function show(io ::IO, object :: Phs)
	println(io, "Port Hamiltonian system")
	println(io, "...", size(object.J,2), " energy variables")
	println(io, "...", size(object.B,2), " inputs/outputs")
	println(io, ".J"); println(io, ".B");
	if isdefined(object, :Q) println(io, ".Q"); end
	println(io, ".D"); 
	if isdefined(object, :Hamiltonian) println(io, ".Hamiltonian"); end
	if isdefined(object, :GradHam) println(io, ".GradHam"); end
	if isdefined(object, :disc_data) println(io, ".disc_data"); end
end