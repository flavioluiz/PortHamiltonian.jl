

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
	Bd :: Array; # distributed control matrix (optional)
	D :: Array; # direct matrix
	Q  # Q matrix, optional
	R :: Array; # damping matrix (optional)
	G  :: Array # constraint matrix (optional)
	G_D :: Array # direct term constraint matrix (optional)
	Hamiltonian :: Function; # Hamiltonian
	GradHam :: Function; # Hamiltonian gradient
	hessian :: Function # Hamiltonian hessian
	disc_data ::  Discrete_domain # discretization data
	TransfMatrix :: Array # transformation matrix after reduction change of variables
							# always start as identity
	StatesNames :: Dict
	InputsNames :: Array
	auxvars :: Dict # auxiliary variables dictionary
	
	function Phs(J :: Array, B :: Array, D :: Array, Ham :: Function)
		this = new()
		this.J = J
		this.B = B
		this.D = D
		this.R = zeros(size(J))
		this.Hamiltonian = Ham
		this.TransfMatrix = eye(size(J,1))
		this.StatesNames = Dict()
		this.auxvars = Dict()
		this.InputsNames = repmat([""],size(B,2))
		this
	end
	function Phs(J :: Array, B :: Array, D :: Array, Q)
		this = new()
		this.J = J;	this.B = B;	this.D = D;	this.Q = Q;
		this.R = zeros(size(J))
		Ham(x :: Array) = (0.5* x'*Q*x)[1]
		this.Hamiltonian = Ham
		GHam(x :: Array) = Q*x
		this.GradHam = GHam
		this.hessian = x->Q
		this.TransfMatrix = eye(size(J,1))
		this.StatesNames = Dict()
		this.auxvars = Dict()
		this.InputsNames = repmat([""],size(B,2))
		this
	end
end
function set_constraint!(ph :: Phs, G, G_D)
	if isdefined(ph, :G) warn("Constraint was already defined! overwriting!"); end
	ph.G = G
	ph.G_D = G_D
end
function set_constraint!(ph ::Phs, G)
	set_constraint!(ph, G, zeros(size(G,2),size(G,2)))
end
function show(io ::IO, object :: Phs)
	println(io, "Port Hamiltonian system")
	println(io, "...", size(object.J,2), " energy variables")
	println(io, "...", size(object.B,2), " inputs/outputs")
	println(io, ".J"); println(io, ".B");
	println(io, ".D"); 
	if isdefined(object, :R) println(io, ".R"); end
	if isdefined(object, :Q) println(io, ".Q"); end
	if isdefined(object, :G) println(io, ".G (constrained system)"); end
	if isdefined(object, :Hamiltonian) println(io, ".Hamiltonian"); end
	if isdefined(object, :GradHam) println(io, ".GradHam"); end
	if isdefined(object, :disc_data) println(io, ".disc_data"); end
	if isdefined(object, :Bd) println(io, ".Bd (distributed ports)"); end	
end