function dermatrix(N :: Integer,a :: Number,b :: Number, order = 1)
	dermatrix(N+1, N, a, b, order)
end

function dermatrix(Nx :: Integer, Nz :: Integer, a :: Number, b :: Number, order :: Integer)
	xi,wi = lgwt(Nx,a,b)
	if Nz > 1
		zi,wiz = lgwt(Nz,a,b)
	else
		zi,wiz = (a+b)/2, (b-a);
	end
	D = dermatrix(xi,zi,order)
	return D,xi,zi
end

function dermatrix(xi :: Array, zi :: Number, order = 1)
	dermatrix(xi, [zi], order)
end

function dermatrix(xi :: Array,zi :: Array, order = 1)
	D = zeros(length(zi),length(xi))
	for j = 1:length(xi)
		g = algo_diff(x-> leg_pol(x,xi[:],j), order)
		for i = 1:length(zi)		
			D[i,j] = g(zi[i])
		end
	end
	D
end

function algo_diff(func :: Function, order :: Integer)
	if order == 1
		g = x-> ForwardDiff.derivative(func,x)
		return g
	elseif order == 2
		g = x->ForwardDiff.hessian(func,x)
		g(x) = g([x])
		return g
		else
			println("no more than 2nd order derivative yet")
	end
end

function dermatrix2(N)
	xi = cos(pi*(1:1:(N+2))/(N+3))
	zi = cos(pi*(1:1:(N))/(N+1))
	D = zeros(length(zi),length(xi))
	for j = 1:length(xi)
		g = algo_diff(x-> leg_pol(x,xi,j),2)
		for i = 1:length(zi)		
			D[i,j] = g([zi[i]])
		end
	end
	return D,xi,zi
end
