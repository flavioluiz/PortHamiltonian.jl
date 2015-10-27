# Legendre-Gauss quadrature weights and nodes
# This code was originally written for Matlab
#  by Greg von Winckel - 02/25/2004
# posted on Mathworks website and licensed under BSD License
# http://fr.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes
# 
# 

function lgwt(N,a,b)

	N=N-1;
	N1=N+1; N2=N+2;
	xu=linspace(-1,1,N1)';

	y=(cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2))';

	L=zeros(N1,N2);

	Lp=zeros(N1,N2);
	y0=2;

	while maximum(abs(y-y0))>eps()				
		L[:,1]=1;
		Lp[:,1]=0;		
		L[:,2]=y;
		for k=2:N1
			L[:,k+1]=( (2*k-1)*y.*L[:,k]-(k-1)*L[:,k-1] )/k;
		end	 
		Lp=(N2)*( L[:,N1]-y.*L[:,N2] )./(1-y.^2);   		
		y0=y;
		y=y0-L[:,N2]./Lp;		
	end

	x=(a*(1-y)+b*(1+y))/2;      

	w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

	return x,w
end