# only linear PHS for now
# coupling using gyrator-type interconnection
function coupled_gyrator(ph1 :: Phs, ports1, ph2 :: Phs, ports2, couple_matrix)
	N1 = size(ph1.J,1); # number of states
	N2 = size(ph2.J,1);
	
	# split the input ports that are not used for interconnection
	# and creates the new B matrix:
	inputs1 = [];
	for input = 1:size(ph1.B,2)
		if ~(input in ports1)
			inputs1 = push!(inputs1,input);
		end
	end
	inputs2 = [];
	for input = 1:size(ph2.B,2)
		if ~(input in ports2)
			inputs2 = push!(inputs2,input);
		end
	end
	Bnew = blkdiag(ph1.B[:,inputs1],ph2.B[:,inputs2]);
	Bnew[1:N1,(length(inputs1)+1):end] = ph1.B[:,ports1]*couple_matrix*ph2.D[ports2,inputs2];
	Bnew[(N1+1):end, 1:length(inputs1)] = -ph2.B[:,ports2]*couple_matrix*ph1.D[ports1,inputs1] ;
	
	Jnew = blkdiag(ph1.J, ph2.J);
	Jnew[1:N1, (N1+1):end] = ph1.B[:,ports1]*couple_matrix*(ph2.B[:,ports2].')
	Jnew[(N1+1):end, 1:(N1)] = -ph2.B[:,ports2]*couple_matrix*(ph1.B[:,ports1].')
	
	Dnew = zeros(length(inputs1)+length(inputs2),length(inputs1)+length(inputs2));
	Dnew[1:length(inputs1), 1:length(inputs1)] = ph1.D[inputs1,inputs1]
	Dnew[length(inputs1)+1:end, length(inputs1)+1:end] = ph2.D[inputs2,inputs2]
	Dnew[1:length(inputs1), length(inputs1)+1:end] = ph1.D[inputs1,ports1]*couple_matrix*ph2.D[ports2, inputs2];
	Dnew[length(inputs1)+1:end, 1:length(inputs1)] = - ph2.D[inputs2, ports2]*couple_matrix*ph1.D[ports1,inputs1];
	Qnew = blkdiag(full(ph1.Q), full(ph2.Q));
	phsnew = Phs(Jnew, Bnew, Dnew, Qnew)
	phsnew.TransfMatrix = blkdiag(ph1.TransfMatrix,ph2.TransfMatrix)
	
	# preserves the names of the states:
	N1 = size(ph1.TransfMatrix,1)
	for key in keys(ph2.StatesNames)
		ph2.StatesNames[key] += N1
	end
	phsnew.StatesNames = merge(ph1.StatesNames, ph2.StatesNames)
	
	# remove the names of the coupled inputs
	inp1 = copy(ph1.InputsNames)
	deleteat!(inp1, ports1);
	inp2 = copy(ph2.InputsNames)
	deleteat!(inp2, ports2);	
	phsnew.InputsNames = [inp1, inp2]
	
	phsnew.R = blkdiag(ph1.R, ph2.R)	
	return phsnew
end

function coupled_transformer(ph1 :: Phs, ports1, ph2 :: Phs, ports2, couple_matrix)
	N1 = size(ph1.J,1); # number of states
	N2 = size(ph2.J,1);
	
	# split the input ports that are not used for interconnection
	# and creates the new B matrix:
	inputs1 = [];
	for input = 1:size(ph1.B,2)
		if ~(input in ports1)
			inputs1 = push!(inputs1,input);
		end
	end
	inputs2 = [];
	for input = 1:size(ph2.B,2)
		if ~(input in ports2)
			inputs2 = push!(inputs2,input);
		end
	end
	Bnew = blkdiag(ph1.B[:,inputs1],ph2.B[:,inputs2]);
	#Bnew[1:N1,(length(inputs1)+1):end] = ph1.B[:,ports1]*couple_matrix*ph2.D[ports2,inputs2];
	#Bnew[(N1+1):end, 1:length(inputs1)] = -ph2.B[:,ports2]*couple_matrix*ph1.D[ports1,inputs1] ;
	
	Jnew = blkdiag(ph1.J, ph2.J);
	#Jnew[1:N1, (N1+1):end] = ph1.B[:,ports1]*couple_matrix*(ph2.B[:,ports2].')
	#Jnew[(N1+1):end, 1:(N1)] = -ph2.B[:,ports2]*couple_matrix*(ph1.B[:,ports1].')
	
	#Dnew = zeros(length(inputs1)+length(inputs2),length(inputs1)+length(inputs2));
	#Dnew[1:length(inputs1), length(inputs1)+1:end] = ph1.D[inputs1,ports1]*couple_matrix*ph2.D[ports2, inputs2];
	#Dnew[length(inputs1)+1:end, 1:length(inputs1)] = - ph2.D[inputs2, ports2]*couple_matrix*ph1.D[ports1,inputs1];
	G = [ph1.B[:,ports1],ph2.B[:,ports2]*couple_matrix]
	GD = ph1.D[ports1,ports1] - ph2.D[ports2,ports2]*couple_matrix
	Dnew = blkdiag(ph1.D[inputs1,inputs1], ph2.D[inputs2,inputs2])
	Qnew = blkdiag(full(ph1.Q), full(ph2.Q));
	phnew = Phs(Jnew, Bnew, Dnew, Qnew)
	set_constraint!(phnew, G, GD)
	phnew.TransfMatrix = blkdiag(ph1.TransfMatrix,ph2.TransfMatrix)
	
	# preserves the names of the states:
	N1 = size(ph1.TransfMatrix,1)
	for key in keys(ph2.StatesNames)
		ph2.StatesNames[key] += N1
	end
	phnew.StatesNames = merge(ph1.StatesNames, ph2.StatesNames)
	phnew.R = blkdiag(ph1.R, ph2.R)
	return phnew
end
