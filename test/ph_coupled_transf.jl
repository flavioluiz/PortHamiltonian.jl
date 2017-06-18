p1 = PortHamiltonian.Phs([0], [1 1], [0 0; 0 0], [1])
p2 = PortHamiltonian.Phs([1], [1], [0], [1])

phct = PortHamiltonian.coupled_transformer(p1,[1], p2, [1],[1])

@test (phct.G == ones(2,1))

@test (phct.J == [0. 0.; 0. 1.])