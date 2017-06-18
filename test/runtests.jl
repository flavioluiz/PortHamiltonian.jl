using PortHamiltonian
using Base.Test


tests = ["ph_discretization",
        "ph_constrained"]
for t in tests
    println(t)
    include("$(t).jl")
end
