using PortHamiltonian
using LinearAlgebra
using Test


tests = ["ph_discretization","ph_discretization2","ph_discretization_closed",
         "ph_constraint_elim", "ph_discretization_weak",
        "ph_constrained", "ph_coupled_transf"]
for t in tests
    println(t)
    include("$(t).jl")
end
