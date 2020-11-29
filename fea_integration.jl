using Proteins

include("fea.jl")
include("interface.jl")

chain = "HHHHPPPPHHHHHHHHHHHHPPPPPPHHHHHHHHHHHHPPPHHHHHHHHHHHHPPPHHHHHHHHHHHHPPPHPPHHPPHHPPHP"
polarity = HP_converter(chain)

@time E, C = folder(chain, latticetype = :triangle, ρ_1 = 0.7,stats = false, sample_limit = 50)
 
n, eidx, dofs = chain2element(C, polarity)

loads = load_nodes(C; load = [0, -1000])

A = 100 * ones(length(eidx))
E = 200e3 * ones(length(eidx))

###### Visualizing
solved = plot(chainvis(C, polarity; link = false), title = "Solved Chain")
linked = plot(chainvis(C, polarity; link = true), title = "H-H Links")
elements = plot(element_plotter(n, eidx, A; lwscale = 2), title = "Truss Conversion")

###### Solving

# First step is to ensure general structural stability, IE add links until det(K) != 0

