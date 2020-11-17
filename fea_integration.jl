include("branchbound.jl")
include("fea.jl")
include("interface.jl")

ch_25 = "PPHPPHHPPPPHHPPPPHHPPPPHH"
polarity = HP_converter(ch_25)

@time E, C = folder(ch_25; latticetype = :triangle, ρ_1 = 0.7,stats = false, sample_limit = 50)
 
n, eidx, dofs = chain2element(C, polarity)

loads = load_nodes(C; load = [0, -1000])

A = 100 * ones(length(eidx))
E = 200e3 * ones(length(eidx))

###### Solving
element_plotter(n, eidx, A)

elements = element_maker(n, eidx)
n_elements = length(elements)
lengths = [elem_length(e) for e in elements]

T = [Γ(θ(e)) for e in elements]
    
n_dof = length(dofs)

k_global = stiffness_init(n, elements, lengths, T, A, E)

idx_expanded = e_dof_expanded.(eidx)

K = build_K(idx_expanded, k_global, n_elements)

F = loadmaker(loads, n, n, n_dof)
