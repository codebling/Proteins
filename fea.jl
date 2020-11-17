using LinearAlgebra, Plots, NLopt, SparseArrays

# Geometric Functions
#Get length of element
function elem_length(element)
    return norm(element[2] - element[1])
end

function element_maker(nodes, element_indices)
    return [[nodes[e_idx[1]], nodes[e_idx[2]]] for e_idx in element_indices]
end


# get angle of element
function θ(element)
    x_global = [1, 0] #unit x vector
    vec_element = element[2] - element[1]

    cos = dot(x_global, vec_element) / (norm(x_global) * norm(vec_element))

    angle = acosd(cos)

    if vec_element[2] < 0
        angle = 360 - angle
    end

    return angle
end

## Transformation matrix
function Γ(ang)
    return [cosd(ang) sind(ang) 0 0;0 0 cosd(ang) sind(ang)]
end

function e_dof_expanded(e_idx)
    x_idx = e_idx .* 2 .- 1
    y_idx = x_idx .+ 1
    
    edof = zeros(Int, 4)
    edof[1:2:end-1] = x_idx
    edof[2:2:end] = y_idx
    
    return edof
end


## Loads
function loadidx(pos, nodes)
    # will create a load value to the closest node at the chosen position
    # pos = [x, y]
    
    load_idx = findmin([sum(abs.(n .- pos)) for n in nodes])[2]
    
    return load_idx
end

function loadmaker(loads, positions, nodes, n_dofs)
    if length(loads) !== length(positions)
        println("load and position vectors must be of equal length.")
        return
    end
    #find the nodal index for each load
    idx = [loadidx(pos, nodes) for pos in positions]
    idx_expanded = vcat([[i*2-1, i*2] for i in idx]...)
    loads_expanded = vcat(loads...)
    
    f = zeros(n_dofs)
    f[idx_expanded] = loads_expanded
    return f
end

##Stiffness

## Local coordinate elemental stiffness
function k_localxy(length, A, E)
    k = A * E / length .* [1 -1; -1 1]
    return k
end

# Convert local stiffness matrices to global stiffness matrices
function k_globalxy(k_local::Array{Float64, 2}, T)
    #Global_xy stiffness matrix
    k = transpose(T) * k_local * T
    
    #Note the indices of the matrix are in the form of:
    # 1: node A, global X
    # 2: node A, global Y
    # 3: node B, global X
    # 4: node B, global Y

    return k
end

# Assemble global stiffness matrix
function build_K(element_node_idx, k_elemental, n)
    # Indices in global stiffness matrix
    idxx = vcat([vcat([ei * ones(Int, 4) for ei in eidx]...) for eidx in element_node_idx]...) 
    idxy = vcat([repeat(ei, 4) for ei in element_node_idx]...)

    kvals = vcat([vcat(reshape(k_elemental[i], (16,1))...) for i = 1:n]...)

    return sparse(idxx, idxy, kvals)
end

## Solving

function displacements(K, F, dof_list, n_dof)
    #Reduce K matrix and force vector to non-zero items:
    K_reduced = K[dof_list, dof_list]
    F_reduced = F[dof_list]
    
    disp = Symmetric(K_reduced) \ F_reduced


    compliance = transpose(F_reduced) * disp #F^T d 
    
    #Extended disp vector
    disp_long = zeros(n_dof)
    disp_long[dof_list] = disp
    
    return disp_long, compliance
end

##Post Processing
function reactions(F_e, n_dofs, element_idx)
    #create empty storage matrix
    F_all = zeros(length(dofs))
    
        #for each element
        for i = 1:length(element_idx)
            #extract the global XY elemental end forecs
            F = F_e[i]
    
            #extract the given element information
            e = element_idx[i]
    
            #Determine the start and end nodes for the given element
            n_start = e[1]
            n_end = e[2]
    
            #Convert to the indexes in the global DOF list 
            idx_start = 2 * n_start - 1 .+ [0,1]
            idx_end = 2 * n_end - 1 .+ [0,1]
    
            #Add the end forces of each element to the appropriate DOF:
            F_all[idx_start] .+= F[1:2]
            F_all[idx_end] .+= F[3:4]
        end
        
    #Return the force values of each DOF that are at REACTIONS
    return F_all .* (dofs .== 0)
end

function forces(n_elements, elements_idx, A, k_global, disp, T, n_dofs)
    #Global element forces
    F_e = [k_global[i] * disp[i] for i = 1:n_elements]
    
    #Convert global element forces t olocal element forces
    f_e = [T[i] * F_e[i] for i = 1:n_elements]
    
    #Single axial force value
    f = [force[2] for force in f_e]
    
    #Stress values
    σ = [f[i] / A[i] for i = 1:n_elements]
    
    #Determing Reaction Forces
    rxns = reactions(F_e, n_dofs, elements_idx)
    
    return F_e, f_e, f, σ, rxns
end


function eq_check(rxns, F, tolerance)
    #Check equilibrium
    if sum(rxns) + sum(F) > tolerance
        println("Equilibrium not met!!")
        return true
    else
        println("Equilibrium met.")
    end
end


function displaced_nodes(nodes, disp; scale_x = 1, scale_y = 1)
    n_disp = []
    for i = 1:length(nodes)
        n = nodes[i]
        n_new = [(n[1] + disp[i*2 - 1]) * scale_x, (n[2] + disp[i*2]) * scale_y]
        push!(n_disp, n_new)
    end
    return n_disp
end

## Secondary
function fem_init(nodes, elements)
    #extract lengths of each element
    lengths = [elem_length(e) for e in elements]
    #extract CC angle (degrees) of each element from global axis
    T = [Γ(θ(e)) for e in elements]
    
    return lengths, T
end

function stiffness_init(nodes, elements, lengths, T, A, E)

    #Stiffness Matrices
    k_element_local = [k_localxy(lengths[i], A[i], E[i]) for i = 1:length(elements)]
    k_element_global = [k_globalxy(k_element_local[i], T[i]) for i = 1:length(elements)]
    
    return k_element_global
end

## Main

## Main solver
function fem2d_solver(nodes, dofs, element_idx, A, E, loads, positions; tol = 1e-5)
    
    elements = element_maker(nodes, element_idx)
    n_elements = length(elements)
    lengths = [elem_length(e) for e in elements]
    
    T = [Γ(θ(e)) for e in elements]
    
    n_dof = length(dofs)

    #Return global stiffness matrix (including 0 rows/columns)
    k_global = stiffness_init(nodes, elements, lengths, T, A, E)
    
    idx_expanded = e_dof_expanded.(eidx)
    
    K = build_K(idx_expanded, k_global, n_elements)
    #Create force vector
    F = loadmaker(loads, positions, nodes, n_dof)
    
    #Displacements in global coordinate system + elemental DOF displacements
    disp, compliance = displacements(K, F, dofs, n_dof)
    
    disp_local = [disp[idx] for idx in idx_expanded]
    
    F_global, f_elemental, f_axial, stress_axial, rxns = forces(n_elements, element_idx, A, k_global, disp_local, T, n_dof)
    
    #If equilibrium is not reached, exit solver
    if eq_check(rxns, F, tol) == true
        return
    end
    
    #Deformed nodes and elements
    new_nodes = displaced_nodes(nodes, disp)
    
    return new_nodes, disp, f_axial, stress_axial, rxns, compliance
end

## Main solver
function analysis(nodes, dofs, element_idx, A, E, loads, positions; tol = 1e-5)
    
    elements = element_maker(nodes, element_idx)
    n_elements = length(elements)
    lengths = [elem_length(e) for e in elements]
    
    T = [Γ(θ(e)) for e in elements]
    
    n_dof = length(dofs)

    #Return global stiffness matrix (including 0 rows/columns)
    k_global = stiffness_init(nodes, elements, lengths, T, A, E)
    
    idx_expanded = e_dof_expanded.(element_idx)
    
    K = build_K(idx_expanded, k_global, n_elements)

    #Create force vector
    F = loadmaker(loads, positions, nodes, n_dof)
    
    #Displacements in global coordinate system + elemental DOF displacements
    disp, compliance = displacements(K, F, dofs, n_dof)
    
    disp_local = [disp[idx] for idx in idx_expanded]
    
    F_global, f_elemental, f_axial, stress_axial, rxns = forces(n_elements, element_idx, A, k_global, disp_local, T, n_dof)
    
    #If equilibrium is not reached, exit solver
    if eq_check(rxns, F, tol) == true
        return
    end
    
    max_stress_index = findmax(abs.(stress_axial))[2]

    σ_max = maximum(stress_axial)
    disp_max = maximum(disp)
    
    return σ_max, disp_max, compliance
end

## For topology Optimization
function topop_compliance(elements, k_glob_sens, idx_expanded, dofs, A, E, F; tol = 1e-5)
    
    n_elements = length(elements)
    
    n_dof = length(dofs)

    #Return global stiffness matrix (including 0 rows/columns)
    k_global = k_glob_sens .* A
    
    K = build_K(idx_expanded, k_global, n_elements)

    #Displacements in global coordinate system + elemental DOF displacements
    disp, compliance = displacements(K, F, dofs, n_dof)
    
    disp_local = [disp[idx] for idx in idx_expanded]

    sens = [-transpose(disp_local[i]) * k_glob_sens[i] * disp_local[i] for i = 1:n_elements]
    
    return compliance, sens
end

function topop_compliance2(nodes, elements, idx_expanded, dofs, A, E, F; tol = 1e-5)
    
    n_elements = length(elements)
    
    n_dof = length(dofs)


    T = [Γ(θ(e)) for e in elements]
    k_global = stiffness_init(nodes, elements, element_lengths, T, A, E)

    K = build_K(idx_expanded, k_global, n_elements)
    
    
    #Displacements in global coordinate system + elemental DOF displacements
    disp, compliance = displacements(K, F, dofs, n_dof)
    
    disp_local = [disp[idx] for idx in idx_expanded]
    
    return compliance
end

##################
# Meshers #
##################
function allmesher(L, H, nx::Int, ny::Int, A, E; fix = :SIMPLE, nodelabel = true)
    #L = length (width) of mesh domain
    #H = height of mesh domain
    #nx = number of nodes in x (L) direction
    #ny = number of nodes in y (H) direction
    
    x_spacing = L/(nx-1) #spacing between nodes in x direction
    y_spacing = H/(ny-1) #spacing between nodes in y direction
    
    #node set
    nodes = []
    dofs = []
    for x = 1:nx
        for y = 1:ny
            #create a node with default released DOF
            push!(nodes, [(x-1)*x_spacing, (y-1)*y_spacing])
            push!(dofs, [true, true])
        end
    end
    
    n_nodes = length(nodes)
    node_idx = collect(1:n_nodes)
    n_dof = 2 * n_nodes

    # Set element fixities
    if fix == :LEFT
        dofs[1:ny] .= [[false, false]]
    elseif fix == :LEFTSYMM
        dofs[1:ny] .= [[false, true]]
        dofs[end - ny + 1] = [true, false]
    elseif fix == :LEFTSYMMFIX
        dofs[1:ny] .= [[false, true]]
        dofs[end - ny + 1] = [false, false]
    elseif fix == :RIGHT
        dofs[end-ny+1:end] .= [[false, false]]
    elseif fix == :SIMPLE
        dofs[1] = [false, false]
        dofs[end - ny + 1] = [true, false]
    elseif fix == :TOP
        dofs[ny:ny:end] .= [[false, false]]
    elseif fix == :BOTTOM
        dofs[1:ny:end] .= [[false, false]]
    elseif fix == :NONE
    else
        println("Unknown fixity.")
        return
    end

    
    #element set
    eidx = []
    for i = 1:n_nodes-1
        e_set = [[node_idx[i], node_idx[j]] for j = i+1:1:n_nodes]
        push!(eidx, e_set)
    end
    
    eidx = vcat(eidx...)
    n_elements = length(eidx)
    
    return nodes, vcat(dofs...), eidx, A * ones(n_elements), E * ones(n_elements)
            
    
end

function gridmesher(L, H, nx::Int, ny::Int, A, E; fix = :SIMPLE, nodelabel = true)
    #L = length (width) of mesh domain
    #H = height of mesh domain
    #nx = number of nodes in x (L) direction
    #ny = number of nodes in y (H) direction
    
    x_spacing = L/(nx-1) #spacing between nodes in x direction
    y_spacing = H/(ny-1) #spacing between nodes in y direction
    
    #node set
    nodes = []
    dofs = []
    for x = 1:nx
        for y = 1:ny
            #create a node with default released DOF
            push!(nodes, [(x-1)*x_spacing, (y-1)*y_spacing])
            push!(dofs, [true, true])
        end
    end
    
    n_nodes = length(nodes)
    node_idx = collect(1:n_nodes)
    n_dof = 2 * n_nodes
    
    # Set element fixities
    if fix == :LEFT
        dofs[1:ny] .= [[false, false]]
    elseif fix == :LEFTSYMM
        dofs[1:ny] .= [[false, true]]
        dofs[end - ny + 1] = [true, false]
    elseif fix == :LEFTSYMMFIX
        dofs[1:ny] .= [[false, true]]
        dofs[end - ny + 1] = [false, false]
    elseif fix == :RIGHT
        dofs[end-ny+1:end] .= [[false, false]]
    elseif fix == :SIMPLE
        dofs[1] = [false, false]
        dofs[end - ny + 1] = [true, false]
    elseif fix == :TOP
        dofs[ny:ny:end] .= [[false, false]]
    elseif fix == :BOTTOM
        dofs[1:ny:end] .= [[false, false]]
    elseif fix == :NONE
    else
        println("Unknown fixity.")
        return
    end


    #element creation
    
    #This reshapes the array of nodes as a matrix where the position of the node value
    #Has the same relative position as the actual physical system ie:
    # node2 node4 node6
    # node1 node3 node5
    # for a 2x3 node system.
    node_idx_matrix = reverse(reshape(node_idx, (ny,nx)), dims = 1)
    
    #empty element set
    elements = []
    
    #Horizontal elements
    for row = 1:ny
        for col = 1:nx-1
            e_start = node_idx_matrix[row,col]
            e_end = node_idx_matrix[row,col+1]
            push!(elements, [e_start, e_end])
        end
    end
    
    #vertical elements
    for col = 1:nx
        for row = 1:ny-1
            e_start = node_idx_matrix[row,col]
            e_end = node_idx_matrix[row+1,col]
            push!(elements, [e_start, e_end])
        end
    end
    
    #\ diagonal elements
    for row = 1:ny-1
        for col = 1:nx-1
            e_start = node_idx_matrix[row,col]
            e_end = node_idx_matrix[row+1,col+1]
            push!(elements, [e_start, e_end])
        end
    end
    
    #/ diagonal elements
    for row = 2:ny
        for col = 1:nx-1
            e_start = node_idx_matrix[row,col]
            e_end = node_idx_matrix[row-1,col+1]
            push!(elements, [e_start, e_end])
        end
    end
    
    n_elements = length(elements)
    
            
    return nodes, vcat(dofs...), elements, A * ones(n_elements), E * ones(n_elements)
            
    
end

function element_plotter(nodes, element_idx, A; cull_limit = 1, show_nodes = true, lwscale = 5)
    elems = element_maker(nodes, element_idx)
    x = [[elems[i][1][1], elems[i][2][1]] for i = 1:length(elems)]
    y = [[elems[i][1][2], elems[i][2][2]] for i = 1:length(elems)]

    # cull procedure
    idx = findall(A .> cull_limit)

    nx = [[node[1]] for node in nodes]
    ny = [[node[2]] for node in nodes]
    
    widths = transpose(A[idx] ./ maximum(A[idx])) * lwscale

    fig = plot(x[idx], y[idx], 
    lw = widths,
    color = :black, 
    framestyle = :none, 
    legend = false, 
    alpha = 0.7,
    aspect_ratio = 1)

    if show_nodes
        scatter!(nx, ny, color = :black)
    end

    return fig
end


############################
######Topology Optimization
############################
function topop(A_initial, nodes, elements_idx, degrees_of_freedom, E, loads, positions, V_max)

    elements = element_maker(nodes, elements_idx)
    element_lengths = elem_length.(elements)
    idx_expanded = e_dof_expanded.(elements_idx)

    if length(A_initial) != length(elements)
        println("A_initial must equal length of elements.")
        return
    end

    n_dims =  length(elements)

    n_dof = length(degrees_of_freedom)
    F = loadmaker(loads, positions, nodes, n_dof)

    T = [Γ(θ(e)) for e in elements]
    k_global_sensitivity = stiffness_init(nodes, elements, element_lengths, T, ones(length(elements)), E)

    local function obj_func(A::Vector, grad::Vector)
    
        comp, sens = topop_compliance(elements, k_global_sensitivity, idx_expanded, degrees_of_freedom, A, E, F)

        if length(grad) > 0
            grad[:] = sens
        end

        return comp
    end
    
    local function cstr_vol(A::Vector, grad::Vector)
        if length(grad) > 0
            grad[:] = element_lengths
        end

        sum(A .* element_lengths) - V_max
    end
    
    opt = Opt(:LD_MMA, n_dims)
    opt.lower_bounds = 1e-3 * ones(n_dims)
    opt.xtol_abs = 1e-4

    opt.min_objective = obj_func
    inequality_constraint!(opt, (x,g) -> cstr_vol(x,g), 1e-8)

    opt_compliance, opt_A, ret = optimize(opt, A_initial)
    
    return opt_compliance, opt_A
end

#####################
# More Plotting
#####################

function loadplotter(loadstore, scale)
    if typeof(loadstore) !== Array{load,1}
        println("No loads specified/incorrect type.")
        return [([0,0], [0,0])]
    end
    #Normalize length of load arrows
    normalizer = maximum([sqrt(load.x_val^2 + load.y_val^2) for load in loadstore])
    
    L = length(loadstore)
    
    x1(i) = loadstore[i].n.x
    x2(i) = loadstore[i].n.x + loadstore[i].x_val / normalizer * scale

    y1(i) = loadstore[i].n.y
    y2(i) = loadstore[i].n.y + loadstore[i].y_val / normalizer * scale
    
    
    loads = [([x1(i), x2(i)], [y1(i), y2(i)]) for i = 1:L]
    
    return loads
end

function bc_plotter(nodes)
    #X-direction BCs
    x_fixed_idx = findall([node.dof_x == 0 for node in nodes])
    #Y-direction BCs
    y_fixed_idx = findall([node.dof_y == 0 for node in nodes])
    
    x_bcs = [([nodes[i].x], [nodes[i].y]) for i in x_fixed_idx]
    y_bcs = [([nodes[i].x], [nodes[i].y]) for i in y_fixed_idx]
    
    return x_bcs, y_bcs
end

function trussplotter(nodes, elements, loads; 
        nodelabels = false,
        node_marker = :circle,
        node_label_color = :red,
        node_label_pos = :right,
        nodesize = 3, 
        nodecolor = :black, 
        elementlabels = false,
        elementcolor = :gray,
        elementstyle = :solid,
        lwscale = 2,
        loadscale = :auto,
        bcscale = 10,
        bccolor = :green)
    
    #initiate plot
    canvas = plot()
    
    #plot Boundary Conditions
    xbc, ybc = bc_plotter(nodes)
    xmarker = [(-1, sqrt(2)/2), (0,0), (-1,-sqrt(2)/2)]
    ymarker = [(-sqrt(2)/2,-1), (0,0), (sqrt(2)/2,-1)]
    
    scatter!(xbc, marker = (Shape(xmarker), bcscale, bccolor))
    scatter!(ybc, marker = (Shape(ymarker), bcscale, bccolor))
    
    #plot elements
    #Linewidth scale
    areas = [element.A for element in elements]
    widths = transpose(areas ./ maximum(areas)) * lwscale

    plot!(elementcollector(elements), 
        color = elementcolor, 
        linestyle = elementstyle, 
        lw = widths,
        label = "")

    #plot element labels
    if nodelabels
        for i = 1:length(nodes)
            annotate!((nodes[i].x, nodes[i].y, Plots.text(string(i), node_label_color, node_label_pos)))
        end
    end

    #plot element labels
    if elementlabels
        for i = 1:length(elements)
            annotate!((element_mid(elements[i])[1], element_mid(elements[i])[2], Plots.text(string(i), :black)))
        end
    end

    #plot loads
    #Determine scale for loads
    if loadscale == :auto
        xrange = maximum(node.x for node in nodes) - minimum(node.x for node in nodes)
        yrange = maximum(node.y for node in nodes) - minimum(node.y for node in nodes)

        scale = 0.10 * max(xrange, yrange)
    else
        scale = loadscale
    end
    
    plot!(loadplotter(loads, scale), 
        color = :red, 
        arrows = true,
        label = "")
    
    #plot nodes
    scatter!(nodecollector(nodes), 
        markershape = node_marker,
        markercolor = nodecolor, 
        label = "")
    
    canvas
    return canvas
end