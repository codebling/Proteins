# Initialize Packages
using LinearAlgebra, Statistics, Plots, DataFrames

#This function determines allowable candidates to place the next node
function M_k(position; dim = 2, lattice = :square, tol = 1e-4)
    last_node = position[end]
    
    # 2 Dimensional Planar
    if dim == 2
        if lattice == :square
            directions = [[1, 0], [0, 1], [-1, 0], [0, -1]]
        elseif lattice ==:triangle
            directions = [[1,0], [1/2, sqrt(3)/2], [-1/2, sqrt(3)/2], [-1, 0], [-1/2, -sqrt(3)/2], [1/2, -sqrt(3)/2]]
        end
    # 3 Dimensional
    elseif dim == 3
        if lattice == :square
            directions = [[-1, 0, 0], [1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
        elseif lattice == :triangle
            u0 = [1.0,0.0,0.0]
            u1 = [1/2, sqrt(3)/2, 0]
            u2 = [1/2, sqrt(3)/6, sqrt(6)/3]
    
            directions = [u0, -u0, u1, -u1, u2, -u2, u0-u1, u1-u0, u0-u2, u2-u0, u1-u2, u2-u1]
        end
    end

    candidates = [last_node + dir for dir in directions]
    permissibility = vcat(transpose(all([norm(existing-cad) >= tol for existing in position, cad in candidates], dims = 1))...)
    
    return candidates[permissibility]
    #Cull the candidates to locations where nodes do not already exist
end

# Determine energy of a given chain
function Energy(chain, polarity; tol = 1e-4)

    n = length(polarity)

    energy = 0

    for i = 1:n-2
        if polarity[i] == 1
            dist = [norm(chain[i] - chain[j]) for j = (i+2):n if polarity[j] == 1]
            energy += length(findall(dist .< (1 + tol)))
        end
    end

    return -energy
end

# Determine mean and minimum energies
function ZU(energies)
    U = minimum(energies)
    Z = mean(energies)
    
    return Z, U
end

# convert between an HP string of text to boolean array
function HP_converter(hp)
    out = []
    for val in hp
        if val == 'H'
            push!(out, 1)
        elseif val == 'P'
            push!(out, 0)
        else
            println("Input must be 'H' or 'P'")
            return
        end
    end
    return out
end

# Visualize statistics of search
function statsvis(stats)
    branchplot = plot(stats.k, [stats.n_best, stats.n_med, stats.n_worst, stats.n_kept], label = ["\$E_k\\leq U\$" "\$U < E_k \\leq Z\$" "\$ Z > E_k\$" "Kept Branches"], color = :black, linestyle = :auto, legend = :topleft, framestyle = :box, ylabel = "# branch candidates", xlabel = "k")

    energyplot = plot(stats.k, [stats.Z, stats.U], label = ["\$Z_k\$" "\$U_k\$"], color = :black, linestyle = :auto, legend = :bottomleft, ylabel = "Energy", framestyle = :box)

    statplot = plot(energyplot, branchplot, layout = (2,1))

    return statplot
end

function links(chain, polarity; tol = 1e-4)

    n = length(polarity)

    pos_start = []
    pos_end = []

    for i = 1:n-2
        if polarity[i] == 1
            for j = (i+2):n
                if polarity[j] == 1 && norm(chain[i] - chain[j]) < (1+tol)
                    push!(pos_start, chain[i])
                    push!(pos_end, chain[j])
                end
            end
        end
    end
    return pos_start, pos_end
end

function linkvis(baseplot, chain, polarity; linkcolor = :black, linkalpha = 0.5)

    if length(chain[1]) == 2
        dims = 2
    elseif length(chain[1]) == 3
        dims = 3
    end

    i, j = links(chain, polarity)

    n = length(i)

    if dims == 2
        for link = 1:n
            x = [i[link][1], j[link][1]]
            y = [i[link][2], j[link][2]]
            plot!(baseplot, x, y, color = linkcolor, alpha = linkalpha, linestyle = :dash)
        end
    elseif dims == 3
        for link = 1:n
            x = [i[link][1], j[link][1]]
            y = [i[link][2], j[link][2]]
            z = [i[link][3], j[link][3]]
            plot!(baseplot, x, y, z, color = linkcolor, alpha = linkalpha, linestyle = :dash)
        end
    end

    return baseplot
end

# Visualize a chain+polarity combination
function chainvis(chain, 
    polarity; 
    title = true, 
    size = 6, 
    link = true,
    linkcolor = :black,
    linkalpha = 0.5)

    if typeof(polarity) == String
        polarity = HP_converter(polarity)
    end


    h = findall(polarity .== 1)
    p = findall(polarity .!= 1)

    if title == true
        tit = "L = " * string(length(polarity))
    else
        tit = ""
    end

    # 2D vis
    if length(chain[1]) == 2
        x = [p[1] for p in chain]
        xh = x[h]
        xp = x[p]

        y = [p[2] for p in chain]
        yh = y[h]
        yp = y[p]

        form = plot(x, y, color = :black, legend = false, framestyle = :none, aspect_ratio = 1, title = tit)

        if link == true
            linkvis(form, chain, polarity; linkcolor = linkcolor, linkalpha = linkalpha)
        end

        scatter!(xh, yh, color = :black, markersize = size)
        scatter!(xp, yp, color = :white, markersize = size)
    
    # 3D vis
    elseif length(chain[1]) == 3
        x = [p[1] for p in chain]
        xh = x[h]
        xp = x[p]

        y = [p[2] for p in chain]
        yh = y[h]
        yp = y[p]

        z = [p[3] for p in chain]
        zh = z[h]
        zp = z[p]

        form = plot(x, y, z, color = :black, legend = false, framestyle = :none, aspect_ratio = 1, title = tit)

        if link == true
            linkvis(form, chain, polarity; linkcolor = linkcolor, linkalpha = linkalpha)
        end

        scatter!(xh, yh, zh, color = :black, markersize = size)
        scatter!(xp, yp, zp, color = :white, markersize = size)
    else
        println("Chain position lengths unrecognized.")
        return
    end

    return form
end


      

##############################
####### Main functions #######
##############################

function bb_fold2d(polarity, E, C, b, best, worst, mid, cands, u, z; state = [[[0,0], [1,0]]], k = 3, ρ0 = 0.1, ρ1 = 0.8, ρ2 = 0.5, lattice = :square, stats = true, threshold = 80e3)

    info1 = length(state)
    println("k = ", k, "/", length(polarity))
    println("There are ", info1, " branches")
    
    ###### Capturing all energies of all potential candidates for new chain of length k
    energies = []
    chains = []

    # If this is the last placement...
    if k == length(polarity)
        for chain in state #for all kept chains of length k-1
        
            # for a given chain of length k-1, find all possible positions for node k
            candidate_positions = M_k(chain; lattice = lattice)
    
            # Create the pseudo chains of length k at each of the candidate positions
            candidate_chains = [vcat(chain, [candidate_position]) for candidate_position in candidate_positions]
            
            # Determine the system energy for each of the pseudo chains of length k
            candidate_energies = [Energy(candidate, polarity[1:k]) for candidate in candidate_chains]
            
            # Store the energies and length-k chains
            push!(E, candidate_energies...)
            push!(C, candidate_chains...)
        end
        println("E_min: ", minimum(E))
        return 
    end

    for chain in state #for all kept chains of length k-1
        
        # for a given chain of length k-1, find all possible positions for node k
        candidate_positions = M_k(chain; lattice = lattice)

        # Create the pseudo chains of length k at each of the candidate positions
        candidate_chains = [vcat(chain, [candidate_position]) for candidate_position in candidate_positions]
        
        # Determine the system energy for each of the pseudo chains of length k
        candidate_energies = [Energy(candidate, polarity[1:k]) for candidate in candidate_chains]
        
        # Store the energies and length-k chains
        push!(energies, candidate_energies...)
        push!(chains, candidate_chains...)
    end

    #Determine the mean (Z) and minimum (U) energies of all possible chains of length k
    Z, U = ZU(energies)
    println("U = ", U)

    ####### Begin bounding process
    kept_chains = []

    # Most promising candidates
    # Extract all the minimum energy conformations
    idx_best = findall(energies .<= U)
    println("#Best: ", length(idx_best))

    if length(idx_best) > threshold
        println("Chain threshold reached.")
        idx_best_culled = rand(idx_best, Int(threshold))
        r0 = rand(length(idx_best_culled))
        push!(kept_chains, chains[idx_best_culled[r0 .> ρ0]]...)
    else
        push!(kept_chains, chains[idx_best]...)
    end

    # Worst candidates
    # Extract all the conformations with worse-than-average energies
    idx_worst = findall(energies .> Z)
    println("#Worst: ", length(idx_worst))

    if length(idx_worst) > threshold
        println("Chain threshold reached.")
        idx_worst_culled = rand(idx_worst, Int(threshold))
        r1 = rand(length(idx_worst_culled))
        push!(kept_chains, chains[idx_worst_culled[r1 .> ρ1]]...)
    else
        r1 = rand(length(idx_worst)) # Assign random probability to each of these candidates
        push!(kept_chains, chains[idx_worst[r1 .> ρ1]]...)
    end

    # Intermediate candidates
    # Extract all the conformations with intermediate energies
    idx_med = findall(U .< energies .<= Z)
    println("#Med: ", length(idx_med))

    if length(idx_med) > threshold
        println("Chain threshold reached.")
        idx_med_culled = rand(idx_med, Int(threshold))
        r2 = rand(length(idx_med_culled))
        push!(kept_chains, chains[idx_med_culled[r2 .> ρ2]]...)
    else
        r2 = rand(length(idx_med)) # Assign random probability to each of these candidates
        push!(kept_chains, chains[idx_med[r2 .> ρ2]]...)
    end

    

    # Determine how many branches are kept
    

    if stats == true
        push!(b, info1)
        push!(u, U)
        push!(z, Z)
        push!(best, length(idx_best))
        push!(worst, length(idx_worst))
        push!(mid, length(idx_med))
        push!(cands, length(kept_chains))
    end


    # Continue expansion with new set of candidates of length k
    bb_fold2d(polarity, E, C, b, best, worst, mid, cands, u, z; state = kept_chains, k = k+1, ρ1 = ρ1, ρ2 = ρ2, lattice = lattice, stats = stats, threshold = threshold)
end

function bb_fold3d(polarity, E, C, b, best, worst, mid, cands, u, z; state = [[[0,0,0], [1,0,0]]], k = 3, ρ0 = 0.1, ρ1 = 0.8, ρ2 = 0.5, lattice = :square, stats = true, threshold = 80e3)

    ###### Capturing all energies of all potential candidates for new chain of length k
    energies = []
    chains = []

    # If this is the last placement...
    if k == length(polarity)
        for chain in state #for all kept chains of length k-1
        
            # for a given chain of length k-1, find all possible positions for node k
            candidate_positions = M_k(chain; dim = 3, lattice = lattice)
    
            # Create the pseudo chains of length k at each of the candidate positions
            candidate_chains = [vcat(chain, [candidate_position]) for candidate_position in candidate_positions]
            
            # Determine the system energy for each of the pseudo chains of length k
            candidate_energies = [Energy(candidate, polarity[1:k]) for candidate in candidate_chains]
            
            # Store the energies and length-k chains
            push!(E, candidate_energies...)
            push!(C, candidate_chains...)
        end
        println("E_min: ", minimum(E))
        return 
    end

    info1 = length(state)
    println("k = ", k, "/", length(polarity))
    println("There are ", info1, " branches")

    for chain in state #for all kept chains of length k-1
        
        # for a given chain of length k-1, find all possible positions for node k
        candidate_positions = M_k(chain; dim = 3, lattice = lattice)

        # Create the pseudo chains of length k at each of the candidate positions
        candidate_chains = [vcat(chain, [candidate_position]) for candidate_position in candidate_positions]
        
        # Determine the system energy for each of the pseudo chains of length k
        candidate_energies = [Energy(candidate, polarity[1:k]) for candidate in candidate_chains]
        
        # Store the energies and length-k chains
        push!(energies, candidate_energies...)
        push!(chains, candidate_chains...)
    end

    #Determine the mean (Z) and minimum (U) energies of all possible chains of length k
    Z, U = ZU(energies)
    println("U = ", U)

    ####### Begin bounding process
    kept_chains = []

    # Most promising candidates
    # Extract all the minimum energy conformations
    idx_best = findall(energies .<= U)
    println("#Best: ", length(idx_best))

    if length(idx_best) > threshold
        println("Chain threshold reached.")
        idx_best_culled = rand(idx_best, Int(threshold))
        r0 = rand(length(idx_best_culled))
        push!(kept_chains, chains[idx_best_culled[r0 .> ρ0]]...)
    else
        push!(kept_chains, chains[idx_best]...)
    end

    # Worst candidates
    # Extract all the conformations with worse-than-average energies
    idx_worst = findall(energies .> Z)
    println("#Worst: ", length(idx_worst))

    if length(idx_worst) > threshold
        println("Chain threshold reached.")
        idx_worst_culled = rand(idx_worst, Int(threshold))
        r1 = rand(length(idx_worst_culled))
        push!(kept_chains, chains[idx_worst_culled[r1 .> ρ1]]...)
    else
        r1 = rand(length(idx_worst)) # Assign random probability to each of these candidates
        push!(kept_chains, chains[idx_worst[r1 .> ρ1]]...)
    end
    

    # Intermediate candidates
    # Extract all the conformations with intermediate energies
    idx_med = findall(U .< energies .<= Z)
    println("#Med: ", length(idx_med))

    if length(idx_med) > threshold
        println("Chain threshold reached.")
        idx_med_culled = rand(idx_med, Int(threshold))
        r2 = rand(length(idx_med_culled))
        push!(kept_chains, chains[idx_med_culled[r2 .> ρ2]]...)
    else
        r2 = rand(length(idx_med)) # Assign random probability to each of these candidates
        push!(kept_chains, chains[idx_med[r2 .> ρ2]]...)
    end

    if stats == true
        push!(b, info1)
        push!(u, U)
        push!(z, Z)
        push!(best, length(idx_best))
        push!(worst, length(idx_worst))
        push!(mid, length(idx_med))
        push!(cands, length(kept_chains))
    end

    # Continue expansion with new set of candidates of length k
    bb_fold3d(polarity, E, C, b, best, worst, mid, cands, u, z; state = kept_chains, k = k+1, ρ1 = ρ1, ρ2 = ρ2, lattice = lattice, stats = stats, threshold = threshold)
end



##############################
###### Interface #############
##############################

function folder(polarity; 
    ρ_1 = 0.8, 
    ρ_2 = 0.5, 
    dims = 2, 
    latticetype = :square, 
    sample_limit = 50e3,
    stats = true)

    # Convert "HPPHPPH" style to [1,0,0,1,0,0,1]
    if typeof(polarity) == String
        polarity = HP_converter(polarity)
    end

    # Storage vectors for statistics
    energies = [] # energies for all final conformations
    conformations = [] # layout of all final conformations
    branch_log = [] # number of candidate branches per iteration
    n_best = [] # number of candidates where E < U
    n_worst = [] # number of candidates where E > Z
    n_med = [] # number of candidates where U < E < Z
    n_kept = [] # number of candidates carried over to next iteration
    u_store = [] # value of U at each iteration
    z_store = [] # value of Z at each iteration

    # 2D HP BB algorithm
    if dims == 2
        bb_fold2d(polarity, 
            energies, 
            conformations, 
            branch_log, 
            n_best, 
            n_worst, 
            n_med, 
            n_kept, 
            u_store, 
            z_store; 
            threshold = sample_limit,
            ρ1 = ρ_1, 
            ρ2 = ρ_2, 
            lattice = latticetype,
            stats = stats)
    # 3D HP BB algorithm
    elseif dims == 3
        bb_fold3d(polarity, 
            energies, 
            conformations, 
            branch_log, 
            n_best, 
            n_worst, 
            n_med, 
            n_kept, 
            u_store, 
            z_store; 
            threshold = sample_limit,
            ρ1 = ρ_1, 
            ρ2 = ρ_2, 
            lattice = latticetype,
            stats = stats)
    else
        println("Dims should be 2 or 3")
    end

    # Processing
    e_min, idx_min = findmin(energies) #minimum energy state
    c_min = conformations[idx_min] #corresponding conformation

    # formplot = chainvis(c_min, polarity) # visualize chain

    k_store = 3:length(polarity)-1 # Storage vector

    # create data frame for iteration-wise information
    if stats == true
        output = DataFrame(k = k_store, 
            U = u_store, 
            Z = z_store, 
            n_branches = branch_log, 
            n_best = n_best, 
            n_med = n_med, 
            n_worst = n_worst, 
            n_kept = n_kept)
        
        return e_min, c_min, output
    else
        return e_min, c_min
    end

end
