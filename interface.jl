function E_idx(chain)
    # this function creates proper node indexing for FEA analysis
    n = length(chain)
    idx = []
    for i = 1:n-1
        push!(idx, [i, i+1])
    end

    return idx
end

function link_idx(chain, polarity; tol = 1e-4)
    n = length(polarity)

    idx = []

    for i = 1:n-2
        if polarity[i] == 1
            for j = (i+2):n
                if polarity[j] == 1 && norm(chain[i] - chain[j]) < (1+tol)
                    push!(idx, [i, j])
                end
            end
        end
    end
    return idx
end

function chain2element(chain, polarity; fix = :INIT)
    n = chain
    eidx = vcat(E_idx(chain), link_idx(chain, polarity))
    dofs = Array(trues(2*length(n)))

    if fix == :INIT
        dofs[1:4] .= false
    end
    return n, eidx, dofs
end

function load_nodes(chain; load = [0, -1])
    return repeat([load], length(chain))
end