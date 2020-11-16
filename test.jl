include("branchbound.jl")


## Chen and Huang Benchmark Proteins
ch_24 = "HHPPHPPHPPHPPHPPHPPHPPHH"
ch_25 = "PPHPPHHPPPPHHPPPPHHPPPPHH"
ch_36 = "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP"
ch_48 = "PPHPPHHPPHHPPPPPHHHHHHHHHHPPPPPPHHPPHHPPHPPHHHHH"
ch_60 = "PPHHHPHHHHHHHHPPPHHHHHHHHHHPHPPPHHHHHHHHHHHHPPPPHHHHHHPHHPHP"
ch_85 = "HHHHPPPPHHHHHHHHHHHHPPPPPPHHHHHHHHHHHHPPPHHHHHHHHHHHHPPPHHHHHHHHHHHHPPPHPPHHPPHHPPHPH"
ch_100 = "PPPHHPPHHHHPPHHHPHHPHHPHHHHPPPPPPPPHHHHHHPPHHHHHHPPPPPPPPPHPHHPHHHHHHHHHHHPPHHHPHHPHPPHPHHHPPPPPPHHH"

ch_huge = "HHPPHPHPHHPHPHPHHHPPHHHPPPHHHPPHHPHPHHHPPPHPHHHPPPPPHHHHHHHPPHHPHPHHPHPHHHHPHPHHHPHPHHPHPHHHPPHPHHHHPPPPHHPHPHPPHPHHPPHHPPPPHHHPPPHPPPHPPPPPPHHHPPHHPPHHHHPPPHHPHPPPHPHPHHPHHPPPPPPPHHHPPHPPPHPPHPHHPPPH"

ch_huge2 = "HHPPHPHPHHPHPHPHHHPPHHHPPPHHHPPHHPHPHHHPPPHPHHHPPPPPHHHHHHHPPHHPHPHHPHPHHHHPHPHHHPHPHHPHPHHHPPHPHHHHPPPPHHPHPHPPHPHHPPHHPPPPHHHPPPHPPPHPPPPPPHPPHHPPPPHHHPPPPHHHHHHHHHHHHHHHHHHHHHHHHHHHHHPPPHPPPHPPPHHHHPPPHHHHHHHPPPPPPPPHHPPPHHHPPHHPPHHHHPPPHHPHPPPHPHPHHPHHPPPPPPPHHHPPHPPPHPPHPHHPPPH"

## Solving:

# 3D Triangular
@time E_3dt, C_3dt = folder(ch_huge2; dims = 3, latticetype = :triangle, ρ_1 = 0.7,stats = false, sample_limit = 50)

# 3D Square
@time E_3ds, C_3ds = folder(ch_huge2; dims = 3, ρ_1 = 0.7,stats = false, sample_limit = 50)

# 2D Triangular
@time E_2dt, C_2dt = folder(ch_huge2; latticetype = :triangle, ρ_1 = 0.7,stats = false, sample_limit = 50)

# 2D Square
@time E_2ds, C_2ds = folder(ch_huge2; ρ_1 = 0.7,stats = false, sample_limit = 50)

### full visualization
growth_all = @animate for i = 1:length(C)
    P1 = plot(chainvis(C_2ds[1:i], ch_huge2[1:i]; size = 5, linkalpha = 0.6), foreground = :white)
    P2 = plot(chainvis(C_2dt[1:i], ch_huge2[1:i]; size = 5, linkalpha = 0.6), foreground = :white)
    P3 = plot(chainvis(C_3ds[1:i], ch_huge2[1:i]; size = 5, linkalpha = 0.3), foreground = :white)
    P4 = plot(chainvis(C_3dt[1:i], ch_huge2[1:i]; size = 5, linkalpha = 0.3), foreground = :white)

    plot(P1, P2, P3, P4, layout = (2,2))
end

pyplot()

chain = C_3dt
polarity = HP_converter(ch_huge2)

h = findall(polarity .== 1)
p = findall(polarity .!= 1)

x = [p[1] for p in chain]
xh = x[h]
xp = x[p]

y = [p[2] for p in chain]
yh = y[h]
yp = y[p]

z = [p[3] for p in chain]
zh = z[h]
zp = z[p]

form = plot(x, y, z, color = :black, legend = false, foreground = :white)
scatter!(xh, yh, zh, color = :black)

formtest = @animate for i = 1:length(chain)
    plot(x[1:i], y[1:i], z[1:i], camera = (i,30), color = :black, legend = false, foreground = :white)
end

gif(formtest, "Figures/formtest.gif", fps = 20)

