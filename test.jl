using Proteins

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

polarity = HP_converter(ch_huge2)
# 3D Triangular
@time E_3dt, C_3dt = folder(ch_huge2; dims = 3, latticetype = :triangle, ρ_1 = 0.7,stats = false, sample_limit = 50)
chainvis(C_3dt, polarity)

# 3D Square
@time E_3ds, C_3ds = folder(ch_huge2; dims = 3, ρ_1 = 0.7,stats = false, sample_limit = 50)

# 2D Triangular
@time E_2dt, C_2dt = folder(ch_huge2; latticetype = :triangle, ρ_1 = 0.7,stats = false, sample_limit = 50)
chainvis(C_2dt, polarity)

# 2D Square
@time E_2ds, C_2ds = folder(ch_huge2; ρ_1 = 0.7,stats = false, sample_limit = 500)
chainvis(C_2ds, polarity)

## full animation
# growth_all = @animate for i = 1:length(ch_huge2)
#     P1 = plot(chainvis(C_2ds[1:i], ch_huge2[1:i]; size = 3, linkalpha = 0.6), foreground = :white)
#     P2 = plot(chainvis(C_2dt[1:i], ch_huge2[1:i]; size = 3, linkalpha = 0.6), foreground = :white)
#     P3 = plot(chainvis(C_3ds[1:i], ch_huge2[1:i]; size = 3, linkalpha = 0.3), camera = (i, 30), foreground = :white)
#     P4 = plot(chainvis(C_3dt[1:i], ch_huge2[1:i]; size = 3, linkalpha = 0.3), camera = (i,30), foreground = :white)

#     plot(P1, P2, P3, P4, layout = (2,2))
# end

# gif(growth_all, "Figures/growth_pyplot.gif", fps = 20)

