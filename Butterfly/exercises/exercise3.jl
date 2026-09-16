using Quantica, GLMakie, StaticArrays
a = 2
N = 4
M=3

function unit_coords(N,a)
    coords = Tuple{Number, Number}[(0,0)]
    for i in 1:N-2
        push!(coords, (0, a*i/(N-1)))
        push!(coords, (a*i/(N-1), 0))
    end
    return coords
end

pts = unit_coords(N,a)
lat = lattice(sublat(pts), bravais=([0,a], [a,0]))
sl = lat |> supercell(region = r -> (0 <= r[1] <= M*a) && (0 <= r[2] <= M*a))
h = sl |> (@onsite((; ϵ = 0.0) -> ϵ) + @hopping((r, dr; t = 1.0, B = 0.0) -> t * cis(-B * r[2] * dr[1])))

Bs = range(-pi/2,pi/2, length=199)
params = (;ϵ=0, t=2, B=0)
b = bands(h, Bs; mapping = B ->ftuple(;params...,B))
qplot(h)
#qplot(b,hide=:nodes)