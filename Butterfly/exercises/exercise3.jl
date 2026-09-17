using Quantica, GLMakie, StaticArrays
a = 1.0
N = 3
M = 15

function unit_coords(N, a)
    coords = Tuple{Float64, Float64}[(0.0, 0.0)]
    for i in 1:(N - 2)
        push!(coords, (0.0, a * i / (N - 1)))
        push!(coords, (a * i / (N - 1), 0.0))
    end
    return coords
end

pts = unit_coords(N, a)
lat = lattice(sublat(pts), bravais = ([0.0, a], [a, 0.0]))
sl = lat |> supercell(region = r -> (0 <= r[1] <= M * a) && (0 <= r[2] <= M * a))
h = sl |> (@onsite((; e = 0.0) -> e) + @hopping((r, dr; t = 1.0, B = 0.0) -> t * cis(-B * r[2] * dr[1])))


flux_ratios = range(0.0, 1.0, length = 301)
params = (; e = 0.0, t = 2.0)

b = bands(h, flux_ratios; mapping = p_ratio -> ftuple(; params..., B = 2π * p_ratio / (a^2)))

fig = qplot(b)
save("hofstadter_butterfly.png", fig)