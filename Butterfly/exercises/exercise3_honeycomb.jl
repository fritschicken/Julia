using Quantica, GLMakie, StaticArrays, LinearAlgebra

d_nn = 1.0
a = sqrt(3) * d_nn
M = 4

A1 = SVector(a * cos(π/3), a * sin(π/3))
A2 = SVector(-a * cos(π/3), a * sin(π/3))

sA = sublat((0.0,  d_nn / 2), name = :A)
sB = sublat((0.0, -d_nn / 2), name = :B)
lat = lattice(sA, sB, bravais = (A1, A2))


hex_center_0 = SVector((A1+A2)/2)

# centers of the M x M grid of hexagons
hex_centers = [hex_center_0 + n1 * A1 + n2 * A2 for n1 in 0:(M - 1), n2 in 0:(M - 1)]

# find the related sites of each center
function is_on_hex_mesh(r, centers, d_nn)
    pt = SVector(r[1], r[2])
    for c in centers
        if isapprox(norm(pt - c), d_nn, atol = 1e-4)
            return true
        end
    end
    return false
end

slat = lat |> supercell(region = r -> is_on_hex_mesh(r, hex_centers, d_nn))

peierls_hopping = @hopping(
    (r, dr; t = 1.0, B = 0.0) -> t * cis(0.5 * B * (r[1] * dr[2] - r[2] * dr[1])),
    range = 1.05,
    sublats = :A => :B
)

h = slat |> (@onsite((; e = 0.0) -> e) + plusadjoint(peierls_hopping))

qplot(h(B = 0.2))


flux_ratios = range(-1.0, 1.0, length = 301)
S_hex = (3 * sqrt(3) / 2) * (d_nn^2)
params = (; e = 0.0, t = 2.0)
b = bands(h, flux_ratios; mapping = p_ratio -> ftuple(; params..., B = 2π * p_ratio / S_hex))
qplot(slat)
fig = qplot(b, hide=:nodes)
save("hex_butterfly.png", fig)
qplot(h(B=0.2))