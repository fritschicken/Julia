using Quantica, GLMakie, StaticArrays, LinearAlgebra

d_nn = 1.0
a = sqrt(3) * d_nn
M = 5
N = 4
A1 = SVector(a * cos(π/3), a * sin(π/3))
A2 = SVector(a * cos(π/3), -a * sin(π/3))
using StaticArrays

function unit_coords(N, d_nn)
    #hexagon part
    startA = SVector(0.0, -d_nn / 2)
    startB = SVector(0.0,  d_nn / 2)

    a1 = SVector(sqrt(3) / 2, -0.5)
    a2 = SVector(sqrt(3) / 2,  0.5)
    a3 = startB-startA

    coords = SVector{2, Float64}[startA, startB]

    m = d_nn .* a1 / (N - 1)
    n = d_nn .* a3 / (N - 1)
    k = d_nn .* a2 / (N - 1)


    for i in 1:(N - 2)
        pt_n = startA + i .* n
        pt_m = startA + i .* m
        pt_k = startB + i .* k
        push!(coords, pt_n, pt_m, pt_k)
    end
    return coords
end

pts = unit_coords(N, d_nn)

lat = lattice(sublat(pts), bravais = (A1, A2))

inv_A = inv(hcat(A1, A2))
function rhomboid_region(r, inv_A, M, d_nn, N)
    pt = SVector(r[1], r[2])
    n = inv_A * pt 
    del = d_nn/(N)
    return (-del <= n[1] <= M + del) && 
           (-del <= n[2] <= M + del) &&
           (abs(r[2])<= (1+(M-1)*1.5))
end


slat = lat |> supercell(region = r -> rhomboid_region(r, inv_A, M, d_nn, N))

peierls_hopping = @hopping(
    (r, dr; t = 1.0, B = 0.0) -> t * cis(0.5 * B * (r[1] * dr[2] - r[2] * dr[1])),
    range = d_nn/(N-1)+0.00002,
)

h = slat |> (@onsite((; e = 0.0) -> e) + peierls_hopping)

flux_ratios = range(0.0, 1.0, length = 301)
S_hex = (3 * sqrt(3) / 2) * (d_nn^2)
params = (; e = 0.0, t = 2.0)

b = bands(h, flux_ratios; mapping = p_ratio -> ftuple(; params..., B = 2π* p_ratio / S_hex))
fig = qplot(b, hide=:bands, color=:black)
save("Figures/hex_butterfly.png", fig)
mesh = qplot(h)
save("Figures/hex_mesh.png", mesh)