using Quantica, GLMakie, StaticArrays
a = 1
N = 4
M=1
L = (N-1) * M

h = LP.square() |> supercell(region = r -> ((r[1] % (N-1) == 0 &&  0<=r[1] <= L && 0 <= r[2] <= L) || 
                           (r[2] % (N-1) == 0 &&  0<=r[2] <= L && 0 <= r[1] <= L)))
h = h |> (@onsite((; ϵ = 0.0) -> ϵ) + @hopping((r, dr; t = 1.0, B = 0.0) -> t * cis(-B * r[2] * dr[1]), range = 1.0))


Bs = range(-pi/2,pi/2, length=199)
params = (;ϵ=0, t=2, B=0)
b = bands(h, Bs; mapping = B ->ftuple(;params...,B))


qplot(h(B=0.1))
qplot(b,hide=:nodes)