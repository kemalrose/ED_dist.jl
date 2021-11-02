

include("ED_module.jl")
using .ED_dist



n_verts = 7
edges = [[1,2], [1,6], [1,7], [2,3], [3,4], [4,6], [5,6], [5,7],[6,7]]
colour = [1,1,2,2,3,3,3]
n_colours = 3
G = Graph(n_verts, edges, colour, n_colours)
A = binom_mod(G)

#EDdeg(reparametrise(A))
#@test EDdeg(A) == 76
B = β_mod(G)
#@test EDdeg(B) == 1


A = reparametrise(B)
EDdeg(A)
u = randn(Float64, 9)
EDdist(A, u)


S,T,U = Matrix{Int64}.(Matrix.(AbstractAlgebra.snf_with_transform(  matrix(ZZ, B)  )))


rk = rank(B)
result = inv(T)[:,1:rk] * S[1:rk,:] * inv(U)
(T*B)[1:rk,:]

B = [1 1 1 0 0 0;
     1 0 0 1 1 0;
     0 1 0 1 0 1;
     0 0 1 0 1 1]

B = [0 0;
     1 0]


B = homogenize([1 1 1 0 0 0 0;
                1 0 0 1 1 0 0;
                0 1 0 1 0 1 0;
                0 0 1 0 1 1 0])
A=[0 1 0 0 1 1 0 1 1 1 1 1;
   0 1 1 1 0 0 1 0 1 1 1 1;
   1 0 0 1 1 0 1 1 0 1 1 1;
   1 0 1 0 0 1 1 1 1 0 1 1;
   0 0 1 1 1 1 1 1 1 1 0 1]
 A=homogenize([0 1 0 0 1 1 0 1 1 1 1 1 0;
    0 1 1 1 0 0 1 0 1 1 1 1 0;
    1 0 0 1 1 0 1 1 0 1 1 1 0;
    1 0 1 0 0 1 1 1 1 0 1 1 0;
    0 0 1 1 1 1 1 1 1 1 0 1 0])
x,q = to_system(A)
λ = randn(Float64, length(q))
u = append!(randn(Float64, length(q)-1), 1)

function crits(x::Vector{Variable}, q::Vector, u::Vector{Float64}=randn(Float64,length(q)), λ::Vector=ones(length(q)), critical_pts=Nothing)
    if critical_pts == Nothing
        critical_pts = criticals(x,q,u,λ)
    end
    g = sum(λ[i] * (q[i] - u[i])^2 for i in 1:length(q))
    q_sys = System(q)
    g_sys = System([g])
    critical_pts_on_X = unique_points(map(s->q_sys(s),critical_pts))
    real_indices = findall(s -> all(abs.(imag.(s)).<1e-8), critical_pts_on_X)
    distances = map(s->real.(g_sys(s)[1]), critical_pts[real_indices])
    distance, index_min = findmin(distances)
    print(critical_pts[real_indices])

    CriticalPointsOnX(critical_pts, critical_pts_on_X, critical_pts_on_X[real_indices], index_min, distance)
end






n_verts = 10
edges = [[1,2], [1,6], [1,7], [1,8], [2,3], [3,4],[3,10], [4,6], [5,6], [5,7],[6,7], [8,9], [8,10], [9,10]]
colour = [1,1,2,2,3,3,3,4,4,4]
n_colours = 4
G = Graph(n_verts, edges, colour, n_colours)
A = binom_mod(G)

#EDdeg(reparametrise(A))
#@test EDdeg(A) == 76
B = β_mod(G)
#@test EDdeg(B) == 1


A = reparametrise(B)
EDdeg(A)
u = randn(Float64, 9)
EDdist(A, u)





function IncidenceCompleteGraph(n:: Int)
    Incidence = ones(Int, n, n)
    for i in 1:n
        Incidence[i, i] = 0
    end
    Incidence
end





IncidenceKarate = [
 0 1 1 1 1 1 1 1 1 0 1 1 1 1 0 0 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 0;
 1 0 1 1 0 0 0 1 0 0 0 0 0 1 0 0 0 1 0 1 0 1 0 0 0 0 0 0 0 0 1 0 0 0;
 1 1 0 1 0 0 0 1 1 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 0;
 1 1 1 0 0 0 0 1 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1;
 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 1 1;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 1 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1;
 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1;
 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1 1;
 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 0 0 0 1 1;
 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 1 0 0 1 0 1 0 1 1 0 0 0 0 0 1 1 1 0 1;
 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1 1 0 0 1 1 1 0 1 1 0 0 1 1 1 1 1 1 1 0
]

IncidenceKarate = IncidenceKarate[1:15,1:15]



A = [0 1 1 1; 1 0 0 0; 1 0 0 0; 1 0 0 0]
Bin  = binom_mod(IncidenceKarate)
colours = rand(1:2, size(IncidenceKarate)[1])
Beta = β_mod(IncidenceKarate, colours)

Edges = edge_list(IncidenceKarate)
print(Beta[end, :])
print([colours[edge[1]] != colours[edge[2]] for edge in Edges])




#H,U = HomotopyContinuation.hnf(transpose(Bin))

#A=transpose(H)
A = reparametrise(Bin)

m,n = size(A)

@var θ[1:m]


q=[prod(θ[i]^A[i,j] for i in 1:m) for j in 1:n]


Q = System(q)

k = length(q)

@var u[1:k]

G = differentiate(q,θ)

λ = randn(ComplexF64, k)

D = diagm(λ)

F = transpose(G)

M = F*D
#we need θ₀ and u₀ to use monodromy_solve
θ₀ = randn(ComplexF64,m)
# define u₀
M₀ = [subs(M[i,j],θ => θ₀) for (i,j) in Iterators.product(1:m ,1:k)]
w = nullspace(convert(Matrix{ComplexF64},M₀))
n = w[:,1]
q₀ = Q(θ₀)
u₀ = q₀ - n

# the system where u are the parameters and θ are the variables
S = System(M*(q-u), variables=θ, parameters=u)
d(θ₁, θ₂) = norm(Q(θ₁) - Q(θ₂), Inf)
@time MS = monodromy_solve(S, [θ₀], u₀,distance = d)

R = solutions(MS)

### test
qs = map(s->Q(s),R)
length(unique_points(qs)) == 36



certify(S, R, target_parameters = u₀)
