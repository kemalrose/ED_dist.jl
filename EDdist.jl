

using HomotopyContinuation

export EDdist, EDdeg, to_system


function homogenize(A::Matrix{Int})
    (n, m) = size(A)
    #Takes an integral matrix A and adds an additional row such that the sum over all columns is equal
    degs = [sum(A[:,i]) for i in 1:m]
    max_deg = maximum(degs)
    degs = max_deg.-degs
    if sum(degs) != 0
        return [A; transpose(degs)]
    end
    A
end

#function EDdist(x, q, u, λ=ones(length(q)) )
#    D = differentiate(q,x)
#    F = System(transpose(D)*(q-u))

#    S = solve(F, only_non_zero=true)
#    q_sys = System(q)
#    R = solutions(S)
#    h = findall(r->all(abs.(r).>1e-8), R)
#    qs = map(s->q_sys(s),R[h])
#    #certify(F,R[h])
#    unique_points(qs)
#end
hm

function EDdist(x::Vector{Variable}, q, u::Vector{Float64}, λ=ones(length(q)) )
    (certification, critical_pts, real_critical_pts, index_min, distance) = criticals(x,q,u,λ)
    distance, real_critical_pts[index_min]
end

function EDdist(A::Matrix, u::Vector{Float64}, λ=ones(length(q)) )
    x, q = to_system(A)
    EDdist(x,q,u,λ)
end

function EDdeg(x, q, generic=false)
    if generic
        λ, u = randn(Float64, length(q)), randn(Float64, length(q))
        (certification, critical_pts_on_X, real_critical_pts, index_min, distance) = criticals(x,q,u,λ)
        return length(critical_pts_on_X)
    end
    (certification, critical_pts_on_X, real_critical_pts, index_min, distance) = criticals(x,q)
    print("There are ",
    length(critical_pts_on_X), " distinct critical points, and ", length(real_critical_pts), " real ones.")
    length(critical_pts_on_X)
end


function to_system(A)
    #Takes an integral n times m matrix and m monomials in n variables
    (n, m) = size(A)
    @var x[1:n]
    q=[prod(x[i]^A[i,j] for i in 1:n) for j in 1:m]
    x, q
end

function EDdeg(A::Matrix, generic=false)
    x, q = to_system(A)
    EDdeg(x,q,generic)
end



function criticals(x::Vector{Variable}, q, u::Vector{Float64}=randn(Float64,length(q)), λ::Vector=ones(length(q)))
    g = sum(λ[i] * (q[i] - u[i])^2 for i in 1:length(q))
    F = System(differentiate(g,x))
    S = solve(F)
    R = solutions(S)
    q_sys = System(q)
    g_sys = System([g])

    nonzero_indices = findall(r->all(abs.(r).>1e-8), R)
    critical_pts = R[nonzero_indices]
    critical_pts_on_X = unique_points(map(s->q_sys(s),critical_pts))
    real_indices = findall(s -> all(abs.(imag.(s)).<1e-8), critical_pts_on_X)
    distances = map(s->real(g_sys(s)[1]), critical_pts[real_indices])

    if length(distances) == 0
        error("No real critical points.")
    end
    distance, index_min = findmin(distances)

    certification = certify(F, critical_pts)
#    if length(real_indices) > nreal_certified(certification)
#        error("Found more real critical_pts than there can be.")
#    end
    (certification, critical_pts_on_X, critical_pts_on_X[real_indices], index_min, distance)
end
