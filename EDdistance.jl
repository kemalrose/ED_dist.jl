
struct CriticalPointsOnX
    critical_pts_on_T
    critical_pts_on_X
    real_critical_pts
    index_min
    distance
end


Base.show(io::IO, X::CriticalPointsOnX) = print(io, "There are ", length(X.critical_pts_on_X),
" distinct critical points on X and ", length(X.real_critical_pts), " of them are real.")

function EDdeg(x::Vector{Variable}, q::Vector, generic::Bool=false)
    if generic
        λ, u = randn(Float64, length(q)), randn(Float64, length(q))
        result = criticals_on_X(x,q,u,λ)
    else
        result = criticals_on_X(x,q)
    end
    result
end


function EDdeg(A::Matrix, generic::Bool=false)

    # Compute the degree of the parametrising morphism.
    deg_φ = index(A)
    n, m = size(A)
    if generic
        λ, u = randn(Float64, m), randn(Float64, m)
    else
        λ, u = ones(m), randn(Float64, m)
    end

    n_crit = length(criticals(A,u,λ))
    if ! (n_crit % deg_φ == 0)
        error("Number of critical points is not divisible by the degree of φ. Critical points: ", n_crit, " degree: ", deg_φ)
    end
    n_crit ÷ deg_φ
end

function criticals(A::Matrix, u::Vector{Float64}=randn(Float64,size(A)[2]), λ::Vector=ones(size(A)[2]))
    θ, q = to_system(A)
    criticals(θ,q,u,λ)
end

function criticals(x::Vector{Variable}, q::Vector, u::Vector{Float64}=randn(Float64,length(q)), λ::Vector=ones(length(q)))
    g = sum(λ[i] * (q[i] - u[i])^2 for i in 1:length(q))
    F = System(differentiate(g,x))
    R = solutions(HomotopyContinuation.solve(F, show_progress=false))
    On_T = filter(r->all(abs.(r).>1e-8), R)
    unique_points(On_T)
end

function criticals_on_X(x::Vector{Variable}, q::Vector, u::Vector{Float64}=randn(Float64,length(q)), λ::Vector=ones(length(q)), critical_pts=Nothing)
    if critical_pts == Nothing
        critical_pts = criticals(x,q,u,λ)
    end
    g = sum(λ[i] * (q[i] - u[i])^2 for i in 1:length(q))
    q_sys = System(q)
    g_sys = System([g])
    critical_pts_on_X = unique_points(map(s->q_sys(s),critical_pts))
    real_indices = findall(s -> all(abs.(imag.(s)).<1e-8), critical_pts_on_X)
    reals = real.(critical_pts_on_X[real_indices])

    distances = map(s-> eucl_dist(s,u,λ), reals)
    distance, index_min = findmin(distances)

    CriticalPointsOnX(critical_pts, critical_pts_on_X, reals, index_min, distance)
end

function criticals_on_X(A::Matrix{Int}, u::Vector{Float64}=randn(Float64,size(A)[2]), λ::Vector=ones(size(A)[2]), critical_pts=criticals(A,u,λ))
    x, q = to_system(A)
    criticals_on_X(x,q,u,λ,critical_pts)
end

function EDdist(x::Vector{Variable}, q::Vector, u::Vector{Float64}, λ=ones(length(q)))
    result = criticals_on_X(x,q,u,λ)
    result.distance, result.real_critical_pts[result.index_min]
end

function EDdist(A::Matrix, u::Vector, λ=ones(size(A)[2]))
    x, q = to_system(A)
    EDdist(x,q,u,λ)
end
