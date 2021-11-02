



function eucl_dist(s,u,λ=similar(u).=1)
    sqrt(  sum((s-u) .^2 .*λ)  )
end



function index(A::Matrix{Int}, S=AbstractAlgebra.snf(matrix(ZZ, A)))
    # Returns the index of the lattice defined by A if it is of full rank and 0 else.
    n, m = size(A)
    if n > m
        return 0
    end
    prod( [S[i,i] for i in 1:n] )
end

function homogenize(A::Matrix{Int})
    #Takes an integral matrix A and adds an additional row such that the sum over all columns is equal
    (n, m) = size(A)
    degs = [sum(A[:,i]) for i in 1:m]
    max_deg = maximum(degs)
    degs = max_deg.-degs
    if sum(degs) != 0
        return [A; transpose(degs)]
    end
    A
end

function to_system(A)
    #Takes an integral n times m matrix and returns m monomials in n variables
    (n, m) = size(A)
    @var θ[1:n]
    q=[prod(x[i]^A[i,j] for i in 1:n) for j in 1:m]
    θ, q
end




function reparametrise(A::Matrix{Int})
    # gives a parametrising map that is finite

    if ! ishomogeneous(A)
        error("A does not define a projective variety. Can not reparametrise.")
    end
    rk = rank(A)
    S,T,U = AbstractAlgebra.snf_with_transform(matrix(ZZ, A))
    T = Matrix{Int64}(Matrix(T))
    result = (T*A)[1:rk,:]
    make_positive(result)
end

function make_positive(A::Matrix{Int})
    n, m = size(A)
    mins_per_row = [minimum(A[i,:]) for i in 1:n]
    Matrix([A[i,j] - mins_per_row[i] for i in 1:n, j in 1:m])
end

function ishomogeneous(A::Matrix{Int})
    (n, m) = size(A)
    v = ones(m)
    rank(A) == rank([A; transpose(v)])
end
