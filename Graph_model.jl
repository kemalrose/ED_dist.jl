


function β_mod(Incidence::Matrix{Int}, colours = ones(Int, size(Incidence)[1]) )
    [binom_mod(Incidence); new_rows(Incidence, colours)]
end

function edge_list(Incidence::Matrix{Int})
    n,m = size(Incidence)
    edges = reshape( [(j, i) for i in 1:n, j in 1:n], n^2 )
    filter(x -> x[1]<x[2] && Incidence[x[1], x[2]] != 0, edges)
end

function binom_mod(Incidence::Matrix{Int})
    n, m = size(Incidence)
    if n != m
        error("Input matrix is not square")
    end

    edges = edge_list(Incidence)
    A = zeros(Int, n, length(edges))
    for k in 1:length(edges)
        (i,j) = edges[k]
        A[i, k] += Incidence[i, j]
        A[j, k] += Incidence[i, j]
    end
    A
end

function new_rows(Incidence::Matrix{Int}, colours)
    edges = edge_list(Incidence)
    n_verts = size(Incidence)[1]
    n_colours = maximum(colours)

    colour_pairs  = reshape( [(j, i) for i in 1:n_colours, j in 1:n_colours], n_colours^2 )
    colour_pairs = filter(x -> x[1]<x[2], colour_pairs)

    New_Rows = zeros(Int, length(colour_pairs), length(edges))

    for j in 1:length(edges)
        edge = edges[j]
        col1, col2 = colours[edge[1]], colours[edge[2]]
        col1, col2 = ( col1 < col2 ? (col1, col2) : (col2, col1) )
        if col1 != col2
            New_Rows[ findfirst(x -> x == (col1, col2), colour_pairs), j] += 1
        end
    end
    New_Rows
end

struct Graph
    n_verts
    edges
    colour
    n_colours
end

function new_rows(G::Graph)
    N = zeros(Int, G.n_colours^2, length(G.edges))
    for j in 1:length(G.edges)
        edge = G.edges[j]
        if G.colour[edge[1]] != G.colour[edge[2]]
            N[G.n_colours*(G.colour[edge[1]] - 1) + G.colour[edge[2]], j] += 1
        end
    end
    a_less_b = filter( s -> s[1] < s[2], [(a,b) for a in 1:G.n_colours, b in 1:G.n_colours])
    allowed_column_indices = [G.n_colours*(a - 1) + b for (a,b) in a_less_b]
    N[allowed_column_indices,:]
end

function β_mod(G::Graph)
    [binom_mod(G); new_rows(G)]
end

function binom_mod(G::Graph)
    A = zeros(Int, G.n_verts, length(G.edges))
    for j in 1:length(G.edges)
        edge = G.edges[j]
        A[edge[1], j] += 1; A[edge[2], j] += 1
    end
    A
end


function clique(n::Int)
    edges = [[i,j] for i in 1:n, j in 1:n]
    edges = filter(edge->(edge[1] < edge[2]), edges)
    Graph(n, edges, Int64.(ones(n)), 1)
end
