



struct Graph
    n_verts
    edges
    colour
    n_colours
end

function β_mod(G)
    A = zeros(Int, G.n_verts + G.n_colours^2, length(G.edges))
    for j in 1:length(G.edges)
        edge = G.edges[j]
        A[edge[1], j] += 1; A[edge[2], j] += 1
        A[n_verts + 3(G.colour[edge[1]] - 1) + G.colour[edge[2]], j] += 1
    end
    A
end


function binom_mod(G)
    A = zeros(Int, G.n_verts, length(G.edges))
    for j in 1:length(G.edges)
        edge = G.edges[j]
        A[edge[1], j] += 1; A[edge[2], j] += 1
    end
    A
end
