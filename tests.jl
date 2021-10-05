
using Test


@testset "EDdegree_examples" begin
    # Jux Cantor example
    A=[0 1 0 0 1 1 0 1 1 1 1 1;
       0 1 1 1 0 0 1 0 1 1 1 1;
       1 0 0 1 1 0 1 1 0 1 1 1;
       1 0 1 0 0 1 1 1 1 0 1 1;
       0 0 1 1 1 1 1 1 1 1 0 1]
    @test EDdeg(A) == 290
    @test EDdeg(A, true) == 630

    distance, x = EDdist(A, ones(12).*2)
    @test distance < 2

    #Octahedron example
    A = [1 1 1 0 0 0;
         1 0 0 1 1 0;
         0 1 0 1 0 1;
         0 0 1 0 1 1]
    @test EDdeg(A) == 28
    @test EDdeg(A, true) == 28
end

@testset "Graph_models" begin
    n_verts = 7
    edges = [[1,2], [1,6], [2,3], [3,4], [4,6], [5,6], [5,7], [6,7]]
    colour = [1,1,2,2,3,3,3]
    n_colours = 3
    G = Graph(n_verts, edges, colour, n_colours)
    A = homogenize(binom_mod(G))
    @test EDdeg(A) == 76
    B = homogenize(β_mod(G))
    @test EDdeg(B) == 1
end
