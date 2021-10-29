
using Test






@testset "utilities" begin
    A=[0 1 0 0 1 1 0 1 1 1 1 1;
       0 1 1 1 0 0 1 0 1 1 1 1;
       1 0 0 1 1 0 1 1 0 1 1 1;
       1 0 1 0 0 1 1 1 1 0 1 1;
       0 0 1 1 1 1 1 1 1 1 0 1]

    @test ishomogeneous(A) == false
    @test ishomogeneous(homogenize(A)) == true

    B = reparametrise(A)
    @test rank(A) == rank(B)
end



@testset "EDdistance" begin
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
    B = [1 1 1 0 0 0;
         1 0 0 1 1 0;
         0 1 0 1 0 1;
         0 0 1 0 1 1]
    @test EDdeg(B) == 28
    @test EDdeg(B, true) == 28
end

@testset "Graph_model" begin
    n_verts = 7
    edges = [[1,2], [1,6], [2,3], [3,4], [4,6], [5,6], [5,7], [6,7]]
    colour = [1,1,2,2,3,3,3]
    n_colours = 3
    G = Graph(n_verts, edges, colour, n_colours)
    A = binom_mod(G)
    #@test EDdeg(A) == 76
    B = β_mod(G)
    #@test EDdeg(B) == 1
end
