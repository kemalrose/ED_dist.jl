

module ED_dist

using HomotopyContinuation, AbstractAlgebra, LinearAlgebra
export EDdist, EDdeg, to_system, Graph, β_mod, binom_mod, criticals, homogenize, criticals_on_X, index,
reparametrise, ishomogeneous, make_positive, new_rows, edge_list, Graph


include("EDdistance.jl")
include("Graph_model.jl")
include("utilities.jl")
end
