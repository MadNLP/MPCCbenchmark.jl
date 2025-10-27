
#=
    Collocation scheme.
=#

"""
    RKScheme

Store Butcher tableau used in the collocation.

"""
struct RKScheme
    c::Vector{Float64}
    b::Vector{Float64}
    a::Matrix{Float64}
end

function GaussLegendre()
    p = √3/6
    return RKScheme([0.5 - p, 0.5 + p], [0.5, 0.5], [0.25 0.25-p; 0.25+p 0.25])
end
ImplicitEuler() = RKScheme([1.0], [1.0], ones(1, 1))
CrankNicolson() = RKScheme([0.0, 1.0], [0.5, 0.5], [0.0 0.0; 0.5 0.5])

