
using Test
using JuMP
using MPCCBenchmark
using ComplementOpt
using Ipopt

@testset "NOSNOC" begin
    include("nosnoc.jl")
end

@testset "Power systems" begin
    include("powersystems.jl")
end

