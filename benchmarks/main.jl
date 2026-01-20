
using DelimitedFiles
using MPCCBenchmark

include(joinpath(@__DIR__, "solvers.jl"))

const RESULTS_DIR = joinpath(@__DIR__, "..", "results")

if !isdir(RESULTS_DIR)
    mkdir(RESULTS_DIR)
end

# benchmark = MPCCBenchmark.PowerFlowBenchmark("/home/fpacaud/dev/matpower/data/")
benchmark = MPCCBenchmark.NOSNOCBenchmark(MPCCBenchmark.RadauIIA(2))
solver = MadNLPCJuMP()
results = MPCCBenchmark.run_benchmark(benchmark, solver)
dump_file = "$(MPCCBenchmark.introduce(benchmark, solver)).csv"
writedlm(joinpath(RESULTS_DIR, dump_file), results)

