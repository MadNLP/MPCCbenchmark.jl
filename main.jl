

using Distributed
using DelimitedFiles

@everywhere begin
    using MPCCBenchmark
    using MadNLPHSL
    include(joinpath(@__DIR__, "benchmarks", "solvers.jl"))
end

const RESULTS_DIR = joinpath(@__DIR__, "results")
const MATPOWER_DIR = "/home/anton/syscop/tools/matpower/data/"

const BENCHMARK_KEYS = Dict{String, MPCCBenchmark.AbstractBenchmarkSetting}(
    "pscc-pf" => MPCCBenchmark.PowerFlowBenchmark(MATPOWER_DIR),
    "pscc-scopf" => MPCCBenchmark.SCOPFBenchmark(MATPOWER_DIR),
    "nosnoc" => MPCCBenchmark.NOSNOCBenchmark(MPCCBenchmark.RadauIIA(2)),
)

const SOLVER_KEYS = Dict{String, MPCCBenchmark.AbstractSolverSetup}(
    "ipopt" => IpoptJuMP(),
    "madnlpc" => MadNLPCJuMP(max_iter=3000, linear_solver=Ma97Solver),#MadNLPCJuMP(linear_solver=Ma97Solver),
)

function parse_args(args::Vector{String})
    # Default options
    solver = SOLVER_KEYS["ipopt"]
    benchmark = BENCHMARK_KEYS["nosnoc"]
    for arg in args
        if startswith(arg, "--solver=")
            solver = SOLVER_KEYS[split(arg, "=")[2]]
        elseif startswith(arg, "--benchmark=")
            benchmark = BENCHMARK_KEYS[split(arg, "=")[2]]
        end
    end
    return solver, benchmark
end

function @main(args::Vector{String})
    solver, benchmark = parse_args(args)

    if !isdir(RESULTS_DIR)
        mkdir(RESULTS_DIR)
    end

    # Launch benchmark
    probs = MPCCBenchmark.get_instances(benchmark)
    retvals = pmap(prob->MPCCBenchmark.evaluate_model(benchmark, solver, prob), probs)

    # Export results
    n = length(probs)
    results = zeros(n, 7)
    for k in 1:n
        results[k, :] .= retvals[k]
    end
    names = MPCCBenchmark.get_name.(probs, Ref(benchmark))
    results = [reshape(MPCCBenchmark.COLS, 1, 8); names results]

    label = MPCCBenchmark.introduce(benchmark, solver)
    writedlm(joinpath(RESULTS_DIR, "$(label).csv"), results)
    return true
end

