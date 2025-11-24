
using DelimitedFiles
using LinearAlgebra
using MPCCBenchmark
using JuMP
using ComplementOpt

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "setup.jl"))
include(joinpath(@__DIR__, "solvers.jl"))

const RESULTS_DIR = joinpath(@__DIR__, "..", "results")

const COLS = String[
    "instance",
    "n",
    "m",
    "ncc",
    "status",
    "obj",
    "it",
    "solve_time",
    "cc_resid",
]

function introduce(benchmark::AbstractBenchmarkSetting, solver::AbstractSolverSetup)
    return "$(get_tag(benchmark))-$(get_solver(solver))-$(get_linear_solver(solver))-$(get_modeler(solver))"
end

function run_benchmark(benchmark, solver)
    instances = get_instances(benchmark)
    m, n = length(instances), 8
    results = zeros(m, n)
    for (k, instance) in enumerate(instances)
        model = load_model(instance, benchmark)
        results[k, :] .= solve_model(model, solver)
    end
    names = parse_name.(instances, Ref(benchmark))
    return [reshape(COLS, 1, n+1); names results]
end

if !isdir(RESULTS_DIR)
    mkdir(RESULTS_DIR)
end
benchmark = NOSNOC()
solver = IpoptJuMP()
results = run_benchmark(benchmark, solver)
dump_file = "$(introduce(benchmark, solver)).csv"
writedlm(joinpath(RESULTS_DIR, dump_file), results)

