
#=
    Use same setup as in MadNLPBenchmark.
=#

abstract type AbstractBenchmarkSetting end
function get_tag end
function load_model end
function get_name end
function get_instances end

abstract type AbstractSolverSetup end
function get_solver end
function solve_model end

const COLS = String[
    "instance",
    "n",
    "m",
    "ncc",
    "status",
    "obj",
    "it",
    "solve_time",
]

function introduce(benchmark::AbstractBenchmarkSetting, solver::AbstractSolverSetup)
    return "$(get_tag(benchmark))-$(get_solver(solver))"
end

function evaluate_model(
    benchmark::AbstractBenchmarkSetting,
    solver::AbstractSolverSetup,
    instance;
    options...,
)
    @info instance
    model = load_model(instance, benchmark)
    results = solve_model(solver, model; options...)
    # N.B. *ALWAYS* call explicitly garbage collector
    #      to avoid annoying memory leak as GC was disable before.
    GC.gc(true)
    return results
end

function run_benchmark(benchmark::AbstractBenchmarkSetting, solver::AbstractSolverSetup; options...)
    instances = get_instances(benchmark)
    m, n = length(instances), 7
    results = zeros(m, n)
    for (k, instance) in enumerate(instances)
        results[k, :] .= evaluate_model(benchmark, solver, instance; options...)
    end
    names = get_name.(instances, Ref(benchmark))
    return [reshape(COLS, 1, n+1); names results]
end

