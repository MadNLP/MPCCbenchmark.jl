
include("utils.jl")

#=
    PowerFlowBenchmark
=#

struct PowerFlowBenchmark <: AbstractBenchmarkSetting
    path_instances::String
end

include("pf.jl")

get_tag(::PowerFlowBenchmark) = "pscc-pf"

function load_model(instance, bench::PowerFlowBenchmark)
    data = ExaPowerIO.parse_matpower(joinpath(bench.path_instances, instance))
    compute_branch_limits!(data)
    return powerflow_model(data)
end

function get_name(instance, ::PowerFlowBenchmark)
    return instance
end

function get_instances(::PowerFlowBenchmark)
    return [
        "case89pegase.m",
        "case118.m",
        "case_ACTIVSg200.m",
        "case300.m",
        "case_ACTIVSg500.m",
        "case1888rte.m",
        "case1951rte.m",
        "case_ACTIVSg2000.m",
        "case2848rte.m",
        "case2868rte.m",
        "case6468rte.m",
        "case6470rte.m",
        "case6495rte.m",
        "case6515rte.m",
        "case_ACTIVSg10k.m",
        "case_ACTIVSg25k.m",
    ]
end


#=
    SCOPFBenchmark
=#

struct SCOPFBenchmark <: AbstractBenchmarkSetting
    path_instances::String
end

include("scopf.jl")

get_tag(::SCOPFBenchmark) = "pscc-scopf"

function load_model(instance, bench::SCOPFBenchmark)
    N = instance[2]
    data = ExaPowerIO.parse_matpower(joinpath(bench.path_instances, "$(instance[1]).m"))
    compute_branch_limits!(data)
    contingencies = CONTINGENCIES[instance[1]][1:N]
    return scopf_model(data, contingencies)
end

function get_name(instance, ::SCOPFBenchmark)
    return "$(instance[1])_$(instance[2])"
end

function get_instances(::SCOPFBenchmark)
    return [
        ("case118", 10),
        ("case118", 100),
        ("case300", 10),
        ("case300", 50),
        ("case_ACTIVSg200", 10),
        ("case_ACTIVSg200", 50),
        ("case_ACTIVSg200", 100),
        ("case_ACTIVSg500", 10),
        ("case_ACTIVSg500", 50),
        ("case_ACTIVSg500", 100),
        ("case1354pegase", 8),
        ("case1354pegase", 16),
        ("case1354pegase", 32),
        ("case_ACTIVSg2000", 8),
        ("case_ACTIVSg2000", 16),
        ("case2869pegase", 8),
        ("case2869pegase", 16),
    ]
end

