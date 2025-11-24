
using MPCCBenchmark
using ExaPowerIO

const MATPOWER_DIR = "/home/fpacaud/dev/matpower/data"

abstract type AbstractBenchmarkSetting end

#=
    Power-flow benchmark.
=#

struct PowerFlow <: AbstractBenchmarkSetting end

get_tag(::PowerFlow) = "powerflow"

function get_instances(::PowerFlow)
    # N.B. pegase instances are commented out as some reactive power generation
    #      are unbounded, leading to a parsing error in ComplementOpt.
    return [
        "case89pegase.m",
        "case118.m",
        "case_ACTIVSg200.m",
        "case300.m",
        "case_ACTIVSg500.m",
        # "case1354pegase.m",
        "case1888rte.m",
        "case1951rte.m",
        "case_ACTIVSg2000.m",
        "case2848rte.m",
        "case2868rte.m",
        # "case2869pegase.m",
        "case6468rte.m",
        "case6470rte.m",
        "case6495rte.m",
        "case6515rte.m",
        # "case8387pegase.m",
        # "case9241pegase.m",
        "case_ACTIVSg10k.m",
        # "case13659pegase.m",
        "case_ACTIVSg25k.m",
    ]
end

function parse_name(instance, ::PowerFlow)
    return split(instance, ".")[1]
end

function load_model(instance, ::PowerFlow)
    file_name = joinpath(MATPOWER_DIR, instance)
    data = ExaPowerIO.parse_matpower(file_name)
    return MPCCBenchmark.powerflow_model(data)
end

#=
    SCOPF benchmark.
=#

struct SCOPF <: AbstractBenchmarkSetting end

get_tag(::SCOPF) = "scopf"

include(joinpath(@__DIR__, "data", "contingencies.jl"))

function get_instances(::SCOPF)
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

function parse_name(instance, ::SCOPF)
    return "$(instance[1])_$(instance[2])"
end

function load_model(instance, ::SCOPF)
    nK = instance[2]
    file_name = joinpath(MATPOWER_DIR, "$(instance[1]).m")
    data = ExaPowerIO.parse_matpower(file_name)
    MPCCBenchmark.compute_branch_limits!(data)
    contingencies = CONTINGENCIES[instance[1]][1:nK]
    return MPCCBenchmark.scopf_model(
        data, contingencies,
    )
end

#=
    NOSNOC benchmark.
=#

struct NOSNOC <: AbstractBenchmarkSetting end

get_tag(::NOSNOC) = "nosnoc"

function get_instances(::NOSNOC)
    return [
        (MPCCBenchmark.nosnoc_sliding_mode_ocp_model, MPCCBenchmark.RadauIIA(3), 6, 3)
        (MPCCBenchmark.nosnoc_sliding_mode_ocp_model, MPCCBenchmark.RadauIIA(3), 37, 3)
        (MPCCBenchmark.nosnoc_sliding_mode_ocp_model, MPCCBenchmark.RadauIIA(3), 50, 3)
        (MPCCBenchmark.nosnoc_sliding_mode_ocp_model, MPCCBenchmark.RadauIIA(3), 63, 3)
        (MPCCBenchmark.nosnoc_schumacher_model, MPCCBenchmark.RadauIIA(4), 50, 2)
        (MPCCBenchmark.nosnoc_schumacher_model, MPCCBenchmark.RadauIIA(4), 80, 2)
        (MPCCBenchmark.nosnoc_motor_with_friction_model, MPCCBenchmark.RadauIIA(2), 27, 3)
        (MPCCBenchmark.nosnoc_motor_with_friction_model, MPCCBenchmark.RadauIIA(2), 30, 3)
        (MPCCBenchmark.nosnoc_motor_with_friction_model, MPCCBenchmark.RadauIIA(2), 50, 3)
        (MPCCBenchmark.nosnoc_motor_with_friction_model, MPCCBenchmark.RadauIIA(2), 53, 3)
        (MPCCBenchmark.nosnoc_cart_pole_with_friction_model, MPCCBenchmark.RadauIIA(2), 30, 2)
        (MPCCBenchmark.nosnoc_cart_pole_with_friction_model, MPCCBenchmark.RadauIIA(2), 33, 2)
        (MPCCBenchmark.nosnoc_cart_pole_with_friction_model, MPCCBenchmark.RadauIIA(2), 50, 2)
        (MPCCBenchmark.nosnoc_cart_pole_with_friction_model, MPCCBenchmark.RadauIIA(2), 57, 2)
    ]
end

function parse_name(instance, ::NOSNOC)
    s = string(instance[1])
    name = match(r"^nosnoc_(.*)_model$", s).captures[1]
    return "$(name)_$(instance[3])_$(instance[4])"
end

function load_model(instance, ::NOSNOC)
    scheme = instance[2]
    N, nfe = instance[3], instance[4]
    return instance[1](N, nfe, scheme)
end

