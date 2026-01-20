
using Test
using MPCCBenchmark
using JuMP
using ComplementOpt
using Ipopt
using ExaPowerIO

function import_data(file_name)
    data = ExaPowerIO.parse_matpower(file_name)
    MPCCBenchmark.compute_branch_limits!(data)
    return data
end

function case14_scopf_model()
    contingencies = [1, 6, 9, 18, 19, 21, 22]
    file_name = joinpath(@__DIR__, "case14.m")
    return MPCCBenchmark.scopf_model(
        import_data(file_name),
        contingencies;
    )
end

function case14_pf_model()
    file_name = joinpath(@__DIR__, "case14.m")
    return MPCCBenchmark.powerflow_model(
        import_data(file_name),
    )
end

@testset "[PowerSystems] Test model $(JuMP.name(model))" for model in [
    case14_scopf_model(),
    case14_pf_model(),
]
    JuMP.set_optimizer(model, () -> ComplementOpt.Optimizer(Ipopt.Optimizer()))
    MOI.set(model, ComplementOpt.RelaxationMethod(), ComplementOpt.ScholtesRelaxation(1e-7))
    JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
    JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
    JuMP.optimize!(model)

    @test JuMP.is_solved_and_feasible(model)
end

