
using Test
using MPCCBenchmark
using JuMP
using ComplementOpt
using Ipopt

@testset "[EPEC] Test model $(JuMP.name(model))" for model in [
    MPCCBenchmark.macepec_ex001_model(),
    MPCCBenchmark.macepec_ex4_model(),
    MPCCBenchmark.macepec_outrata3_model(),
    MPCCBenchmark.macepec_outrata4_model(),
    MPCCBenchmark.macepec_electric002_model(),
    # MPCCBenchmark.macepec_electric004_model(), # this instance does not converge with Ipopt
]
    ind_cc1, ind_cc2 = ComplementOpt.reformulate_to_vertical!(JuMP.backend(model))
    ComplementOpt.reformulate_as_nonlinear_program!(JuMP.backend(model), ComplementOpt.ScholtesRelaxation(1e-5))

    JuMP.set_optimizer(model, Ipopt.Optimizer)
    JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
    # Use very aggressive bound_push for better convergence
    JuMP.set_optimizer_attribute(model, "bound_push", 1e0)
    JuMP.set_silent(model)
    JuMP.optimize!(model)

    @test JuMP.is_solved_and_feasible(model)
end
