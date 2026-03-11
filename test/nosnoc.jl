
@testset "[NOSNOC] Test problem $(JuMP.name(model))" for model in [
    MPCCBenchmark.nosnoc_stewart_anitescu_model(10, MPCCBenchmark.RadauIIA(1)),
]
    JuMP.set_optimizer(model, () -> ComplementOpt.Optimizer(Ipopt.Optimizer()))
    MOI.set(model, ComplementOpt.RelaxationMethod(), ComplementOpt.ScholtesRelaxation(1e-5))
    JuMP.set_silent(model)
    JuMP.optimize!(model)

    @test JuMP.is_solved_and_feasible(model)
end

