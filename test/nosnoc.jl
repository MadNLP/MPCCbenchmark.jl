
@testset "[NOSNOC] Test problem $(JuMP.name(model))" for model in [
    MPCCBenchmark.nosnoc_stewart_anitescu_model(10, MPCCBenchmark.CrankNicolson()),
    MPCCBenchmark.nosnoc_schumacher_model(10, 3, MPCCBenchmark.CrankNicolson(); step_eq=:heuristic_mean),
    MPCCBenchmark.nosnoc_schumacher_model(10, 3, MPCCBenchmark.CrankNicolson(); step_eq=:lcc)
]
    ind_cc1, ind_cc2 = ComplementOpt.reformulate_to_vertical!(JuMP.backend(model))
    ComplementOpt.reformulate_as_nonlinear_program!(JuMP.backend(model); relaxation=1e-5)

    JuMP.set_optimizer(model, Ipopt.Optimizer)
    JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
    JuMP.set_silent(model)
    JuMP.optimize!(model)

    @test JuMP.is_solved_and_feasible(model)
end

