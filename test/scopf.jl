

using PowerModels

PowerModels.silence()

function import_data(file_name)
    data = PowerModels.parse_file(file_name)
    PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)
    return PowerModels.build_ref(data)[:it][:pm][:nw][0]
end

function case14_scopf_model()
    contingencies = [1, 6, 9, 18, 19, 21, 22]
    file_name = joinpath(@__DIR__, "case14.m")
    return MPCCBenchmark.scopf_model(
        import_data(file_name),
        contingencies;
    )
end

@testset "[PowerSystems] Test model $(JuMP.name(model))" for model in [
    case14_scopf_model(),
]
    ind_cc1, ind_cc2 = ComplementOpt.reformulate_to_vertical!(JuMP.backend(model))
    ComplementOpt.reformulate_as_nonlinear_program!(JuMP.backend(model); relaxation=1e-5)

    JuMP.set_optimizer(model, Ipopt.Optimizer)
    JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
    JuMP.set_silent(model)
    JuMP.optimize!(model)

    @test JuMP.is_solved_and_feasible(model)
end

