
using JuMP
using ComplementOpt
using Ipopt
using HSL_jll
using MadNLP
using MadNLPHSL
using MadNCL
using CCOpt
using NLPModels
using NLPModelsIpopt
using ExaModels
using LinearAlgebra

# Dirty hacks for ExaModels
NLPModels.jac_nln_structure!(nlp::ExaModels.ExaModel, rows, cols) = NLPModels.jac_structure!(nlp, rows, cols)
NLPModels.jac_nln_coord!(nlp::ExaModels.ExaModel, x, jac) = NLPModels.jac_coord!(nlp, x, jac)
NLPModels.jac_lin_coord!(nlp::ExaModels.ExaModel, x, jac) = nothing

function get_complementarity_residual(model::JuMP.Model, ind_cc1, ind_cc2)
    moi_model = JuMP.backend(model)
    x1 = MOI.get.(moi_model, MOI.VariablePrimal(), ind_cc1)
    x2 = MOI.get.(moi_model, MOI.VariablePrimal(), ind_cc2)

    bounds = MOI.Utilities.get_bounds.(moi_model, Float64, ind_cc2)
    lb = [b[1] for b in bounds]
    ub = [b[2] for b in bounds]

    resid = max.(min.(x1, x2 .- lb), x2 .- ub)
    return norm(resid, Inf)
end

function get_complementarity_residual(nlp::AbstractNLPModel, solution::AbstractVector, ind_cc1, ind_cc2)
    lb = NLPModels.get_lvar(nlp)
    ub = NLPModels.get_uvar(nlp)
    x1 = solution[ind_cc1]
    x2 = solution[ind_cc2]
    resid = min.(x1 .- lb[ind_cc1], x2 .- lb[ind_cc2])
    return norm(resid, Inf)
end


#=
    Reference : Ipopt with JuMP
=#

@kwdef struct IpoptJuMP <: MPCCBenchmark.AbstractSolverSetup
    linear_solver::String = "ma57"
    max_iter::Int = 3000
    relaxation::Float64 = 1e-8
end

MPCCBenchmark.get_solver(solver::IpoptJuMP) = "ipopt-scholtes"

function MPCCBenchmark.solve_model(solver::IpoptJuMP, model)
    JuMP.set_optimizer(model, () -> ComplementOpt.Optimizer(Ipopt.Optimizer()))
    MOI.set(model, ComplementOpt.RelaxationMethod(), ComplementOpt.ScholtesRelaxation(solver.relaxation))
    JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
    JuMP.set_optimizer_attribute(model, "max_iter", solver.max_iter)
    JuMP.set_optimizer_attribute(model, "linear_solver", solver.linear_solver)
    JuMP.set_optimizer_attribute(model, "hsllib", HSL_jll.libhsl)
    JuMP.set_optimizer_attribute(model, "bound_push", 1e-1)
    JuMP.set_silent(model)
    JuMP.optimize!(model)

    return (
        JuMP.num_variables(model),
        JuMP.num_constraints(model; count_variable_in_set_constraints=false),
        MOI.get(model, MOI.NumberOfConstraints{MOI.VectorOfVariables,MOI.Complements}()),
        Int(JuMP.is_solved_and_feasible(model)),
        JuMP.objective_value(model),
        JuMP.barrier_iterations(model),
        JuMP.solve_time(model),
    )
end

#=
    CCOpt Relaxation solver
=#

@kwdef struct CCOptRelaxation <: MPCCBenchmark.AbstractSolverSetup
    linear_solver = Ma27Solver
    max_iter::Int = 3000
    tol::Float64 = 1e-8
end

MPCCBenchmark.get_solver(solver::CCOptRelaxation) = "ccopt-relaxation"

function MPCCBenchmark.solve_model(config::CCOptRelaxation, model)
    model = MPCCBenchmark.reformulate_to_vertical!(JuMP.backend(model))
    ind_cc1, ind_cc2 = MPCCBenchmark.reformulate_to_standard_form!(model)
    ind_x1 = getfield.(ind_cc1, :value)
    ind_x2 = getfield.(ind_cc2, :value)

    nlp = ExaModel(model)
    mpcc = CCOpt.MPCCModelVarVar(nlp, ind_x1, ind_x2)

    madnlpc_opts = CCOpt.RelaxationOptions(
        ;
        print_level=MadNLP.INFO,
        relaxation=CCOpt.ScholtesRelaxation,
		relaxation_update=CCOpt.RolloffRelaxationUpdate(),
        use_magic_step=false,
        use_specialized_barrier_update=false,
    )
    solver = CCOpt.RelaxationSolver(
        mpcc;
        solver_opts=madnlpc_opts,
        print_level=MadNLP.ERROR,
        bound_relax_factor=0.0,
        max_iter=config.max_iter,
        tol=config.tol,
        linear_solver=config.linear_solver,
        # ma57_automatic_scaling=true,
        # barrier=MadNLP.QualityFunctionUpdate(mu_max = 1.0, max_gs_iter=12),
    )
    stats = CCOpt.solve_homotopy!(solver)
    # TODO: fix CC resid
    cc_resid = get_complementarity_residual(nlp, stats.solution, ind_x1, ind_x2)
    return (
        NLPModels.get_nvar(nlp),
        NLPModels.get_ncon(nlp),
        length(ind_cc1),
        Int(stats.status),
        stats.objective,
        stats.iter,
        stats.counters.counters.total_time,
    )
end

#=
    MadNLP homotopy
=#

@kwdef struct MadNLPHomotopy <: MPCCBenchmark.AbstractSolverSetup
    linear_solver = Ma57Solver
    max_iter::Int = 3000
end

MPCCBenchmark.get_solver(solver::MadNLPHomotopy) = "ccopt-homotopy"

function MPCCBenchmark.solve_model(config::MadNLPHomotopy, model)
    model = MPCCBenchmark.reformulate_to_vertical!(JuMP.backend(model))
    ind_cc1, ind_cc2 = MPCCBenchmark.reformulate_to_standard_form!(model)

    ind_x1 = getfield.(ind_cc1, :value)
    ind_x2 = getfield.(ind_cc2, :value)

    nlp = ExaModel(model)
    mpcc = CCOpt.MPCCModelVarVar(nlp, ind_x1, ind_x2)

    # homotopy_opts = CCOpt.HomotopySolverOptions(max_inner_iter=config.max_iter)
    # homotopy_opts.nlp_solver_options = Dict(:bound_relax_factor=>0.0,
    #                                         :print_level=>MadNLP.ERROR,
    #                                         :linear_solver=>config.linear_solver,
    #                                         :max_iter=>config.max_iter,
    #                                         :barrier=>MadNLP.QualityFunctionUpdate())
    # solver = CCOpt.HomotopySolver(mpcc, MadNLP.MadNLPSolver, homotopy_opts)
	homotopy_opts = CCOpt.HomotopySolverOptions(max_inner_iter=1000)
	solver = CCOpt.HomotopySolver(mpcc, NLPModelsIpopt.IpoptSolver, homotopy_opts)

    stats = CCOpt.solve!(solver)
    # TODO: fix CC resid
    cc_resid = get_complementarity_residual(nlp, stats.solution, ind_x1, ind_x2)
    return (
        NLPModels.get_nvar(nlp),
        NLPModels.get_ncon(nlp),
        length(ind_cc1),
        Int(stats.status),
        stats.objective,
        stats.iter,
        stats.wall_time,
    )
end

#=
    MadNCL
=#

@kwdef struct MadNCLSolver <: MPCCBenchmark.AbstractSolverSetup
    linear_solver = Ma57Solver
    max_iter::Int = 3000
    tol::Float64 = 1e-5
end

MPCCBenchmark.get_solver(solver::MadNCLSolver) = "madncl"

function MPCCBenchmark.solve_model(config::MadNCLSolver, model)
    MPCCBenchmark.reformulate_to_vertical!(JuMP.backend(model))
    ind_cc1, ind_cc2 = MPCCBenchmark.reformulate_to_standard_form!(JuMP.backend(model))
    ncc = length(ind_cc1)
    MOI.add_constraint(JuMP.backend(model), [ind_cc1; ind_cc2] , MOI.Complements(ncc))
    MPCCBenchmark.reformulate_to_nonlinear!(JuMP.backend(model), ComplementOpt.ScholtesRelaxation(config.tol))

    nlp = ExaModel(model)
    ncl_options = MadNCL.NCLOptions{Float64}(;
        opt_tol=config.tol,
        feas_tol=config.tol,
        scaling=true,
        scaling_max_gradient=100.0,
        extrapolation=true,
        verbose=true,
    )

    stats = MadNCL.madncl(
        nlp;
        ncl_options=ncl_options,
        linear_solver=Ma57Solver,
        print_level=MadNLP.ERROR,
        richardson_tol=1e-12,
        richardson_max_iter=20,
        max_iter=1000,
        ma57_automatic_scaling=true,
        kkt_system=MadNCL.K2rAuglagKKTSystem,
    )

    # cc_resid = get_complementarity_residual(nlp, stats.solution, ind_x1, ind_x2)
    return (
        NLPModels.get_nvar(nlp),
        NLPModels.get_ncon(nlp),
        length(ind_cc1),
        Int(stats.status),
        stats.objective,
        stats.iter,
        stats.counters.total_time,
    )
end

