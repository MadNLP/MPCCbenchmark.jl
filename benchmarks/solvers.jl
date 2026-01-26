
using JuMP
using ComplementOpt
using Ipopt
using HSL_jll
using MadNLP
using MadNLPHSL
using MadMPEC
using NLPModels
using NLPModelsJuMP
using LinearAlgebra

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
    MadNLPC
=#

@kwdef struct MadNLPCJuMP <: MPCCBenchmark.AbstractSolverSetup
    linear_solver = Ma57Solver
    max_iter::Int = 3000
end

MPCCBenchmark.get_solver(solver::MadNLPCJuMP) = "madnlpc"

function MPCCBenchmark.solve_model(config::MadNLPCJuMP, model)
    ind_cc1, ind_cc2 = MPCCBenchmark.reformulate_to_vertical!(model)

    println(ind_cc1)
    ind_x1 = getfield.(ind_cc1, :value)
    ind_x2 = getfield.(ind_cc2, :value)
    println(ind_x1)

    nlp = MathOptNLPModel(model)
    mpcc = MadMPEC.MPCCModelVarVar(nlp, ind_x1, ind_x2)

    println(mpcc.meta.ind_cc1)
    madnlpc_opts = MadMPEC.MadNLPCOptions(
        ;
        print_level=MadNLP.INFO,
        relaxation=MadMPEC.ScholtesRelaxation,
        #relaxation_update=MadMPEC.RelaxLBUpdate(),
        use_magic_step=false,
        use_specialized_barrier_update=true,
        center_complementarities=true,
    )
    solver = MadMPEC.MadNLPCSolver(
        mpcc;
        solver_opts=madnlpc_opts,
        print_level=MadNLP.ERROR,
        bound_relax_factor=0.0,
        max_iter=config.max_iter,
        tol=1e-8,
        linear_solver=config.linear_solver,
        barrier=MadNLP.QualityFunctionUpdate(mu_max = 1.0, max_gs_iter=12),
    )
    stats = MadMPEC.solve_homotopy!(solver)
    # TODO: fix CC resid
    cc_resid = get_complementarity_residual(nlp, stats.solution, ind_x1, ind_x2)
    println("status = $(stats.status)")
    return (
       \ NLPModels.get_nvar(nlp),
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

@kwdef struct MadNLPHomotopyJump <: MPCCBenchmark.AbstractSolverSetup
    linear_solver = Ma27Solver
    max_iter::Int = 1000
end

MPCCBenchmark.get_solver(solver::MadNLPCJuMP) = "madnlpc"

function MPCCBenchmark.solve_model(config::MadNLPCJuMP, model)
    ind_cc1, ind_cc2 = MPCCBenchmark.reformulate_to_vertical!(model)

    println(ind_cc1)
    ind_x1 = getfield.(ind_cc1, :value)
    ind_x2 = getfield.(ind_cc2, :value)
    println(ind_x1)

    nlp = MathOptNLPModel(model)
    mpcc = MadMPEC.MPCCModelVarVar(nlp, ind_x1, ind_x2)

    println(mpcc.meta.ind_cc1)
    madnlpc_opts = MadMPEC.MadNLPCOptions(
        ;
        print_level=MadNLP.INFO,
        relaxation=MadMPEC.ScholtesRelaxation,
        #relaxation_update=MadMPEC.RelaxLBUpdate(),
        use_magic_step=false,
        use_specialized_barrier_update=true,
        #center_complementarities=true,
    )
    solver = MadMPEC.MadNLPCSolver(
        mpcc;
        solver_opts=madnlpc_opts,
        print_level=MadNLP.ERROR,
        bound_relax_factor=0.0,
        max_iter=config.max_iter,
        tol=1e-8,
        linear_solver=config.linear_solver,
        #barrier=MadNLP.MonotoneUpdate(mu_init=1.0),
        #barrier=MadNLP.QualityFunctionUpdate(mu_max = 1.0, max_gs_iter=12),
    )
    stats = MadMPEC.solve_homotopy!(solver)
    # TODO: fix CC resid
    cc_resid = get_complementarity_residual(nlp, stats.solution, ind_x1, ind_x2)
    println("status = $(stats.status)")
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
