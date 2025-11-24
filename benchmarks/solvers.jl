
using Ipopt
using HSL_jll
using MadNLP
using MadNLPHSL
using MadMPEC
using NLPModels
using NLPModelsJuMP

abstract type AbstractSolverSetup end

#=
    Reference : Ipopt with JuMP
=#

@kwdef struct IpoptJuMP <: AbstractSolverSetup
    linear_solver::String = "ma27"
    max_iter::Int = 1000
    relaxation::Float64 = 1e-8
end

get_solver(solver::IpoptJuMP) = "ipopt-scholtes"
get_linear_solver(solver::IpoptJuMP) = solver.linear_solver
get_modeler(solver::IpoptJuMP) = "jump"

function solve_model(model, solver::IpoptJuMP)
    ind_cc1, ind_cc2 = ComplementOpt.reformulate_to_vertical!(JuMP.backend(model))
    ComplementOpt.reformulate_as_nonlinear_program!(JuMP.backend(model), ComplementOpt.ScholtesRelaxation(solver.relaxation))

    JuMP.set_optimizer(model, Ipopt.Optimizer)
    JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
    JuMP.set_optimizer_attribute(model, "max_iter", solver.max_iter)
    JuMP.set_optimizer_attribute(model, "linear_solver", solver.linear_solver)
    JuMP.set_optimizer_attribute(model, "hsllib", HSL_jll.libhsl)
    JuMP.set_optimizer_attribute(model, "bound_push", 1e-1)
    # JuMP.set_silent(model)
    JuMP.optimize!(model)

    return (
        JuMP.num_variables(model),
        JuMP.num_constraints(model; count_variable_in_set_constraints=false),
        length(ind_cc1),
        Int(JuMP.is_solved_and_feasible(model)),
        JuMP.objective_value(model),
        JuMP.barrier_iterations(model),
        JuMP.solve_time(model),
        get_complementarity_residual(model, ind_cc1, ind_cc2),
    )
end

#=
    MadNLPC
=#

@kwdef struct MadNLPCJuMP <: AbstractSolverSetup
    linear_solver = Ma27Solver
    max_iter::Int = 1000
end

get_solver(solver::MadNLPCJuMP) = "madnlpc"
get_linear_solver(solver::MadNLPCJuMP) = MadNLP.introduce(solver.linear_solver)
get_modeler(solver::MadNLPCJuMP) = "jump"


function solve_model(model, config::MadNLPCJuMP)
    ind_cc1, ind_cc2 = ComplementOpt.reformulate_to_vertical!(JuMP.backend(model))
    cc_cons = MOI.get(model, MOI.ListOfConstraintIndices{MOI.VectorOfVariables, MOI.Complements}())[1]
    MOI.delete(JuMP.backend(model), cc_cons)

    ind_x1 = getfield.(ind_cc1, :value)
    ind_x2 = getfield.(ind_cc2, :value)

    nlp = MathOptNLPModel(model)
    mpcc = MadMPEC.MPCCModelVarVar(nlp, ind_x1, ind_x2)

    madnlpc_opts = MadMPEC.MadNLPCOptions(;
        print_level=MadNLP.INFO,
        relaxation=MadMPEC.ScholtesRelaxation,
        use_magic_step=false,
        use_specialized_barrier_update=false,
    )
    solver = MadMPEC.MadNLPCSolver(
        mpcc;
        solver_opts=madnlpc_opts,
        print_level=MadNLP.INFO,
        bound_relax_factor=0.0,
        max_iter=config.max_iter,
        tol=1e-8,
        linear_solver=config.linear_solver,
    )
    stats = MadMPEC.solve_homotopy!(solver)

    cc_resid = get_complementarity_residual(nlp, stats.solution, ind_x1, ind_x2)
    return (
        NLPModels.get_nvar(nlp),
        NLPModels.get_ncon(nlp),
        length(ind_cc1),
        Int(stats.status),
        stats.stats.objective,
        stats.stats.iter,
        stats.stats.counters.total_time,
        cc_resid,
    )
end

