using JuMP
using MPCCBenchmark
using ComplementOpt
using MadNLP
using MadNLPHSL
using MadMPEC
using NLPModelsJuMP
using Plots


N = 57
nfe = 2
collocation = MPCCBenchmark.RadauIIA(2)
######
# Import problem from MPCCBenchmark
######
model = MPCCBenchmark.nosnoc_cart_pole_with_friction_model(N, nfe, collocation; step_eq=:heuristic_mean)

######
# Reformulate problem in vertical form with ComplementOpt
######
ind_cc1, ind_cc2 = MPCCBenchmark.reformulate_to_vertical!(model)

######
# Convert problem in MadMPEC format
######
# Get evaluator as NLPModels
nlp = MathOptNLPModel(model)
# Build a MPCC problem
mpcc = MadMPEC.MPCCModelVarVar(nlp, getfield.(ind_cc1, :value), getfield.(ind_cc2, :value))

######
# Solve problem with MadMPEC
#####
madnlpc_opts = MadMPEC.MadNLPCOptions(
    ;
    print_level=MadNLP.INFO,
    relaxation=MadMPEC.ScholtesRelaxation,
    use_magic_step=false,
    use_specialized_barrier_update=false,
    center_complementarities=false,
    relaxation_update=MadMPEC.ProportionalRelaxationUpdate(sigma_mu_ratio=1.0),
)
solver = MadMPEC.MadNLPCSolver(
    mpcc;
    solver_opts=madnlpc_opts,
    print_level=MadNLP.INFO,
    bound_relax_factor=0.0,
    tol=1e-8,
    linear_solver=Ma97Solver,
    #barrier=MadNLP.MonotoneUpdate(mu_init=1.0),
    barrier=MadNLP.QualityFunctionUpdate(mu_max = 1.0, max_gs_iter=12),
)
stats = MadMPEC.solve_homotopy!(solver)

println("Termination status: ", stats.status)
println("Objective value: ", stats.objective)

## again
######
# Import problem from MPCCBenchmark
######
model = MPCCBenchmark.nosnoc_cart_pole_with_friction_model(N, nfe, collocation; step_eq=:heuristic_mean)

######
# Reformulate problem in vertical form with ComplementOpt
######
ind_cc1, ind_cc2 = MPCCBenchmark.reformulate_to_vertical!(model)

######
# Convert problem in MadMPEC format
######
# Get evaluator as NLPModels
nlp = MathOptNLPModel(model)
# Build a MPCC problem
mpcc = MadMPEC.MPCCModelVarVar(nlp, getfield.(ind_cc1, :value), getfield.(ind_cc2, :value))

######
# Solve problem with MadMPEC
#####
madnlpc_opts = MadMPEC.MadNLPCOptions(
    ;
    print_level=MadNLP.INFO,
    relaxation=MadMPEC.ScholtesRelaxation,
    use_magic_step=false,
    use_specialized_barrier_update=false,
    center_complementarities=false,
    relaxation_update=MadMPEC.ProportionalRelaxationUpdate(sigma_mu_ratio=1.0),
)
solver = MadMPEC.MadNLPCSolver(
    mpcc;
    solver_opts=madnlpc_opts,
    print_level=MadNLP.INFO,
    bound_relax_factor=0.0,
    tol=1e-8,
    linear_solver=Ma97Solver,
    #barrier=MadNLP.MonotoneUpdate(mu_init=1.0),
    barrier=MadNLP.QualityFunctionUpdate(mu_max = 1.0, max_gs_iter=12),
)
stats = MadMPEC.solve_homotopy!(solver)

println("Termination status: ", stats.status)
println("Objective value: ", stats.objective)
