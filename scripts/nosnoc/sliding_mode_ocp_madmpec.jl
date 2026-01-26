using JuMP
using MPCCBenchmark
using ComplementOpt
using MadNLP
using MadNLPHSL
using MadMPEC
using NLPModelsJuMP
using Plots

######
# Import problem from MPCCBenchmark
######
N = 60
nfe = 3
collocation = MPCCBenchmark.ImplicitEuler()
model = MPCCBenchmark.nosnoc_sliding_mode_ocp_model(N, nfe, collocation; step_eq=:heuristic_mean, rho_h=1e4)

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
    center_complementarities=true
)
solver = MadMPEC.MadNLPCSolver(
    mpcc;
    solver_opts=madnlpc_opts,
    print_level=MadNLP.INFO,
    bound_relax_factor=0.0,
    tol=1e-8,
    linear_solver=Ma97Solver,
)
stats = MadMPEC.solve_homotopy!(solver)

println("Termination status: ", stats.status)
println("Objective value: ", stats.objective)
