
using JuMP
using MPCCBenchmark
using ComplementOpt
using MadNLP
using MadNLPHSL
using CCOpt
using NLPModelsJuMP
using Plots

######
# Import problem from MPCCBenchmark
######
N = 100
nfe = 2
collocation = MPCCBenchmark.RadauIIA(1)
model = MPCCBenchmark.nosnoc_schumacher_model(N, nfe, collocation; step_eq=:lcc, big_M=1e-2)

######
# Reformulate problem in vertical form with ComplementOpt
######
ind_cc1, ind_cc2 = reformulate_to_vertical!(model)

######
# Convert problem in CCOpt format
######
# Get evaluator as NLPModels
nlp = MathOptNLPModel(model)
# Build a MPCC problem
mpcc = CCOpt.MPCCModelVarVar(nlp, getfield.(ind_cc1, :value), getfield.(ind_cc2, :value))

######
# Solve problem with CCOpt
#####
madnlpc_opts = CCOpt.RelaxationOptions(;
    print_level=MadNLP.INFO,
    relaxation=CCOpt.ScholtesRelaxation,
    use_magic_step=false,
    use_specialized_barrier_update=false,
)
solver = CCOpt.RelaxationSolver(
    mpcc;
    solver_opts=madnlpc_opts,
    print_level=MadNLP.INFO,
    bound_relax_factor=0.0,
    tol=1e-8,
    linear_solver=Ma27Solver,
)
stats = CCOpt.solve_homotopy!(solver)

######
# Display solution
######
nh = N * nfe
qx = stats.solution[1:((nh+1)*3)][1:(nh+1)]
qy = stats.solution[((nh+1)*3+1):(6*(nh+1))][1:(nh+1)]

plot()
plot!(qx, qy, label="Trajectory")
xlabel!("qx")
ylabel!("qy")
xx = 0:0.1:3π
plot!(xx, MPCCBenchmark._schumacher_track.(xx) .+ 0.25, label="Track lb")
plot!(xx, MPCCBenchmark._schumacher_track.(xx) .- 0.25, label="Track ub")
