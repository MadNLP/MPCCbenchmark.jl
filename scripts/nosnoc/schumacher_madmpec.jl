
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
N = 100
nfe = 3
collocation = MPCCBenchmark.CrankNicolson()
model = MPCCBenchmark.nosnoc_schumacher_model(N, nfe, collocation; step_eq=:heuristic_mean)

######
# Reformulate problem in vertical form with ComplementOpt
######
ind_cc1, ind_cc2 = ComplementOpt.reformulate_to_vertical!(JuMP.backend(model))

######
# Convert problem in MadMPEC format
######
# Remove complementarity constraints
cc_cons = MOI.get(model, MOI.ListOfConstraintIndices{MOI.VectorOfVariables, MOI.Complements}())[1]
MOI.delete(JuMP.backend(model), cc_cons)
# Get evaluator as NLPModels
nlp = MathOptNLPModel(model)
# Build a MPCC problem
mpcc = MadMPEC.MPCCModelVarVar(nlp, getfield.(ind_cc1, :value), getfield.(ind_cc2, :value))

######
# Solve problem with MadMPEC
#####
madnlpc_opts = MadMPEC.MadNLPCOptions(;
    print_level=MadNLP.INFO,
    use_mpecopt=false,
    phase_I_oracle=:naive,
    eps_proj=1e-6,
    relaxation=MadMPEC.ScholtesRelaxation,
    use_magic_step=false,
    use_specialized_barrier_update=false,
)
solver = MadMPEC.MadNLPCSolver(
    mpcc;
    solver_opts=madnlpc_opts,
    print_level=MadNLP.INFO,
    bound_relax_factor=0.0,
    tol=1e-8,
    linear_solver=Ma27Solver,
)
stats = MadMPEC.solve_homotopy!(solver)

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
