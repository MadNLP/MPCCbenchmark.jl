
using Test
using MPCCBenchmark
using JuMP
using ComplementOpt
using Ipopt
using ExaPowerIO
using LinearAlgebra

using MadNLP, MadNLPHSL
using CCOpt
using NLPModelsJuMP

const MOIU = MOI.Utilities


contingencies = [1, 6, 9, 18, 19, 21, 22]
filename = joinpath(@__DIR__, "..", "case14.m")

data = ExaPowerIO.parse_matpower(filename)
MPCCBenchmark.compute_branch_limits!(data)

model = MPCCBenchmark.scopf_model(data, contingencies)

######
# Reformulate problem in standard form
######
ind_cc1, ind_cc2 = MPCCBenchmark.reformulate_to_standard_form!(JuMP.backend(model))

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
