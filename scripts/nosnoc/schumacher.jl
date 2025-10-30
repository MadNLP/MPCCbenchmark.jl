
using JuMP
using MPCCBenchmark
using ComplementOpt
using Ipopt
using Plots


######
# Import problem from MPCCBenchmark
######
N = 100
nh = 3
collocation = MPCCBenchmark.CrankNicolson()
model = MPCCBenchmark.nosnoc_schumacher_model(N, nh, collocation; step_eq=:heuristic_mean)

######
# Reformulate problem as a nonlinear program with ComplementOpt
######
ind_cc1, ind_cc2 = ComplementOpt.reformulate_to_vertical!(JuMP.backend(model))
ComplementOpt.reformulate_as_nonlinear_program!(JuMP.backend(model); relaxation=1e-5)

######
# Solve problem with Ipopt
#####
JuMP.set_optimizer(model, Ipopt.Optimizer)
JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
JuMP.optimize!(model)

######
# Display solution
######
qx = JuMP.value.(model[:qx])[:, 1]
qy = JuMP.value.(model[:qy])[:, 1]

plot()
plot!(qx, qy, label="Trajectory")
xlabel!("qx")
ylabel!("qy")
xx = 0:0.1:3π
plot!(xx, MPCCBenchmark._schumacher_track.(xx) .+ 0.25, label="Track lb")
plot!(xx, MPCCBenchmark._schumacher_track.(xx) .- 0.25, label="Track ub")
