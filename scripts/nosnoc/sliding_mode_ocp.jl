
using JuMP
using MPCCBenchmark
using ComplementOpt
using Ipopt
using Plots

######
# Import problem from MPCCBenchmark
######
N = 60
nfe = 3
collocation = MPCCBenchmark.ImplicitEuler()
model = MPCCBenchmark.nosnoc_sliding_mode_ocp_model(N, nfe, collocation; step_eq=:heuristic_mean)

######
# Reformulate problem as a nonlinear program with ComplementOpt
######
ind_cc1, ind_cc2 = ComplementOpt.reformulate_to_vertical!(JuMP.backend(model))
ComplementOpt.reformulate_as_nonlinear_program!(JuMP.backend(model), ComplementOpt.ScholtesRelaxation(1e-6))

######
# Solve problem with Ipopt
######
JuMP.set_optimizer(model, Ipopt.Optimizer)
JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
JuMP.set_optimizer_attribute(model, "bound_push", 1e-1)
JuMP.optimize!(model)

######
# Display solution
######
x = JuMP.value.(model[:x1])[:, 1]
y = JuMP.value.(model[:x2])[:, 1]
vx = JuMP.value.(model[:v1])[:, 1]
vy = JuMP.value.(model[:v2])[:, 1]
u1 = JuMP.value.(model[:u1])
u2 = JuMP.value.(model[:u2])
θ1 = JuMP.value.(model[:θ])[:, 1, :, 1]
θ2 = JuMP.value.(model[:θ])[:, 1, :, 2]
λ1 = JuMP.value.(model[:λ])[:, 1, :, 1].data
λ2 = JuMP.value.(model[:λ])[:, 1, :, 2].data
tm = [0.0; cumsum(JuMP.value.(model[:h])[:])]

fig = plot(
    layout=(3,2),
    sharex=true,
    guidefontsize=7,
    legendfontsize=4,
    titlefontsize=6,
    tickfontsize=4,
    labelfontsize=4,
)
xlb, xub = minimum(x), maximum(x)
ylb, yub = minimum(x), maximum(x)
xvals = xlb:0.01:xub
yvals = ylb:0.01:yub
q1(y) = 0.15*y^2
q2(x) = 0.05*x^3

plot!(q1.(yvals), yvals, label="c1", subplot=1, lw=0.5)
plot!(xvals, q2.(xvals), label="c2", subplot=1, lw=0.5)
plot!(x, y, label="Trajectory", subplot=1, lw=2.0)
title!("Trajectory", subplot=1)

plot!(tm, vx, label="vx", subplot=3)
plot!(tm, vy, label="vy", subplot=3)
title!("Speed", subplot=3)
xlabel!("Time", subplot=3)

plot!(tm[1:nfe:end-1], u1, label="u1", subplot=5)
plot!(tm[1:nfe:end-1], u2, label="u2", subplot=5)
title!("Control", subplot=5)
xlabel!("Time", subplot=5)

plot!(tm[1:end-1], θ1, label="θ1", subplot=2)
plot!(tm[1:end-1], θ2, label="θ2", subplot=2)
title!("Stewart indicator variable θ", subplot=2)
xlabel!("Time", subplot=2)

plot!(tm[1:end-1], λ1, subplot=4, label="λ1")
plot!(tm[1:end-1], λ2, subplot=4, label="λ2")
title!("Stewart multipliers λ", subplot=4)
xlabel!("Time", subplot=4)

plot!(tm[1:end-1], h, label=nothing, subplot=6)
title!("Timestep Δh", subplot=6)
xlabel!("Time", subplot=6)
