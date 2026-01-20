
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
model = MPCCBenchmark.nosnoc_motor_with_friction_model(N, nfe, collocation; step_eq=:heuristic_mean)

######
# Solve problem with Ipopt
######
JuMP.set_optimizer(model, () -> ComplementOpt.Optimizer(Ipopt.Optimizer()))
MOI.set(model, ComplementOpt.RelaxationMethod(), ComplementOpt.ScholtesRelaxation(1e-6))
JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
JuMP.set_optimizer_attribute(model, "bound_push", 1e-2)
JuMP.optimize!(model)

######
# Display solution
######
x = JuMP.value.(model[:x1])[:, 1]
y = JuMP.value.(model[:x2])[:, 1]
vx = JuMP.value.(model[:v1])[:, 1]
vy = JuMP.value.(model[:v2])[:, 1]
I = JuMP.value.(model[:I])[:, 1]
U = JuMP.value.(model[:U])
θ = JuMP.value.(model[:θ])[:, 1, :]
λ = JuMP.value.(model[:λ]).data[:, 1, :]
g = JuMP.value.(model[:g])[:, 2]
μ = JuMP.value.(model[:μ][:, 0]).data
h = JuMP.value.(model[:h])[:]
Δ = cumsum(h)
tm = [0.0; Δ]

fig = plot(
    layout=(3,2),
    sharex=true,
    guidefontsize=7,
    legendfontsize=4,
    titlefontsize=6,
    tickfontsize=4,
    labelfontsize=4,
)
plot!(tm, x, label="x", subplot=1)
plot!(tm, y, label="y", subplot=1)
title!("Trajectory", subplot=1)
xlabel!("Time", subplot=1)

plot!(tm, vx, label="vx", subplot=3)
plot!(tm, vy, label="vy", subplot=3)
title!("Speed", subplot=3)
xlabel!("Time", subplot=3)

plot!(tm[1:nfe:end-1], U, label="U", subplot=5)
title!("Control", subplot=5)
xlabel!("Time", subplot=5)

plot!(tm[1:end-1], θ, labels=["θ1"  "θ2"], subplot=2)
title!("Stewart indicator variable θ", subplot=2)
xlabel!("Time", subplot=2)

plot!(tm[1:end-1], JuMP.value.(model[:λs])[:, :], subplot=4, labels=["λ1"  "λ2"])
title!("Stewart multipliers λ", subplot=4)
xlabel!("Time", subplot=4)

plot!(tm[1:end-1], h, label=nothing, subplot=6)
title!("Timestep Δh", subplot=6)
xlabel!("Time", subplot=6)
