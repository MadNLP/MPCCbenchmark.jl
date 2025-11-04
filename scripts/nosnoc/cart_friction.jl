
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
model = MPCCBenchmark.nosnoc_cart_pole_with_friction_model(N, nfe, collocation; step_eq=:heuristic_mean)

######
# Reformulate problem as a nonlinear program with ComplementOpt
######
ind_cc1, ind_cc2 = ComplementOpt.reformulate_to_vertical!(JuMP.backend(model))
ComplementOpt.reformulate_as_nonlinear_program!(JuMP.backend(model), ComplementOpt.ScholtesRelaxation(1e-8))

######
# Solve problem with Ipopt
#####
JuMP.set_optimizer(model, Ipopt.Optimizer)
JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
JuMP.set_optimizer_attribute(model, "bound_push", 1e-1)
JuMP.optimize!(model)

######
# Display solution
######
qx = JuMP.value.(model[:q])[:, 1, 1]
qy = JuMP.value.(model[:q])[:, 1, 2]
vx = JuMP.value.(model[:v])[:, 1, 1]
vy = JuMP.value.(model[:v])[:, 1, 2]
θ = JuMP.value.(model[:θ])[:, 1, :]
u = JuMP.value.(model[:u])
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
plot!(qx, qy, subplot=1, label=nothing)
xlabel!("qx", subplot=1)
ylabel!("qy", subplot=1)
title!("Trajectory", subplot=1)

plot!(tm, vx, label="vx", subplot=3)
plot!(tm, vy, label="vy", subplot=3)
title!("Speed", subplot=3)
xlabel!("Time", subplot=3)

plot!(tm[1:nfe:end-1], u, label="u", subplot=5)
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
