
using JuMP
using MPCCBenchmark
using ComplementOpt
using Ipopt
using Plots


######
# Import problem from MPCCBenchmark
######
N = 100
nfe = 2
collocation = MPCCBenchmark.RadauIIA(1)
model = MPCCBenchmark.nosnoc_schumacher_model(N, nfe, collocation; step_eq=:lcc, big_M=1e-2)

######
# Solve problem with ComplementOpt+Ipopt
#####

JuMP.set_optimizer(model, () -> ComplementOpt.Optimizer(Ipopt.Optimizer()))
MOI.set(model, ComplementOpt.RelaxationMethod(), ComplementOpt.ScholtesRelaxation(1e-5))
JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
JuMP.set_optimizer_attribute(model, "bound_push", 1e-1)
JuMP.optimize!(model)

######
# Display solution
######
qx = JuMP.value.(model[:qx])[:, 1]
qy = JuMP.value.(model[:qy])[:, 1]
vx = JuMP.value.(model[:vx])[:, 1]
vy = JuMP.value.(model[:vy])[:, 1]
θ = JuMP.value.(model[:θ])[:, 1, :]
λ = JuMP.value.(model[:λ]).data[:, 1, :]
λs = JuMP.value.(model[:λs])[:, :]
g = JuMP.value.(model[:g])[:, 2]
μ = JuMP.value.(model[:μ][:, 0]).data
h = JuMP.value.(model[:h])[:]
Δ = cumsum(h)
τ = JuMP.value(model[:sot])
tm = [0.0; Δ .* τ]
u1 = JuMP.value.(model[:a])
u2 = JuMP.value.(model[:s])

fig = plot(
    layout=(3,2),
    sharex=true,
    guidefontsize=7,
    legendfontsize=4,
    titlefontsize=6,
    tickfontsize=4,
    labelfontsize=4,
)
plot!(qx, qy, subplot=1)
xlabel!("qx", subplot=1)
ylabel!("qy", subplot=1)
xx = 0:0.1:3π
plot!(xx, MPCCBenchmark._schumacher_track.(xx) .+ 0.25, label="Track lb", subplot=1)
plot!(xx, MPCCBenchmark._schumacher_track.(xx) .- 0.25, label="Track ub", subplot=1)
title!("Trajectory", subplot=1)

plot!(tm, vx, label="vx", subplot=3)
plot!(tm, vy, label="vy", subplot=3)
title!("Speed", subplot=3)
xlabel!("Time", subplot=3)

plot!(tm[1:nfe:end-1], u1, label="Acceleration", subplot=5)
plot!(tm[1:nfe:end-1], u2, label="Wheel", subplot=5)
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

