using JuMP
using MPCCBenchmark
using ComplementOpt
using Ipopt
using Plots

######
# Import problem from MPCCBenchmark
######
N = 30
nfe = 2
#collocation = MPCCBenchmark.ImplicitEuler()  # 3 stages in MATLAB (n_s = 3)
collocation = MPCCBenchmark.RadauIIA(2)
model = MPCCBenchmark.nosnoc_bilinear_spring_damper_model(N, nfe, collocation; step_eq=:heuristic_mean)


######
# Solve problem with Ipopt
######
JuMP.set_optimizer(model, () -> ComplementOpt.Optimizer(Ipopt.Optimizer()))
MOI.set(model, ComplementOpt.RelaxationMethod(), ComplementOpt.ScholtesRelaxation(1e-3))
JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
JuMP.set_optimizer_attribute(model, "bound_push", 1e-2)
JuMP.set_optimizer_attribute(model, "max_iter", 3000)
JuMP.optimize!(model)

######
# Display solution status
######
println("Termination status: ", JuMP.termination_status(model))
println("Objective value: ", JuMP.objective_value(model))

######
# Extract solution
######
x1 = JuMP.value.(model[:x1])[:, 1]  # position
x2 = JuMP.value.(model[:x2])[:, 1]  # velocity
U = JuMP.value.(model[:U])
θ = JuMP.value.(model[:θ])[:, 1, :]
λ = JuMP.value.(model[:λ]).data[:, 1, :]
g = JuMP.value.(model[:g])[:, 2]
μ = JuMP.value.(model[:μ][:, 0]).data
h = JuMP.value.(model[:h])[:]
Δ = cumsum(h)
tm = [0.0; Δ]

# Reference values for plotting
x_ref = [1.5, 0.0]
u_max = 5.0

######
# Create plots
######
fig = plot(
    layout=(3,1),
    sharex=true,
    guidefontsize=9,
    legendfontsize=7,
    titlefontsize=9,
    tickfontsize=7,
    labelfontsize=8,
    size=(800, 700)
)

# Plot 1: Position
plot!(tm, x1, label="q(t)", subplot=1, linewidth=2, color=:blue)
hline!([x_ref[1]], label="q_ref", subplot=1, linestyle=:dash, color=:black, linewidth=1.5)
xlabel!("t", subplot=1)
ylabel!("q", subplot=1)
title!("Position", subplot=1)

# Plot 2: Velocity
plot!(tm, x2, label="v(t)", subplot=2, linewidth=2, color=:red)
hline!([x_ref[2]], label="v_ref", subplot=2, linestyle=:dash, color=:black, linewidth=1.5)
xlabel!("t", subplot=2)
ylabel!("v", subplot=2)
title!("Velocity", subplot=2)

# Plot 3: Control input
plot!(tm[1:nfe:end-1], U, label="u(t)", subplot=3, linewidth=2, color=:green, seriestype=:steppost)
hline!([u_max, -u_max], label="u_max", subplot=3, linestyle=:dash, color=:black, linewidth=1.5)
xlabel!("t", subplot=3)
ylabel!("u", subplot=3)
title!("Control Input", subplot=3)

display(fig)

# Optional: Additional diagnostic plots
fig2 = plot(
    layout=(2,1),
    sharex=true,
    guidefontsize=9,
    legendfontsize=7,
    titlefontsize=9,
    tickfontsize=7,
    labelfontsize=8,
    size=(800, 500)
)

# Plot Stewart indicator variables
plot!(tm[1:end-1], θ, labels=["θ1 (x<0)"  "θ2 (x>0)"], subplot=1, linewidth=2)
title!("Stewart Indicator Variables θ", subplot=1)
xlabel!("Time", subplot=1)
ylabel!("θ", subplot=1)
#grid!(subplot=1)

# Plot timesteps
plot!(tm[1:end-1], h, label=nothing, subplot=2, linewidth=2)
title!("Timestep Δh", subplot=2)
xlabel!("Time", subplot=2)
ylabel!("h", subplot=2)

display(fig2)

######
# Print final state
######
println("\nFinal state:")
println("  Position: ", x1[end])
println("  Velocity: ", x2[end])
println("\nReference state:")
println("  Position: ", x_ref[1])
println("  Velocity: ", x_ref[2])
