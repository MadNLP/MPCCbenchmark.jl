using JuMP
using MPCCBenchmark
using ComplementOpt
using Ipopt
using Plots

######
# Import problem from MPCCBenchmark
######
dt = 0.05
T = 5.0
N = Int(T/dt)  # Number of control intervals
nfe = 3        # Number of finite elements per control interval
collocation = MPCCBenchmark.ImplicitEuler()  # RADAU_IIA with n_s=2 in MATLAB
model = MPCCBenchmark.nosnoc_acrobot_model(N, nfe, collocation; step_eq=:heuristic_mean)

######
# Reformulate problem as a nonlinear program with ComplementOpt
######
ind_cc1, ind_cc2 = ComplementOpt.reformulate_to_vertical!(JuMP.backend(model))
ComplementOpt.reformulate_as_nonlinear_program!(JuMP.backend(model), ComplementOpt.ScholtesRelaxation(1e-2))

######
# Solve problem with Ipopt
######
JuMP.set_optimizer(model, Ipopt.Optimizer)
JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
JuMP.set_optimizer_attribute(model, "bound_push", 1e-2)
JuMP.set_optimizer_attribute(model, "max_iter", 5000)
JuMP.set_optimizer_attribute(model, "tol", 1e-6)
JuMP.set_optimizer_attribute(model, "print_level", 5)
JuMP.optimize!(model)

######
# Display solution status
######
println("Termination status: ", JuMP.termination_status(model))
println("Objective value: ", JuMP.objective_value(model))

######
# Extract solution
######
q1 = JuMP.value.(model[:q1])[:, 1]      # joint 1 angle
q2 = JuMP.value.(model[:q2])[:, 1]      # joint 2 angle
q1_dot = JuMP.value.(model[:q1_dot])[:, 1]  # joint 1 velocity
q2_dot = JuMP.value.(model[:q2_dot])[:, 1]  # joint 2 velocity
U = JuMP.value.(model[:U])              # control inputs [tau1, tau2]
θ1 = JuMP.value.(model[:θ1])[:, 1, :]  # Stewart indicators for subsystem 1
θ2 = JuMP.value.(model[:θ2])[:, 1, :]  # Stewart indicators for subsystem 2
h = JuMP.value.(model[:h])[:]
Δ = cumsum(h)
tm = [0.0; Δ]

# Reference values for plotting
x_ref = [π, 0.0, 0.0, 0.0]

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

# Plot 1: Joint angles
plot!(tm, q1, label="q₁(t)", subplot=1, linewidth=2, color=:blue)
plot!(tm, q2, label="q₂(t)", subplot=1, linewidth=2, color=:red)
hline!([x_ref[1]], label="q₁_ref", subplot=1, linestyle=:dash, color=:blue, linewidth=1.5, alpha=0.6)
hline!([x_ref[2]], label="q₂_ref", subplot=1, linestyle=:dash, color=:red, linewidth=1.5, alpha=0.6)
xlabel!("t", subplot=1)
ylabel!("q", subplot=1)
title!("Joint Angles", subplot=1)

# Plot 2: Joint velocities
plot!(tm, q1_dot, label="q₁_dot(t)", subplot=2, linewidth=2, color=:blue)
plot!(tm, q2_dot, label="q₂_dot(t)", subplot=2, linewidth=2, color=:red)
hline!([x_ref[3]], label="q₁_dot_ref", subplot=2, linestyle=:dash, color=:blue, linewidth=1.5, alpha=0.6)
hline!([x_ref[4]], label="q₂_dot_ref", subplot=2, linestyle=:dash, color=:red, linewidth=1.5, alpha=0.6)
xlabel!("t", subplot=2)
ylabel!("v", subplot=2)
title!("Joint Velocities", subplot=2)

# Plot 3: Control torques
plot!(tm[1:nfe:end-1], U[:, 1], label="τ₁(t)", subplot=3, linewidth=2, color=:blue, seriestype=:steppost)
plot!(tm[1:nfe:end-1], U[:, 2], label="τ₂(t)", subplot=3, linewidth=2, color=:red, seriestype=:steppost)
xlabel!("t", subplot=3)
ylabel!("u", subplot=3)
title!("Control Torques", subplot=3)

display(fig)

# Optional: Additional diagnostic plots
fig2 = plot(
    layout=(3,1),
    sharex=true,
    guidefontsize=9,
    legendfontsize=7,
    titlefontsize=9,
    tickfontsize=7,
    labelfontsize=8,
    size=(800, 700)
)

# Plot Stewart indicator variables for joint 1
plot!(tm[1:end-1], θ1, labels=["θ₁₁ (q₁_dot>0)"  "θ₁₂ (q₁_dot<0)"], subplot=1, linewidth=2)
title!("Stewart Indicators θ₁ (Joint 1 Friction)", subplot=1)
xlabel!("Time", subplot=1)
ylabel!("θ₁", subplot=1)

# Plot Stewart indicator variables for joint 2
plot!(tm[1:end-1], θ2, labels=["θ₂₁ (q₂_dot>0)"  "θ₂₂ (q₂_dot<0)"], subplot=2, linewidth=2)
title!("Stewart Indicators θ₂ (Joint 2 Friction)", subplot=2)
xlabel!("Time", subplot=2)
ylabel!("θ₂", subplot=2)

# Plot timesteps
plot!(tm[1:end-1], h, label=nothing, subplot=3, linewidth=2)
title!("Timestep Δh", subplot=3)
xlabel!("Time", subplot=3)
ylabel!("h", subplot=3)

display(fig2)

######
# Print final state
######
println("\nFinal state:")
println("  Joint 1 angle: ", q1[end])
println("  Joint 2 angle: ", q2[end])
println("  Joint 1 velocity: ", q1_dot[end])
println("  Joint 2 velocity: ", q2_dot[end])
println("\nReference state:")
println("  Joint 1 angle: ", x_ref[1], " (", rad2deg(x_ref[1]), " deg)")
println("  Joint 2 angle: ", x_ref[2], " (", rad2deg(x_ref[2]), " deg)")
println("  Joint 1 velocity: ", x_ref[3])
println("  Joint 2 velocity: ", x_ref[4])