using JuMP
using MPCCBenchmark
using ComplementOpt
using Ipopt
using Plots

######
# Import problem from MPCCBenchmark
######
N   = 100
nfe = 3
collocation = MPCCBenchmark.ImplicitEuler()

model = MPCCBenchmark.nosnoc_liquid_gas_tank_model(
    N,
    nfe,
    collocation;
    step_eq = :heuristic_mean,
)

######
# Reformulate problem as a nonlinear program with ComplementOpt
######
ind_cc1, ind_cc2 = ComplementOpt.reformulate_to_vertical!(JuMP.backend(model))
ComplementOpt.reformulate_as_nonlinear_program!(
    JuMP.backend(model),
    ComplementOpt.ScholtesRelaxation(1e0),
)

######
# Solve with Ipopt
######
JuMP.set_optimizer(model, Ipopt.Optimizer)
JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
JuMP.set_optimizer_attribute(model, "bound_push", 1e-1)
JuMP.optimize!(model)

######
# Extract solution
######
xval = JuMP.value.(model[:x])           # (nh+1) × (nc+1) × 2
u    = JuMP.value.(model[:u])
h    = JuMP.value.(model[:h])[:]

Δ  = cumsum(h)
tm = [0.0; Δ]                           # length nh+1

MG = xval[:, 1, 1]                      # M_G at FE boundaries
ML = xval[:, 1, 2]                      # M_L at FE boundaries

# Parameters (for post-processing)
V     = 10.0
V_s   = 5.0
T_gas = 300.0
ρ_L   = 50.0
Rgas  = 0.0821

P = MG .* Rgas .* T_gas ./ (V .- ML ./ ρ_L)
c = ML ./ ρ_L .- V_s

######
# Plots
######
plt = plot(layout = (4, 1), size = (800, 900))

# M_G
plot!(
    plt[1],
    tm[1:end-1],
    MG[1:end-1],
    label = "M_G",
)
xlabel!(plt[1], "Time")
ylabel!(plt[1], "M_G")

# M_L
plot!(
    plt[2],
    tm[1:end-1],
    ML[1:end-1],
    label = "M_L",
)
xlabel!(plt[2], "Time")
ylabel!(plt[2], "M_L")

# Pressure P
plot!(
    plt[3],
    tm[1:end-1],
    P[1:end-1],
    label = "P",
)
xlabel!(plt[3], "Time")
ylabel!(plt[3], "P")

# Control u (stair / piecewise-constant per stage)
plot!(
    plt[4],
    tm[1:nfe:end-1],
    u,
    label    = "u",
    linetype = :steppost,
)
xlabel!(plt[4], "Time")
ylabel!(plt[4], "u")
ylims!(plt[4], (0.07, 0.12))

display(plt)

# c(x(t)) = M_L/ρ_L - V_s
plt_c = plot(
    tm[1:end-1],
    c[1:end-1],
    label = "c(x(t))",
)
xlabel!("Time")
ylabel!("c(x(t))")
display(plt_c)
