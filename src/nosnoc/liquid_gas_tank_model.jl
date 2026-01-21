# Ideal gas–liquid tank with pressure control valve
#
# Based on:
# Saif R. Kazi, Mandar Thombre, Lorenz Biegler,
# "Globally Convergent Method for Optimal Control of Hybrid Dynamical Systems"
# 12th IFAC Symposium on Advanced Control of Chemical Processes, 2024.

using JuMP
#const MOI = JuMP.MathOptInterface

"""
    nosnoc_liquid_gas_tank_model(N, nfe, rk::RKScheme; big_M=1e5, step_eq=:lcc)

Hybrid OCP for the ideal gas–liquid tank with a pressure control valve.

Two modes (PSS/Stewart):
- Mode 1 (c > 0): liquid mode, dynamics `f₁ = [F_G; F_L - L]`
- Mode 2 (c < 0): gas mode,    dynamics `f₂ = [F_G - G; F_L]`

Switching function:
    c(x) = M_L/ρ_L - V_s
"""
function nosnoc_liquid_gas_tank_model(N, nfe, rk::RKScheme; big_M=1e5, step_eq=:heuristic_mean, rho_h=1e0)
    nh = N * nfe      # total number of finite elements
    nf = 2            # number of nonsmooth modes
    nc = length(rk.c) # collocation points per FE

    # Horizon and timestep bounds (MATLAB: N_stages = 100, DT = 0.25)
    T = 25.0
    hmin, hmax = 0.0 / nh, 2.0 * T / nh

    # Parameters
    F_L   = 2.5
    F_G   = 0.1
    V     = 10.0
    V_s   = 5.0
    T_gas = 300.0
    P_out = 1.0
    ρ_L   = 50.0
    k_L   = 1.0
    k_G   = 1.0
    Rgas  = 0.0821

    # Initial state (M_G, M_L)
    M_G0 = 6.83
    M_L0 = 260.0
    x0   = [M_G0, M_L0]

    u_max = 0.15
    u_ref = 0.10

    # Simple linear guess of state between initial and some nominal value (here just constant)
    x_guess = [x0 for _ in 1:(nh+1)]

    model = Model()
    JuMP.set_name(model, "NOSNOC-LIQUID-GAS-TANK")

    # State: x = [M_G, M_L]
    @variable(model, x[i=1:(nh+1), j=1:(nc+1), k=1:2], start = x_guess[i][k])
    # Control
    @variable(model, 0.0 <= u[1:N] <= u_max, start = u_ref)
    # Step sizes (single vector of length nh)
    @variable(model, hmin <= h[1:nh] <= hmax, start = T / nh)

    # Stewart multipliers
    @variable(model, 0.0 <= θ[1:nh, 1:nc, 1:nf], start = 0.5)
    @variable(model, 0.0 <= λ[1:nh, 0:nc, 1:nf], start = 1.0)
    @variable(model, μ[1:nh, 0:nc])

    # Switch detection
    @variable(model, 0.0 <= λs[1:nh, 1:nf], start = 0.5)
    @variable(model, 0.0 <= θs[1:nh, 1:nf], start = 1.0)
    if step_eq == :lcc
        @variable(model, 0.0 <= πλ[1:(nh-1), 1:nf])
        @variable(model, 0.0 <= πθ[1:(nh-1), 1:nf])
        @variable(model, 0.0 <= τ1[1:(nh-1), 1:nf])
        @variable(model, 0.0 <= τ2[1:(nh-1), 1:nf])
        @variable(model, 0.0 <= vcc[1:(nh-1), 1:nf])
    end

    # Stage sampling of state at FE boundaries
    @expression(
        model,
        x_stage[i=1:(N+1)],
        [
            x[(i-1)*nfe+1, 1, 1],
            x[(i-1)*nfe+1, 1, 2]
        ]
    )

    # Cost: 100*(u - u_ref)^2 + terminal M_L
    @expression(
        model,
        terminal_cost,
        (1 / N) * sum(100.0 * (u[i] - u_ref)^2 for i in 1:N) +
        x_stage[N+1][2]
    )

    # Initial condition
    @constraints(model, begin
        x_init1, x[1, 1, 1] == x0[1]
        x_init2, x[1, 1, 2] == x0[2]
    end)

    # Stewart switching function: c(x) = M_L/ρ_L - V_s
    @expression(
        model,
        g[i=1:nh, j=1:(nc+1)],
        x[i, j, 2] / ρ_L - V_s
    )

    # Vector fields f = θ₁ f₁ + θ₂ f₂
    # f₁ = [F_G; F_L - L],  f₂ = [F_G - G; F_L]
    # L = k_L * u * (P - P_out),  G = k_G * u * (P - P_out)
    # P = M_G * R * T / (V - M_L/ρ_L)
    @expression(
        model,
        f[i=1:nh, j=1:nc],
        begin
            stage = div(i - 1, nfe) + 1
            P_ij = x[i, j, 1] * Rgas * T_gas / (V - x[i, j, 2] / ρ_L)
            ΔP   = P_ij - P_out
            L_ij = k_L * u[stage] * ΔP
            G_ij = k_G * u[stage] * ΔP
            [
                F_G - θ[i, j, 2] * G_ij ;
                F_L - θ[i, j, 1] * L_ij
            ]
        end
    )

    @expressions(
        model,
        begin
            dx[i=1:nh, j=1:nc, k=1:2], f[i, j][k]
        end
    )

    # Collocation equations
    @constraints(
        model,
        begin
            # within each FE
            con_dx[i=1:nh, j=2:(nc+1), k=1:2],
            x[i, j, k] == x[i, 1, k] + h[i] * sum(rk.a[j-1, ℓ] * dx[i, ℓ, k] for ℓ in 1:nc)
        end
    )

    # Continuity between elements
    @constraints(
        model,
        begin
            [i=1:nh, k=1:2],
            x[i+1, 1, k] == x[i, 1, k] + h[i] * sum(rk.b[ℓ] * dx[i, ℓ, k] for ℓ in 1:nc)
        end
    )

    # Stewart constraints
    @constraints(model, begin
        [i=1:nh, j=1:nc], sum(θ[i, j, k] for k in 1:nf) == 1.0
        [i=1:nh, j=0:nc],  g[i, j+1] - λ[i, j, 1] - μ[i, j] == 0
        [i=1:nh, j=0:nc], -g[i, j+1] - λ[i, j, 2] - μ[i, j] == 0
    end)

    # Switch detection + complementary slackness
    @constraints(
        model,
        begin
            [t=1:nh, j=1:nf], sum(θ[t, i, j] for i in 1:nc) == θs[t, j]
            [t=1:nh, j=1:nf], sum(λ[t, i, j] for i in 0:nc) == λs[t, j]
            [t=1:nh, i=1:nc, j=1:nf], [θ[t, i, j], λs[t, j]] ∈ MOI.Complements(2)
        end
    )

    # Optional logical step-equilibration constraints
    if step_eq == :lcc
        @constraints(
            model,
            begin
                # logical constraints
                [t=1:(nh-1), j=1:nf], πλ[t, j] + λs[t, j] == λs[t+1, j]
                [t=1:(nh-1), j=1:nf], πθ[t, j] + θs[t+1, j] == θs[t, j]
                [t=1:(nh-1), j=1:nf], θs[t, j] - πλ[t, j] - τ1[t, j] == 0
                [t=1:(nh-1), j=1:nf], λs[t, j] - πθ[t, j] - τ2[t, j] == 0
                [t=1:(nh-1), j=1:nf], [τ1[t, j], τ2[t, j]] ∈ MOI.Complements(2)
            end
        )
        # step equilibration via big-M
        @constraints(
            model,
            begin
                [t=1:(nh-1)], -big_M * sum(vcc[t, j] for j in 1:nf) <= (h[t+1] - h[t])
                [t=1:(nh-1)],  (h[t+1] - h[t]) <= big_M * sum(vcc[t, j] for j in 1:nf)
            end
        )
        @expression(model, step_eq_cost, 0.0)
    elseif step_eq == :heuristic_mean
        @expression(model, step_eq_cost, rho_h * sum((h .- (T / nh)) .^ 2))
    else
        @expression(model, step_eq_cost, 0.0)
    end

    # Per-stage time consistency: sum over FE steps per stage
    @constraint(
        model,
        [j=1:N],
        sum(h[(j-1)*nfe + 1 : j*nfe]) == T / N
    )

    @objective(model, Min, terminal_cost + step_eq_cost)
    return model
end
