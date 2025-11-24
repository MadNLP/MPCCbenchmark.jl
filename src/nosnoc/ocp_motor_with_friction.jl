# Adapted from NOSBENCH
# https://github.com/nosnoc/nosbench
#
# Nurkanović, A., Pozharskiy, A. & Diehl, M.
# Solving Mathematical Programs with Complementarity Constraints Arising in Nonsmooth Optimal Control.
# Vietnam J. Math. 53, 659–697 (2025).
# https://doi.org/10.1007/s10013-024-00704-z

"""
    nosnoc_motor_with_friction_model

Optimal control of a linear voice-coil motor.

# Reference

Christiansen, Bahne, Helmut Maurer, and Oliver Zirn.
"Optimal control of a voice-coil-motor with coulombic friction."
2008 47th IEEE conference on decision and control. IEEE, 2008.

"""
function nosnoc_motor_with_friction_model(N, nfe, rk::RKScheme; big_M=1e5, rho_h=1e0, step_eq=:heuristic_mean)
    nh = N * nfe      # total number of finite elements
    nf = 2            # total number of nonsmooth modes
    nc = length(rk.c) # total number of intermediate integration points

    T = 0.08
    hmin, hmax = 0.5*T/nh, 10.0*T/nh

    # Parameters
    m1 = 1.03 # slide mass
    m2 = 0.56 # load mass
    k = 2.4e3 #  spring constant N/m
    c = 0.00 # damping
    U_max = 5.0 # voltage Back-EMF, U = K_s*v_1;1
    R  = 2.0 # coil resistance ohm
    L = 2e-3 # inductivity, henry
    K_F = 12 # force constant N/A ; F_L = K_F*I; % Lorenz force
    K_S = 12 # Vs/m (not provided in the paper above)
    F_R = 2. # guide friction force, N

    model = Model()
    JuMP.set_name(model, "NOSNOC-MWFOCP")
    # State
    @variable(model, x1[i=1:(nh+1), j=1:(nc+1)])
    @variable(model, x2[i=1:(nh+1), j=1:(nc+1)])
    @variable(model, v1[i=1:(nh+1), j=1:(nc+1)])
    @variable(model, v2[i=1:(nh+1), j=1:(nc+1)])
    @variable(model, I[i=1:(nh+1), j=1:(nc+1)])
    @variable(model, -U_max <= U[1:N] <= U_max)
    # Step in numerical time
    @variable(model, hmin <= h[1:nfe, 1:N] <= hmax, start=5*T/nh)
    # Stewart's multipliers
    @variable(model, 0.0 <= θ[1:nh, 1:nc, 1:nf], start=0.5)
    @variable(model, 0.0 <= λ[1:nh, 0:nc, 1:nf], start=1.0)
    @variable(model, μ[1:nh, 0:nc])
    # Switch detection
    @variable(model, 0.0 <= λs[1:nh, 1:nf], start=1.0)
    @variable(model, 0.0 <= θs[1:nh, 1:nf], start=0.5)
    if step_eq == :lcc
        @variable(model, 0.0 <= πλ[1:(nh-1), 1:nf])
        @variable(model, 0.0 <= πθ[1:(nh-1), 1:nf])
        @variable(model, 0.0 <= τ1[1:(nh-1), 1:nf])
        @variable(model, 0.0 <= τ2[1:(nh-1), 1:nf])
        @variable(model, 0.0 <= vcc[1:(nh-1), 1:nf])
    end

    @expression(
        model,
        cost,
        1/N*sum(U[i]^2 for i in 1:N)
    )

    # Initial and final positions
    @constraints(model, begin
        x1[1, 1] == 0.0
        x2[1, 1] == 0.0
        v1[1, 1] == 0.0
        v2[1, 1] == 0.0
        I[1, 1] == 0.0
        x1[nh+1, 1] == 0.01
        x2[nh+1, 1] == 0.01
        v1[nh+1, 1] == 0.0
        v2[nh+1, 1] == 0.0
        I[nh+1, 1] == 0.0
    end)
    # Stewart position function
    @expression(
        model,
        g[i=1:nh, j=1:(nc+1)],
        v1[i, j]
    )
    # Initial position for multipliers
    @constraint(model, λ[1, 0, 1] ==  g[1, 1] - μ[1, 0])
    @constraint(model, λ[1, 0, 2] == -g[1, 1] - μ[1, 0])

    # Dynamics
    @expressions(
        model,
        begin
            dx1[i=1:nh, j=1:nc], v1[i, j]
            dx2[i=1:nh, j=1:nc], v2[i, j]
            dv1[i=1:nh, j=1:nc], -k/m1*x1[i, j] - c/m1*v1[i, j] + k/m1*x2[i, j] + c/m1*v2[i, j] + K_F/m1*I[i, j] - F_R/m1 * (θ[i, j, 1] - θ[i, j, 2])
            dv2[i=1:nh, j=1:nc], k/m2*x1[i, j] + c/m2*v1[i, j] - k/m2*x2[i, j] - c/m2*v2[i, j]
            dI[i=1:nh, j=1:nc], -K_S/L*v1[i, j] - R/L*I[i, j] + 1/L * U[div(i-1, nfe)+1]
        end
    )
    # Collocations
    @constraints(
        model,
        begin
            con_dx1[i=1:nh, j=2:(nc+1)],
            x1[i, j] == x1[i, 1] + h[i] * sum(rk.a[j-1, k] * dx1[i, k] for k in 1:nc)
            con_dx2[i=1:nh, j=2:(nc+1)],
            x2[i, j] == x2[i, 1] + h[i] * sum(rk.a[j-1, k] * dx2[i, k] for k in 1:nc)
            con_dv1[i=1:nh, j=2:(nc+1)],
            v1[i, j] == v1[i, 1] + h[i] * sum(rk.a[j-1, k] * dv1[i, k] for k in 1:nc)
            con_dv2[i=1:nh, j=2:(nc+1)],
            v2[i, j] == v2[i, 1] + h[i] * sum(rk.a[j-1, k] * dv2[i, k] for k in 1:nc)
            con_dI[i=1:nh, j=2:(nc+1)],
            I[i, j] == I[i, 1] + h[i] * sum(rk.a[j-1, k] * dI[i, k] for k in 1:nc)
        end
    )
    # Continuity
    # w.r.t primal variables
    @constraints(
        model,
        begin
            [i=1:nh], x1[i+1, 1] == x1[i, 1] + h[i] * sum(rk.b[k] * dx1[i, k] for k in 1:nc)
            [i=1:nh], x2[i+1, 1] == x2[i, 1] + h[i] * sum(rk.b[k] * dx2[i, k] for k in 1:nc)
            [i=1:nh], v1[i+1, 1] == v1[i, 1] + h[i] * sum(rk.b[k] * dv1[i, k] for k in 1:nc)
            [i=1:nh], v2[i+1, 1] == v2[i, 1] + h[i] * sum(rk.b[k] * dv2[i, k] for k in 1:nc)
            [i=1:nh], I[i+1, 1] == I[i, 1] + h[i] * sum(rk.b[k] * dI[i, k] for k in 1:nc)
        end
    )
    # w.r.t dual variables
    @constraints(model, begin
        [i=1:(nh-1), j=1:nf], λ[i, nc, j] == λ[i+1, 0, j]
    end)
    # Stewart model
    @constraints(model, begin
        [i=1:nh, j=1:nc], sum(θ[i, j, k] for k in 1:nf) == 1.0
        [i=1:nh, j=0:nc],  g[i, j+1] - λ[i, j, 1] - μ[i, j] == 0
        [i=1:nh, j=0:nc], -g[i, j+1] - λ[i, j, 2] - μ[i, j] == 0
    end)
    # Switch detection
    @constraints(
        model,
        begin
            [t=1:nh, j=1:nf], sum(θ[t, i, j] for i in 1:nc) == θs[t, j]
            [t=1:nh, j=1:nf], sum(λ[t, i, j] for i in 0:nc) == λs[t, j]
            [t=1:nh, i=1:nc, j=1:nf], [θ[t, i, j], λs[t, j]] ∈ MOI.Complements(2)
        end
    )

    if step_eq == :lcc
        # Logical constraints
        @constraints(
            model,
            begin
                # OR
                [t=1:(nh-1), j=1:nf], πλ[t, j] >= λs[t, j]
                [t=1:(nh-1), j=1:nf], πλ[t, j] >= λs[t+1, j]
                [t=1:(nh-1), j=1:nf], πλ[t, j] <= λs[t, j] + λs[t+1, j]
                # OR
                [t=1:(nh-1), j=1:nf], πθ[t, j] >= θs[t, j]
                [t=1:(nh-1), j=1:nf], πθ[t, j] >= θs[t+1, j]
                [t=1:(nh-1), j=1:nf], πθ[t, j] <= θs[t, j] + θs[t+1, j]
                # AND
                [t=1:(nh-1), j=1:nf], vcc[t, j] <= πλ[t, j]
                [t=1:(nh-1), j=1:nf], vcc[t, j] <= πθ[t, j]
                [t=1:(nh-1), j=1:nf], vcc[t, j] >= πθ[t, j] - τ1[t, j]
                [t=1:(nh-1), j=1:nf], πθ[t, j] - πλ[t, j] == τ1[t, j] - τ2[t, j]
                [t=1:(nh-1), j=1:nf], [τ1[t, j], τ2[t, j]] ∈ MOI.Complements(2)
            end
        )
        # Step equilibration
        @constraints(
            model,
            begin
                [t=1:(nh-1)], -big_M * sum(vcc[t, j] for j in 1:nf) <= (h[t+1] - h[t])
                [t=1:(nh-1)], (h[t+1] - h[t]) <= big_M * sum(vcc[t, j] for j in 1:nf)
            end
        )
        @expression(model, step_eq_cost, 0)
    elseif step_eq == :heuristic_mean
        @expression(model, step_eq_cost, rho_h * sum(1*(h .- (T/nh)) .^ 2))
    end
    @constraint(model, [j=1:N], sum(h[i, j] for i in 1:nfe) == T / N)

    @objective(model, Min, cost + step_eq_cost)
    return model
end

