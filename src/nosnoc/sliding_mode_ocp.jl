# Adapted from NOSBENCH
# https://github.com/nosnoc/nosbench
#
# Nurkanović, A., Pozharskiy, A. & Diehl, M.
# Solving Mathematical Programs with Complementarity Constraints Arising in Nonsmooth Optimal Control.
# Vietnam J. Math. 53, 659–697 (2025).
# https://doi.org/10.1007/s10013-024-00704-z


"""
    nosnoc_sliding_mode_ocp_model

Optimal control problem with multiple sliding modes.

# Reference

Nurkanović, Armin, et al.
"Finite elements with switch detection for direct optimal control of nonsmooth systems."
Numerische Mathematik 156.3 (2024): 1115-1162.

"""
function nosnoc_sliding_mode_ocp_model(N, nfe, rk::RKScheme; big_M=1e5, step_eq=:lcc)
    nh = N * nfe      # total number of finite elements
    nf = 2            # total number of nonsmooth modes
    nc = length(rk.c) # total number of intermediate integration points
    nsys = 2          # total number of nonsmooth systems

    T = 4.0

    hmin, hmax = 0.1/nh, 10.0*T/nh

    p = 2; a = 0.15; a1 = 0
    b = -0.05; q = 3

    # Parameters
    u_max = 10.0
    x0 = [2*pi/3, pi/3, 0.0, 0.0];
    x_target = [-pi/6, -pi/4]

    x1_guess = [x0[1] + i/(nh+1) * (x_target[1] - x0[1]) for i in 1:nh+1]
    x2_guess = [x0[2] + i/(nh+1) * (x_target[2] - x0[2]) for i in 1:nh+1]

    model = Model()
    JuMP.set_name(model, "NOSNOC-Motor-Friction")
    # State
    @variable(model, x1[i=1:(nh+1), j=1:(nc+1)], start=x1_guess[i])
    @variable(model, x2[i=1:(nh+1), j=1:(nc+1)], start=x2_guess[i])
    @variable(model, v1[i=1:(nh+1), j=1:(nc+1)])
    @variable(model, v2[i=1:(nh+1), j=1:(nc+1)])
    @variable(model, -u_max <= u1[1:N] <= u_max)
    @variable(model, -u_max <= u2[1:N] <= u_max)
    # Step in numerical time
    @variable(model, hmin <= h[1:nfe, 1:N] <= hmax, start=T/nh)
    # Stewart's multipliers
    @variable(model, 0.0 <= θ[1:nh, 1:nc, 1:nf, 1:nsys], start=0.5)
    @variable(model, 0.0 <= λ[1:nh, 0:nc, 1:nf, 1:nsys], start=0.0)
    @variable(model, μ[1:nh, 0:nc, 1:nsys])
    # Switch detection
    @variable(model, 0.0 <= λs[1:nh, 1:nf, 1:nsys], start=0.0)
    @variable(model, 0.0 <= θs[1:nh, 1:nf, 1:nsys], start=0.5)
    # TODO
    if step_eq == :lcc
        @variable(model, 0.0 <= πλ[1:(nh-1), 1:nf])
        @variable(model, 0.0 <= πθ[1:(nh-1), 1:nf])
        @variable(model, 0.0 <= τ1[1:(nh-1), 1:nf])
        @variable(model, 0.0 <= τ2[1:(nh-1), 1:nf])
        @variable(model, 0.0 <= vcc[1:(nh-1), 1:nf])
    end

    # TODO
    @expression(
        model,
        terminal_cost,
        1/N*sum(v1[(i-1)*nfe+1, 1]^2 + v2[(i-1)*nfe+1, 1]^2 for i in 1:N)
    )

    # Initial and final positions
    @constraints(model, begin
        x1[1, 1] == x0[1]
        x2[1, 1] == x0[2]
        v1[1, 1] == x0[3]
        v2[1, 1] == x0[4]
        x1[nh+1, 1] == x_target[1]
        x2[nh+1, 1] == x_target[2]
    end)
    # Stewart position function
    ## system 1
    @expression(
        model,
        g1[i=1:nh, j=1:(nc+1)],
        x1[i, j] + a*(x2[i, j] - a1)^p
    )
    ## system 2
    @expression(
        model,
        g2[i=1:nh, j=1:(nc+1)],
        x2[i, j] + b*x1[i, j]^q
    )
    # Initial position for multipliers
    for (s, g) in enumerate([g1, g2])
        @constraint(model, λ[1, 0, 1, s] ==  g[1, 1] - μ[1, 0, s])
        @constraint(model, λ[1, 0, 2, s] == -g[1, 1] - μ[1, 0, s])
    end

    # Dynamics
    @expressions(
        model,
        begin
            dx1[i=1:nh, j=1:nc], v1[i, j] + (θ[i, j, 2, 1] - θ[i, j, 1, 1])
            dx2[i=1:nh, j=1:nc], v2[i, j] + (θ[i, j, 2, 2] - θ[i, j, 1, 2])
            dv1[i=1:nh, j=1:nc], u1[div(i-1, nfe) + 1]
            dv2[i=1:nh, j=1:nc], u2[div(i-1, nfe) + 1]
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
        end
    )
    # w.r.t dual variables
    @constraints(model, begin
        [i=1:(nh-1), j=1:nf, s=1:nsys], λ[i, nc, j, s] == λ[i+1, 0, j, s]
    end)
    # Stewart model
    for (s, g) in enumerate([g1, g2])
        @constraints(model, begin
            [i=1:nh, j=1:nc], sum(θ[i, j, k, s] for k in 1:nf) == 1.0
            [i=1:nh, j=0:nc],  g[i, j+1] - λ[i, j, 1, s] - μ[i, j, s] == 0
            [i=1:nh, j=0:nc], -g[i, j+1] - λ[i, j, 2, s] - μ[i, j, s] == 0
        end)
        # Switch detection
        @constraints(
            model,
            begin
                [t=1:nh, j=1:nf], sum(θ[t, i, j, s] for i in 1:nc) == θs[t, j, s]
                [t=1:nh, j=1:nf], sum(λ[t, i, j, s] for i in 0:nc) == λs[t, j, s]
                [t=1:nh, i=1:nc, j=1:nf], [θ[t, i, j, s], λs[t, j, s]] ∈ MOI.Complements(2)
            end
        )
    end

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
        @expression(model, step_eq_cost, sum(1*(h .- (T/nh)) .^ 2))
    end
    @constraint(model, [j=1:N], sum(h[i, j] for i in 1:nfe) == T / N)

    @objective(model, Min, terminal_cost + step_eq_cost)
    return model
end

