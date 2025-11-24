# Adapted from NOSBENCH
# https://github.com/nosnoc/nosbench
#
# Nurkanović, A., Pozharskiy, A. & Diehl, M.
# Solving Mathematical Programs with Complementarity Constraints Arising in Nonsmooth Optimal Control.
# Vietnam J. Math. 53, 659–697 (2025).
# https://doi.org/10.1007/s10013-024-00704-z

"""
    nosnoc_cart_pole_with_friction_model

Inverted pendulum on a cart with Coulomb friction.

# Reference

Howell, Taylor A., et al.
"Trajectory optimization with optimization-based dynamics."
IEEE Robotics and Automation Letters 7.3 (2022): 6750-6757.

"""
function nosnoc_cart_pole_with_friction_model(N, nfe, rk::RKScheme; big_M=1e5, rho_h=1e0, step_eq=:heuristic_mean)
    nh = N * nfe      # total number of finite elements
    nf = 2            # total number of nonsmooth modes
    nc = length(rk.c) # total number of intermediate integration points

    T = 4.0
    hmin, hmax = 0.1/nh, 2.0*T/nh

    m1 = 1.0
    m2 = 0.1
    ll = 1.0 # link length
    gravity = 9.81
    μ_fric = 2.0

    q0 = [1.0, 0.0]
    v0 = [0.0, 0.0]
    x_ref = [0.0, pi, 0.0, 0.0]
    Q = diagm([1.0, 100.0, 1.0, 1.0])
    Q_terminal = diagm([100.0, 100.0, 10.0, 10.0])
    R = 1.0
    q_max = [5.0, 240.0/180.0*pi]
    q_min = [0.0, -240.0/180.0*pi]
    v_max = 20.0
    u_max = 30.0
    u_ref = 0.0

    q_guess = [q0 .+ i/(nh+1) .* (x_ref[1:2] .- q0) for i in 1:nh+1]

    model = Model()
    JuMP.set_name(model, "NOSNOC-CPWF")
    # State
    @variable(model, q_min[k] <= q[i=1:(nh+1), j=1:(nc+1), k=1:2] <= q_max[k], start=q_guess[i][k])
    @variable(model, -v_max <= v[i=1:(nh+1), j=1:(nc+1), k=1:2] <= v_max, start=0.0)
    @variable(model, -u_max <= u[1:N] <= u_max, start=1.0)
    # Step in numerical time
    @variable(model, hmin <= h[1:nfe, 1:N] <= hmax, start=T/nh)
    # Stewart's multipliers
    @variable(model, 0.0 <= θ[1:nh, 1:nc, 1:nf], start=0.5)
    @variable(model, 0.0 <= λ[1:nh, 0:nc, 1:nf], start=1.0)
    @variable(model, μ[1:nh, 0:nc])
    # Switch detection
    @variable(model, 0.0 <= λs[1:nh, 1:nf], start=0.5)
    @variable(model, 0.0 <= θs[1:nh, 1:nf], start=1.0)
    if step_eq == :lcc
        @variable(model, 0.0 <= πλ[1:(nh-1), 1:nf])
        @variable(model, 0.0 <= πθ[1:(nh-1), 1:nf])
        @variable(model, 0.0 <= τ1[1:(nh-1), 1:nf])
        @variable(model, 0.0 <= τ2[1:(nh-1), 1:nf])
        @variable(model, 0.0 <= vcc[1:(nh-1), 1:nf])
    end

    @expression(model, x[i=1:(N+1)], [q[(i-1)*nfe+1, 1, 1], q[(i-1)*nfe+1, 1, 2], v[(i-1)*nfe+1, 1, 1], v[(i-1)*nfe+1, 1, 2]])

    @expression(
        model,
        terminal_cost,
        1/N*sum(dot(x[i] .- x_ref, Q, x[i] .- x_ref) + R * (u[i] - u_ref)^2 for i in 1:N)
        + dot(x[N+1] .- x_ref, Q_terminal, x[N+1] .- x_ref)
    )

    # Initial position
    @constraints(model, begin
        q0, q[1, 1, :] .== q0
        v0, v[1, 1, :] .== v0
    end)
    # Stewart position function
    @expression(
        model,
        g[i=1:nh, j=1:(nc+1)],
        v[i, j, 1]
    )
    # Initial position for multipliers
    @constraint(model, λ[1, 0, 1] ==  g[1, 1] - μ[1, 0])
    @constraint(model, λ[1, 0, 2] == -g[1, 1] - μ[1, 0])
    # Dynamics
    @expression(
        model,
        invM[i=1:nh, j=1:nc],
        1 / ((m1+m2)*m2*ll^2 - m2^2*ll^2*cos(q[i, j, 2])) .* [m2 * ll^2    -m2*ll*cos(q[i, j, 2]);
                                                            -m2*ll*cos(q[i, j, 2])   m1+m2]
    )
    @expression(
        model,
        force[i=1:nh, j=1:nc],
        [
            u[div(i-1, nfe)+1] - m2*ll*v[i, j, 2]^2*sin(q[i, j, 2]) - μ_fric * (θ[i, j, 1] - θ[i, j, 2]) ;
            -m2*gravity*ll*sin(q[i, j, 2])
        ]
    )
    @expressions(
        model,
        begin
            dq[i=1:nh, j=1:nc, l=1:2], v[i, j, l]
            dv[i=1:nh, j=1:nc, l=1:2], invM[i, j][l, 1] * force[i, j][1] + invM[i, j][l, 2] * force[i, j][2]
        end
    )
    # Collocations
    @constraints(
        model,
        begin
            con_dq[i=1:nh, j=2:(nc+1), l=1:2],
            q[i, j, l] == q[i, 1, l] + h[i] * sum(rk.a[j-1, k] * dq[i, k, l] for k in 1:nc)
            con_dv[i=1:nh, j=2:(nc+1), l=1:2],
            v[i, j, l] == v[i, 1, l] + h[i] * sum(rk.a[j-1, k] * dv[i, k, l] for k in 1:nc)
        end
    )
    # Continuity
    # w.r.t primal variables
    @constraints(
        model,
        begin
            [i=1:nh, l=1:2],
            q[i+1, 1, l] == q[i, 1, l] + h[i] * sum(rk.b[k] * dq[i, k, l] for k in 1:nc)
            [i=1:nh, l=1:2],
            v[i+1, 1, l] == v[i, 1, l] + h[i] * sum(rk.b[k] * dv[i, k, l] for k in 1:nc)
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

    @objective(model, Min, terminal_cost + step_eq_cost)
    return model
end

