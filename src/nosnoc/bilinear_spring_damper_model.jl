# Bilinear spring-damper system with state and input constraints
# 
# A mass-spring-damper system with different spring constants 
# in positive and negative displacement regions

"""
    nosnoc_bilinear_spring_damper_model

Optimal control of a bilinear spring-damper system.

# System Description
A mass-spring-damper system where the spring constant switches 
between two values depending on the sign of the position:
- k1 when x < 0 (negative region)
- k2 when x > 0 (positive region)

# Parameters
- N: number of control intervals
- nfe: number of finite elements per control interval
- rk: Runge-Kutta scheme
- big_M: big-M value for logical constraints
- step_eq: step equilibration method (:lcc or :heuristic_mean)

"""
function nosnoc_bilinear_spring_damper_model(N, nfe, rk::RKScheme; big_M=1e5, step_eq=:heuristic_mean)
    nh = N * nfe      # total number of finite elements
    nf = 2            # total number of nonsmooth modes
    nc = length(rk.c) # total number of intermediate integration points

    T = 5.0
    hmin, hmax = 0.5*T/nh, 10.0*T/nh

    # System parameters
    m = 1.0       # mass (kg)
    k1 = 0.5      # spring constant in negative region (N/m)
    k2 = 3.0      # spring constant in positive region (N/m)
    c = 0.7       # damping coefficient (Ns/m)
    
    # Constraint bounds
    x1_max = 3.0  # position bound
    x2_max = 4.0  # velocity bound
    u_max = 5.0   # control input bound
    
    # Reference state and control
    x_ref = [1.5, 0.0]  # reference state [position; velocity]
    u_ss2 = -k2/m*x_ref[1] - c/m*x_ref[2]
    u_ref = -u_ss2
    
    # Cost matrices
    Q = [1.0 0.0; 0.0 0.5]      # state tracking weight
    R = 0.1                      # control input weight
    P = [5.0 0.0; 0.0 2.5]      # terminal cost weight
    
    # Initial condition
    x0 = [-1.0, 0.0]

    model = Model()
    JuMP.set_name(model, "NOSNOC-BilinearSpringDamper")
    
    # State variables: [position, velocity]
    @variable(model, -x1_max <= x1[i=1:(nh+1), j=1:(nc+1)] <= x1_max)
    @variable(model, -x2_max <= x2[i=1:(nh+1), j=1:(nc+1)] <= x2_max)
    
    # Control input
    @variable(model, -u_max <= U[1:N] <= u_max)
    
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

    # Cost function
    @expression(
        model,
        stage_cost,
        1/N * sum(
            (x1[i*nfe+1, 1] - x_ref[1])^2 * Q[1,1] + 
            (x2[i*nfe+1, 1] - x_ref[2])^2 * Q[2,2] + 
            R * (U[i] - u_ref)^2 
            for i in 1:N
        )
    )
    
    @expression(
        model,
        terminal_cost,
        (x1[nh+1, 1] - x_ref[1])^2 * P[1,1] + 
        (x2[nh+1, 1] - x_ref[2])^2 * P[2,2]
    )

    # Initial conditions
    @constraints(model, begin
        x1[1, 1] == x0[1]
        x2[1, 1] == x0[2]
    end)
    
    # Stewart position function (switching based on position)
    @expression(
        model,
        g[i=1:nh, j=1:(nc+1)],
        x1[i, j]
    )
    
    # Initial position for multipliers
    @constraint(model, λ[1, 0, 1] == -g[1, 1] - μ[1, 0])  # c < 0 region
    @constraint(model, λ[1, 0, 2] ==  g[1, 1] - μ[1, 0])  # c > 0 region

    # Dynamics
    @expressions(
        model,
        begin
            # Mode 1: x < 0 (negative region, spring constant k1)
            dx1_mode1[i=1:nh, j=1:nc], x2[i, j]
            dx2_mode1[i=1:nh, j=1:nc], -k1/m*x1[i, j] - c/m*x2[i, j] + 1/m*U[div(i-1, nfe)+1]
            
            # Mode 2: x > 0 (positive region, spring constant k2)
            dx1_mode2[i=1:nh, j=1:nc], x2[i, j]
            dx2_mode2[i=1:nh, j=1:nc], -k2/m*x1[i, j] - c/m*x2[i, j] + 1/m*U[div(i-1, nfe)+1]
            
            # Combined dynamics using Stewart multipliers
            dx1[i=1:nh, j=1:nc], θ[i, j, 1] * dx1_mode1[i, j] + θ[i, j, 2] * dx1_mode2[i, j]
            dx2[i=1:nh, j=1:nc], θ[i, j, 1] * dx2_mode1[i, j] + θ[i, j, 2] * dx2_mode2[i, j]
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
        end
    )
    
    # Continuity w.r.t primal variables
    @constraints(
        model,
        begin
            [i=1:nh], x1[i+1, 1] == x1[i, 1] + h[i] * sum(rk.b[k] * dx1[i, k] for k in 1:nc)
            [i=1:nh], x2[i+1, 1] == x2[i, 1] + h[i] * sum(rk.b[k] * dx2[i, k] for k in 1:nc)
        end
    )
    
    # Continuity w.r.t dual variables
    @constraints(model, begin
        [i=1:(nh-1), j=1:nf], λ[i, nc, j] == λ[i+1, 0, j]
    end)
    
    # Stewart model
    @constraints(model, begin
        [i=1:nh, j=1:nc], sum(θ[i, j, k] for k in 1:nf) == 1.0
        [i=1:nh, j=0:nc], -g[i, j+1] - λ[i, j, 1] - μ[i, j] == 0
        [i=1:nh, j=0:nc],  g[i, j+1] - λ[i, j, 2] - μ[i, j] == 0
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
        @expression(model, step_eq_cost, sum(1*(h .- (T/nh)) .^ 2))
    end
    
    @constraint(model, [j=1:N], sum(h[i, j] for i in 1:nfe) == T / N)

    @objective(model, Min, stage_cost + terminal_cost + step_eq_cost)
    
    return model
end