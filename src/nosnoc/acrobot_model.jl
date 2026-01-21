# Acrobot (double pendulum) with Coulomb friction
#
# A two-link underactuated robot with viscous and Coulomb friction at both joints

"""
    nosnoc_acrobot_model

Optimal control of a double pendulum (acrobot) with Coulomb friction.

# System Description
A two-link pendulum with:
- Viscous friction at both joints (proportional to angular velocity)
- Coulomb friction at both joints (switches sign with velocity direction)
- Two independent switching functions (one per joint)

# Parameters
- N: number of control intervals
- nfe: number of finite elements per control interval
- rk: Runge-Kutta scheme
- big_M: big-M value for logical constraints
- step_eq: step equilibration method (:lcc or :heuristic_mean)

"""
function nosnoc_acrobot_model(N, nfe, rk::RKScheme; big_M=1e5, step_eq=:heuristic_mean, rho_h=1e2)
    nh = N * nfe      # total number of finite elements
    nf1 = 2           # number of modes for subsystem 1 (joint 1)
    nf2 = 2           # number of modes for subsystem 2 (joint 2)
    nc = length(rk.c) # total number of intermediate integration points

    T = 4.0
    hmin, hmax = 0.5*T/nh, 10.0*T/nh

    # System parameters
    m1 = 1.0      # mass of first link (kg)
    m2 = 1.2      # mass of second link (kg)
    l1 = 1.0      # length of first link (m)
    l2 = 1.2      # length of second link (m)
    g = 9.81      # gravity (m/s^2)
    b1 = 0.3      # viscous friction coefficient at joint 1
    b2 = 0.3      # viscous friction coefficient at joint 2
    mu1 = 0.7     # Coulomb friction coefficient at joint 1
    mu2 = 0.7     # Coulomb friction coefficient at joint 2

    # Control bounds
    tau_max = 10.0

    # Reference state and control
    x_ref = [π, 0.0, 0.0, 0.0]  # upright position
    u_ref = [0.0, 0.0]

    # Cost matrices
    Q = [10.0 0.0 0.0 0.0;
         0.0 10.0 0.0 0.0;
         0.0 0.0 1.0 0.0;
         0.0 0.0 0.0 1.0]
    R = [0.1 0.0;
         0.0 0.1]
    P = 10.0 * Q  # terminal cost

    # Initial condition
    x0 = [0.0, 0.0, 0.0, 0.0]

    model = Model()
    JuMP.set_name(model, "NOSNOC-Acrobot")

    # State variables: [q1, q2, q1_dot, q2_dot]
    @variable(model, -π <= q1[i=1:(nh+1), j=1:(nc+1)] <= π)
    @variable(model, -π <= q2[i=1:(nh+1), j=1:(nc+1)] <= π)
    @variable(model, q1_dot[i=1:(nh+1), j=1:(nc+1)])
    @variable(model, q2_dot[i=1:(nh+1), j=1:(nc+1)])

    # Control inputs: [tau1, tau2]
    @variable(model, -tau_max <= U[1:N, 1:2] <= tau_max)

    # Step in numerical time
    @variable(model, hmin <= h[1:nfe, 1:N] <= hmax, start=5*T/nh)

    # Stewart's multipliers for subsystem 1 (joint 1 friction)
    @variable(model, 0.0 <= θ1[1:nh, 1:nc, 1:nf1], start=0.5)
    @variable(model, 0.0 <= λ1[1:nh, 0:nc, 1:nf1], start=1.0)
    @variable(model, μ1[1:nh, 0:nc])

    # Stewart's multipliers for subsystem 2 (joint 2 friction)
    @variable(model, 0.0 <= θ2[1:nh, 1:nc, 1:nf2], start=0.5)
    @variable(model, 0.0 <= λ2[1:nh, 0:nc, 1:nf2], start=1.0)
    @variable(model, μ2[1:nh, 0:nc])

    # Switch detection for subsystem 1
    @variable(model, 0.0 <= λs1[1:nh, 1:nf1], start=1.0)
    @variable(model, 0.0 <= θs1[1:nh, 1:nf1], start=0.5)

    # Switch detection for subsystem 2
    @variable(model, 0.0 <= λs2[1:nh, 1:nf2], start=1.0)
    @variable(model, 0.0 <= θs2[1:nh, 1:nf2], start=0.5)

    if step_eq == :lcc
        @variable(model, 0.0 <= πλ1[1:(nh-1), 1:nf1])
        @variable(model, 0.0 <= πθ1[1:(nh-1), 1:nf1])
        @variable(model, 0.0 <= τ1_1[1:(nh-1), 1:nf1])
        @variable(model, 0.0 <= τ1_2[1:(nh-1), 1:nf1])
        @variable(model, 0.0 <= vcc1[1:(nh-1), 1:nf1])

        @variable(model, 0.0 <= πλ2[1:(nh-1), 1:nf2])
        @variable(model, 0.0 <= πθ2[1:(nh-1), 1:nf2])
        @variable(model, 0.0 <= τ2_1[1:(nh-1), 1:nf2])
        @variable(model, 0.0 <= τ2_2[1:(nh-1), 1:nf2])
        @variable(model, 0.0 <= vcc2[1:(nh-1), 1:nf2])
    end

    # Cost function
    @expression(
        model,
        stage_cost,
        1/N * sum(
            (q1[i*nfe+1, 1] - x_ref[1])^2 * Q[1,1] +
            (q2[i*nfe+1, 1] - x_ref[2])^2 * Q[2,2] +
            (q1_dot[i*nfe+1, 1] - x_ref[3])^2 * Q[3,3] +
            (q2_dot[i*nfe+1, 1] - x_ref[4])^2 * Q[4,4] +
            R[1,1] * (U[i, 1] - u_ref[1])^2 +
            R[2,2] * (U[i, 2] - u_ref[2])^2
            for i in 1:N
        )
    )

    @expression(
        model,
        terminal_cost,
        (q1[nh+1, 1] - x_ref[1])^2 * P[1,1] +
        (q2[nh+1, 1] - x_ref[2])^2 * P[2,2] +
        (q1_dot[nh+1, 1] - x_ref[3])^2 * P[3,3] +
        (q2_dot[nh+1, 1] - x_ref[4])^2 * P[4,4]
    )

    # Initial conditions
    @constraints(model, begin
        q1[1, 1] == x0[1]
        q2[1, 1] == x0[2]
        q1_dot[1, 1] == x0[3]
        q2_dot[1, 1] == x0[4]
    end)

    # Stewart position functions (switching based on joint velocities)
    @expression(
        model,
        g1[i=1:nh, j=1:(nc+1)],
        q1_dot[i, j]
    )

    @expression(
        model,
        g2[i=1:nh, j=1:(nc+1)],
        q2_dot[i, j]
    )

    # Initial position for multipliers
    @constraint(model, λ1[1, 0, 1] ==  g1[1, 1] - μ1[1, 0])  # q1_dot > 0
    @constraint(model, λ1[1, 0, 2] == -g1[1, 1] - μ1[1, 0])  # q1_dot < 0

    @constraint(model, λ2[1, 0, 1] ==  g2[1, 1] - μ2[1, 0])  # q2_dot > 0
    @constraint(model, λ2[1, 0, 2] == -g2[1, 1] - μ2[1, 0])  # q2_dot < 0

    # Mass matrix, Coriolis, gravity, viscous friction expressions
    @expressions(
        model,
        begin
            M11[i=1:nh, j=1:nc], m1*l1^2/3 + m2*(l1^2 + l2^2/3 + l1*l2*cos(q2[i, j]))
            M12[i=1:nh, j=1:nc], m2*(l2^2/3 + l1*l2*cos(q2[i, j])/2)
            M22[i=1:nh, j=1:nc], m2*l2^2/3

            C1[i=1:nh, j=1:nc], -m2*l1*l2*sin(q2[i, j])*(2*q1_dot[i, j]*q2_dot[i, j] + q2_dot[i, j]^2)/2
            C2[i=1:nh, j=1:nc], m2*l1*l2*sin(q2[i, j])*(q1_dot[i, j]^2)/2

            G1[i=1:nh, j=1:nc], -g*(m1*l1/2 + m2*l1)*sin(q1[i, j]) - g*m2*l2/2*sin(q1[i, j] + q2[i, j])
            G2[i=1:nh, j=1:nc], -g*m2*l2/2*sin(q1[i, j] + q2[i, j])

            F_v1[i=1:nh, j=1:nc], -b1*q1_dot[i, j]
            F_v2[i=1:nh, j=1:nc], -b2*q2_dot[i, j]
        end
    )

    # Compute determinant for matrix inversion
    @expression(model, det_M[i=1:nh, j=1:nc], M11[i, j]*M22[i, j] - M12[i, j]^2)

    # Control input indexing
    @expression(model, tau1[i=1:nh], U[div(i-1, nfe)+1, 1])
    @expression(model, tau2[i=1:nh], U[div(i-1, nfe)+1, 2])

    # Base dynamics (without Coulomb friction)
    @expressions(
        model,
        begin
            # Accelerations from M \ (tau - C - G - F_v)
            q1_ddot_base[i=1:nh, j=1:nc], (M22[i, j]*(tau1[i] - C1[i, j] - G1[i, j] - F_v1[i, j]) -
                                            M12[i, j]*(tau2[i] - C2[i, j] - G2[i, j] - F_v2[i, j])) / det_M[i, j]
            q2_ddot_base[i=1:nh, j=1:nc], (M11[i, j]*(tau2[i] - C2[i, j] - G2[i, j] - F_v2[i, j]) -
                                            M12[i, j]*(tau1[i] - C1[i, j] - G1[i, j] - F_v1[i, j])) / det_M[i, j]
        end
    )

    # Coulomb friction contributions
    @expressions(
        model,
        begin
            # Joint 1 Coulomb friction modes
            f1_coulomb_mode1[i=1:nh, j=1:nc], -mu1  # q1_dot > 0
            f1_coulomb_mode2[i=1:nh, j=1:nc],  mu1  # q1_dot < 0

            # Joint 2 Coulomb friction modes
            f2_coulomb_mode1[i=1:nh, j=1:nc], -mu2  # q2_dot > 0
            f2_coulomb_mode2[i=1:nh, j=1:nc],  mu2  # q2_dot < 0
        end
    )

    # Transform Coulomb friction to acceleration space using M^(-1)
    @expressions(
        model,
        begin
            # Effect of joint 1 Coulomb friction on accelerations
            q1_ddot_coulomb1_mode1[i=1:nh, j=1:nc], M22[i, j] * f1_coulomb_mode1[i, j] / det_M[i, j]
            q1_ddot_coulomb1_mode2[i=1:nh, j=1:nc], M22[i, j] * f1_coulomb_mode2[i, j] / det_M[i, j]
            q2_ddot_coulomb1_mode1[i=1:nh, j=1:nc], -M12[i, j] * f1_coulomb_mode1[i, j] / det_M[i, j]
            q2_ddot_coulomb1_mode2[i=1:nh, j=1:nc], -M12[i, j] * f1_coulomb_mode2[i, j] / det_M[i, j]

            # Effect of joint 2 Coulomb friction on accelerations
            q1_ddot_coulomb2_mode1[i=1:nh, j=1:nc], -M12[i, j] * f2_coulomb_mode1[i, j] / det_M[i, j]
            q1_ddot_coulomb2_mode2[i=1:nh, j=1:nc], -M12[i, j] * f2_coulomb_mode2[i, j] / det_M[i, j]
            q2_ddot_coulomb2_mode1[i=1:nh, j=1:nc], M11[i, j] * f2_coulomb_mode1[i, j] / det_M[i, j]
            q2_ddot_coulomb2_mode2[i=1:nh, j=1:nc], M11[i, j] * f2_coulomb_mode2[i, j] / det_M[i, j]
        end
    )

    # Combined dynamics with both subsystems
    @expressions(
        model,
        begin
            dq1[i=1:nh, j=1:nc], q1_dot[i, j]
            dq2[i=1:nh, j=1:nc], q2_dot[i, j]
            dq1_dot[i=1:nh, j=1:nc], q1_ddot_base[i, j] +
                                      θ1[i, j, 1] * q1_ddot_coulomb1_mode1[i, j] +
                                      θ1[i, j, 2] * q1_ddot_coulomb1_mode2[i, j] +
                                      θ2[i, j, 1] * q1_ddot_coulomb2_mode1[i, j] +
                                      θ2[i, j, 2] * q1_ddot_coulomb2_mode2[i, j]
            dq2_dot[i=1:nh, j=1:nc], q2_ddot_base[i, j] +
                                      θ1[i, j, 1] * q2_ddot_coulomb1_mode1[i, j] +
                                      θ1[i, j, 2] * q2_ddot_coulomb1_mode2[i, j] +
                                      θ2[i, j, 1] * q2_ddot_coulomb2_mode1[i, j] +
                                      θ2[i, j, 2] * q2_ddot_coulomb2_mode2[i, j]
        end
    )

    # Collocations
    @constraints(
        model,
        begin
            con_dq1[i=1:nh, j=2:(nc+1)],
            q1[i, j] == q1[i, 1] + h[i] * sum(rk.a[j-1, k] * dq1[i, k] for k in 1:nc)
            con_dq2[i=1:nh, j=2:(nc+1)],
            q2[i, j] == q2[i, 1] + h[i] * sum(rk.a[j-1, k] * dq2[i, k] for k in 1:nc)
            con_dq1_dot[i=1:nh, j=2:(nc+1)],
            q1_dot[i, j] == q1_dot[i, 1] + h[i] * sum(rk.a[j-1, k] * dq1_dot[i, k] for k in 1:nc)
            con_dq2_dot[i=1:nh, j=2:(nc+1)],
            q2_dot[i, j] == q2_dot[i, 1] + h[i] * sum(rk.a[j-1, k] * dq2_dot[i, k] for k in 1:nc)
        end
    )

    # Continuity w.r.t primal variables
    @constraints(
        model,
        begin
            [i=1:nh], q1[i+1, 1] == q1[i, 1] + h[i] * sum(rk.b[k] * dq1[i, k] for k in 1:nc)
            [i=1:nh], q2[i+1, 1] == q2[i, 1] + h[i] * sum(rk.b[k] * dq2[i, k] for k in 1:nc)
            [i=1:nh], q1_dot[i+1, 1] == q1_dot[i, 1] + h[i] * sum(rk.b[k] * dq1_dot[i, k] for k in 1:nc)
            [i=1:nh], q2_dot[i+1, 1] == q2_dot[i, 1] + h[i] * sum(rk.b[k] * dq2_dot[i, k] for k in 1:nc)
        end
    )

    # Continuity w.r.t dual variables
    @constraints(model, begin
        [i=1:(nh-1), j=1:nf1], λ1[i, nc, j] == λ1[i+1, 0, j]
        [i=1:(nh-1), j=1:nf2], λ2[i, nc, j] == λ2[i+1, 0, j]
    end)

    # Stewart model for subsystem 1
    @constraints(model, begin
        [i=1:nh, j=1:nc], sum(θ1[i, j, k] for k in 1:nf1) == 1.0
        [i=1:nh, j=0:nc],  g1[i, j+1] - λ1[i, j, 1] - μ1[i, j] == 0
        [i=1:nh, j=0:nc], -g1[i, j+1] - λ1[i, j, 2] - μ1[i, j] == 0
    end)

    # Stewart model for subsystem 2
    @constraints(model, begin
        [i=1:nh, j=1:nc], sum(θ2[i, j, k] for k in 1:nf2) == 1.0
        [i=1:nh, j=0:nc],  g2[i, j+1] - λ2[i, j, 1] - μ2[i, j] == 0
        [i=1:nh, j=0:nc], -g2[i, j+1] - λ2[i, j, 2] - μ2[i, j] == 0
    end)

    # Switch detection for subsystem 1
    @constraints(
        model,
        begin
            [t=1:nh, j=1:nf1], sum(θ1[t, i, j] for i in 1:nc) == θs1[t, j]
            [t=1:nh, j=1:nf1], sum(λ1[t, i, j] for i in 0:nc) == λs1[t, j]
            [t=1:nh, i=1:nc, j=1:nf1], [θ1[t, i, j], λs1[t, j]] ∈ MOI.Complements(2)
        end
    )

    # Switch detection for subsystem 2
    @constraints(
        model,
        begin
            [t=1:nh, j=1:nf2], sum(θ2[t, i, j] for i in 1:nc) == θs2[t, j]
            [t=1:nh, j=1:nf2], sum(λ2[t, i, j] for i in 0:nc) == λs2[t, j]
            [t=1:nh, i=1:nc, j=1:nf2], [θ2[t, i, j], λs2[t, j]] ∈ MOI.Complements(2)
        end
    )

    if step_eq == :lcc
        # Logical constraints for subsystem 1
        @constraints(
            model,
            begin
                [t=1:(nh-1), j=1:nf1], πλ1[t, j] >= λs1[t, j]
                [t=1:(nh-1), j=1:nf1], πλ1[t, j] >= λs1[t+1, j]
                [t=1:(nh-1), j=1:nf1], πλ1[t, j] <= λs1[t, j] + λs1[t+1, j]
                [t=1:(nh-1), j=1:nf1], πθ1[t, j] >= θs1[t, j]
                [t=1:(nh-1), j=1:nf1], πθ1[t, j] >= θs1[t+1, j]
                [t=1:(nh-1), j=1:nf1], πθ1[t, j] <= θs1[t, j] + θs1[t+1, j]
                [t=1:(nh-1), j=1:nf1], vcc1[t, j] <= πλ1[t, j]
                [t=1:(nh-1), j=1:nf1], vcc1[t, j] <= πθ1[t, j]
                [t=1:(nh-1), j=1:nf1], vcc1[t, j] >= πθ1[t, j] - τ1_1[t, j]
                [t=1:(nh-1), j=1:nf1], πθ1[t, j] - πλ1[t, j] == τ1_1[t, j] - τ1_2[t, j]
                [t=1:(nh-1), j=1:nf1], [τ1_1[t, j], τ1_2[t, j]] ∈ MOI.Complements(2)
            end
        )

        # Logical constraints for subsystem 2
        @constraints(
            model,
            begin
                [t=1:(nh-1), j=1:nf2], πλ2[t, j] >= λs2[t, j]
                [t=1:(nh-1), j=1:nf2], πλ2[t, j] >= λs2[t+1, j]
                [t=1:(nh-1), j=1:nf2], πλ2[t, j] <= λs2[t, j] + λs2[t+1, j]
                [t=1:(nh-1), j=1:nf2], πθ2[t, j] >= θs2[t, j]
                [t=1:(nh-1), j=1:nf2], πθ2[t, j] >= θs2[t+1, j]
                [t=1:(nh-1), j=1:nf2], πθ2[t, j] <= θs2[t, j] + θs2[t+1, j]
                [t=1:(nh-1), j=1:nf2], vcc2[t, j] <= πλ2[t, j]
                [t=1:(nh-1), j=1:nf2], vcc2[t, j] <= πθ2[t, j]
                [t=1:(nh-1), j=1:nf2], vcc2[t, j] >= πθ2[t, j] - τ2_1[t, j]
                [t=1:(nh-1), j=1:nf2], πθ2[t, j] - πλ2[t, j] == τ2_1[t, j] - τ2_2[t, j]
                [t=1:(nh-1), j=1:nf2], [τ2_1[t, j], τ2_2[t, j]] ∈ MOI.Complements(2)
            end
        )

        # Step equilibration
        @constraints(
            model,
            begin
                [t=1:(nh-1)], -big_M * (sum(vcc1[t, j] for j in 1:nf1) + sum(vcc2[t, j] for j in 1:nf2)) <= (h[t+1] - h[t])
                [t=1:(nh-1)], (h[t+1] - h[t]) <= big_M * (sum(vcc1[t, j] for j in 1:nf1) + sum(vcc2[t, j] for j in 1:nf2))
            end
        )
        @expression(model, step_eq_cost, 0)
    elseif step_eq == :heuristic_mean
        @expression(model, step_eq_cost, rho_h * sum(1*(h .- (T/nh)) .^ 2))
    end

    @constraint(model, [j=1:N], sum(h[i, j] for i in 1:nfe) == T / N)

    @objective(model, Min, stage_cost + terminal_cost + step_eq_cost)

    return model
end
