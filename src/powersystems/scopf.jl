
"""
    scopf_model(data, contingencies)

Implement corrective AC-SCOPF with complementarity constraints from [1].

# Reference

[1] F Pacaud, A Nurkanović, A Pozharskiy, A Montoison, S Shin.
"An Augmented Lagrangian Method on GPU for Security-Constrained AC Optimal Power Flow", 2025.

"""
function scopf_model(
    data,
    contingencies;
    scale_obj=1e-4,
    load_factor=1.0,
)
    nbus = length(data.bus)
    ngen = length(data.gen)
    nbranch = length(data.branch)
    narc = length(data.arc)
    alpha = ones(ngen)

    pv_buses = [k for k in 1:nbus if data.bus[k].type == 2]
    ref_buses = [k for k in 1:nbus if data.bus[k].type == 3]
    npv = length(pv_buses)

    # Pre-processing
    connected_arcs = [Int[] for k in 1:nbus]
    for (k, arc) in enumerate(data.arc)
        b = arc.bus
        push!(connected_arcs[b], k)
    end
    connected_gens = [Int[] for k in 1:nbus]
    qmin = zeros(nbus)
    qmax = zeros(nbus)
    qg0 = zeros(nbus)
    for (k, gen) in enumerate(data.gen)
        b = gen.bus
        push!(connected_gens[b], k)
        # Aggregate reactive power generations
        qmin[b] += gen.qmin
        qmax[b] += gen.qmax
        qg0[b] += gen.qg
    end

    # Total number of scenarios
    K = length(contingencies) + 1

    # Build model
    model = JuMP.Model()
    JuMP.set_name(model, "Corrective-AC-SCOPF")

    JuMP.@variable(model, va[i in 1:nbus, 1:K], start=deg2rad(data.bus[i].va))
    JuMP.@variable(model, data.bus[i].vmin <= vm[i in 1:nbus, 1:K] <= data.bus[i].vmax, start=data.bus[i].vm)
    JuMP.@variable(model, data.gen[i].pmin <= pg[i in 1:ngen, 1:K] <= data.gen[i].pmax, start=data.gen[i].pg)
    # N.B. Aggregate the reactive power generations for the PV/PQ switches
    JuMP.@variable(model, qmin[i] <= qg[i in [ref_buses; pv_buses], 1:K] <= qmax[i], start=qg0[i])
    JuMP.@variable(model, -data.arc[i].rate_a <= p[i in 1:narc, 1:K] <= data.arc[i].rate_a)
    JuMP.@variable(model, -data.arc[i].rate_a <= q[i in 1:narc, 1:K] <= data.arc[i].rate_a)
    # Automatic adjustment of generators
    JuMP.@variable(model, Δ[1:K-1])
    # Switches
    JuMP.@variable(model, γ[i in 1:ngen, 2:K])
    JuMP.@variable(model, λ[i in [ref_buses; pv_buses], 2:K])

    JuMP.@objective(model, Min, scale_obj * sum(data.gen[i].c[1]*pg[i, 1]^2 + data.gen[i].c[2]*pg[i, 1] + data.gen[i].c[3] for i in 1:ngen))

    for k in 2:K
        # Droop control
        for i in 1:ngen
            @constraint(model, pg[i, k] == pg[i, 1] + alpha[i] * Δ[k-1] + γ[i, k])
            @constraint(model, [γ[i, k], pg[i, k]] in MOI.Complements(2))
        end
        # PV/PQ switches
        # N.B. : we allow the reference buses to switch as well
        for b in [ref_buses; pv_buses]
            @constraint(model, vm[b, k] == vm[b, 1] + λ[b, k])
            @constraint(model, [λ[b, k], qg[b, k]] in MOI.Complements(2))
        end
    end

    for k in 1:K
        for i in ref_buses
            JuMP.@constraint(model, va[i, k] == 0)
        end

        for b in 1:nbus
            JuMP.@constraint(model,
                sum(p[a, k] for a in connected_arcs[b]) ==
                sum(pg[g, k] for g in connected_gens[b])
                - load_factor * data.bus[b].pd
                - data.bus[b].gs * vm[b, k]^2
            )
            if data.bus[b].type == 1
                JuMP.@constraint(model,
                    sum(q[a, k] for a in connected_arcs[b]) ==
                    - load_factor * data.bus[b].qd
                    + data.bus[b].bs * vm[b, k]^2
                )
            else
                JuMP.@constraint(model,
                    sum(q[a, k] for a in connected_arcs[b]) ==
                    qg[b, k]
                    - load_factor * data.bus[b].qd
                    + data.bus[b].bs * vm[b, k]^2
                )
            end
        end

        # Branch power flow physics and limit constraints
        for l in data.branch
            f_idx = l.f_idx
            t_idx = l.t_idx

            p_fr = p[f_idx, k]
            q_fr = q[f_idx, k]
            p_to = p[t_idx, k]
            q_to = q[t_idx, k]

            if (k >= 2) && l.i == contingencies[k-1]
                # Set flux to 0
                @constraint(model, p_fr == 0.0)
                @constraint(model, q_fr == 0.0)
                @constraint(model, p_to == 0.0)
                @constraint(model, q_to == 0.0)
            else

                vm_fr = vm[l.f_bus, k]
                vm_to = vm[l.t_bus, k]
                va_fr = va[l.f_bus, k]
                va_to = va[l.t_bus, k]

                # FR
                JuMP.@constraint(
                    model,
                    p_fr - l.c5*vm_fr^2
                    - l.c3 * (vm_fr * vm_to * cos(va_fr - va_to)) - l.c4 * (vm_fr * vm_to * sin(va_fr - va_to))
                    == 0.0
                )
                JuMP.@constraint(
                    model,
                    q_fr + l.c6*vm_fr^2
                    + l.c4 * (vm_fr * vm_to * cos(va_fr - va_to))
                    - l.c3 * (vm_fr * vm_to * sin(va_fr - va_to))
                    == 0.0
                )

                # TO
                JuMP.@constraint(
                    model,
                    p_to - l.c7*vm_to^2
                    - l.c1 * (vm_to * vm_fr * cos(va_to - va_fr))
                    - l.c2 * (vm_to * vm_fr * sin(va_to - va_fr))
                    == 0.0
                )
                JuMP.@constraint(
                    model,
                    q_to + l.c8*vm_to^2
                    + l.c2 * (vm_to * vm_fr * cos(va_to - va_fr))
                    - l.c1 * (vm_to * vm_fr * sin(va_to - va_fr))
                    == 0.0
                )

                # Apparent power limit, from side and to side
                JuMP.@constraint(model, p_fr^2 + q_fr^2 <= l.rate_a^2)
                JuMP.@constraint(model, p_to^2 + q_to^2 <= l.rate_a^2)
            end
        end
    end

    return model
end
